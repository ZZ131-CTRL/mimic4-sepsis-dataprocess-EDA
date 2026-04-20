###################################################
# 脓毒症患者你CT-HMM数据EDA分析
# DATE:20260325
# 1.基础设置和function定义
###################################################
{
  rm(list = objects())
  library(dplyr)     # 数据处理与操作，如过滤、分组、汇总、排序等（管道操作非常常用 %>%）
  library(ggplot2)   # 数据可视化，构建灵活、层次化的图形（如折线图、散点图、柱状图等）
  library(ggpubr)    # 美化 ggplot 图形，添加显著性标注、排列图形、生成出版级图表
  library(GGally)    # 扩展 ggplot2 功能，常用于画成对关系图（如 ggpairs() 多变量可视化）
  library(caTools)   # 包含各种实用工具，如数据拆分（split）、移动平均、AUC 计算等
  library(stringr)   # 字符串处理，如提取、替换、模式匹配等（基于 stringi，更直观易用）
  library(xgxr)      # 辅助 PK/PD 数据建模可视化
  library(ggbreak)   # 为 ggplot 添加断轴（break axis）功能，适用于跨度大或极端值的可视化
  library(tidyr)     # 数据整理（如列拆分 separate()，宽长格式转换等）
  library(patchwork) # 组合多个 ggplot 图形为一个图，语法直观简洁
  library(ggforce)   # 分页
  library(ggplot2)
  library(cowplot)   # 用于ggplot2图形的增强和组合。
  library(openxlsx)
  library(corrplot)
  library(tibble)
  library(readxl)
  library(purrr)
  library(lubridate)
  library(lme4)      # 用于计算整体显著性
  library(lmerTest)
  
}
# 一、数据处理 #----------------------------------------------------------------
### 1.基础设置和function定义 ####-----------------------------------------------
{
  Totalworkpath <- "D:\\Lab_project\\2026work\\sepsis"
  raw_data <- "DATA\\sepsis\\original_data"
  derived_data <- "DATA\\sepsis\\derived_data"
  density <- "PLOT\\sepsis\\density"
  CT_curve  <- "PLOT\\sepsis\\CT_curve"
  correlation  <- "PLOT\\sepsis\\correlation"
  
  Input_derived_data <- file.path(Totalworkpath,raw_data)
  Output_derived_data <- file.path(Totalworkpath,derived_data)
  Output_density <- file.path(Totalworkpath,density)           # 协变量信息输出路径
  Output_CTcurve <- file.path(Totalworkpath,CT_curve)  # CT曲线输出路径
  Output_correlation <- file.path(Totalworkpath,correlation)
  
  # 用于绘制相关性总图
  lowerFn <- function(data, mapping, method = "lm", ...) {
    p <- ggplot(data = data, mapping = mapping) +
      geom_point(colour = "blue",alpha = 0.3,size = 0.1) +
      geom_smooth(method = method, color = "red", size = 0.5,...)
    p
  }
  options(scipen = 999)  # 不使用科学计数法
}
{
  theme_cor <-theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    panel.border = element_rect(colour = "black", fill=NA, linewidth =1.5),
                    axis.title = element_text(size = 14, face = "bold"),
                    axis.title.x = element_text(vjust = -0.8),
                    axis.title.y = element_text(vjust = 1.5),
                    axis.text = element_text(size = 12, face = "bold", color = "black"),
                    axis.ticks = element_line(color = "black", linewidth = 1),
                    axis.ticks.length = unit(0.15, "cm"),
                    legend.position = "right",  # 放置图例到右侧
                    legend.justification = c(0, 0.5),  # 图例左对齐，垂直居中
                    legend.text = element_text(size = 14),  # 设置图例文本大小
                    legend.title = element_text(size = 16, face = "bold"),
                    plot.margin = unit(c(0.75, 1.5, 0.75, 0.75), "lines"),
                    legend.key.height = unit(0.8, "cm")  # 增加这个值可以拉开行间距
                    # legend.key.width = unit(0.01, "cm")   # 如果需要也可以调整宽度
  )
  
}

# {
# 
#   library(usethis)
#   git_sitrep()
#   
#   usethis::create_github_token()
#   gitcreds::gitcreds_set()
#   
#   usethis::create_github_token()
#   gitcreds::gitcreds_set()
# }
### 2.数据整理 ####-------------------------------------------------------------
setwd(Output_derived_data)
df_lab_log <- readRDS("mimic4_Sepsis住院天数大于2小于40的患者脓毒症相关指标log.rds")
 
{
  library(dplyr)
  # ,"Creatine Kinase, MB Isoenzyme"
  # 1. 定义你需要的目标指标名称 (请务必确保这里的名称与 df_lab$label 中的完全一致)
  target_all <- c("Troponin T", "Creatine Kinase, MB Isoenzyme")
  # "CK-MB Index" 
  # "Creatine Kinase, MB Isoenzyme"
  # 2. 执行筛选逻辑
  # 我们先通过分组计算，找出符合条件的 stay_id（或 subject_id）
  eligible_ids <- df_lab_log %>%
    # filter(TIME <= 7 ) %>% 
    # 只保留这三项指标的数据，且数值不为空
    filter(label %in% target_all, !is.na(valuenum)) %>%
    # 按患者和指标分组，统计每人每项指标测了多少次
    group_by(hadm_id, label) %>%
    summarise(measure_count = n(), .groups = "drop") %>%
    # 筛选出测量次数 >= 3 的记录
    filter(measure_count >= 3) %>%
    # 再次按患者分组，看该患者现在还剩下几种指标
    group_by(hadm_id) %>%
    summarise(label_count = n_distinct(label)) %>%
    # 只保留那些依然拥有全部 3 种指标的患者
    filter(label_count == length(target_all)) %>%
    pull(hadm_id)
  
  # 3. 使用这些 ID 过滤原始数据集
  df_final <- df_lab_log %>%
    filter(hadm_id %in% eligible_ids)
  
  # 查看筛选结果
  print(paste("符合条件的患者数量:", length(eligible_ids)))
}
 
final_df4 <- df_final
 
setwd(Output_derived_data)
saveRDS(final_df4,"mimic4_sepsis_ckmb_ctnt_GE3_相关指标.rds")
df_target_id <- final_df4 %>% distinct(hadm_id) %>% select(hadm_id) 
write.csv(df_target_id,"df_target_id_heart.csv",row.names = FALSE,quote = FALSE)
df_cs_id <- final_df4 %>% filter(Cardiogenic_shock == 1) %>% distinct(hadm_id) %>% select(hadm_id)  
length(df_cs_id$hadm_id)
# [1] 120
# 检视用药信息 床边监护 诊断 #-------------------------------------------------------
# 用药信息
setwd(Input_derived_data)
df_input <- read.xlsx("mimic4_heart_sepsis_用药信息.xlsx")
length(unique(df_input$hadm_id))
# [1] 1541

length(unique(subset(df_input, hadm_id %in% df_cs_id$hadm_id)$hadm_id))
# [1] 119

colnames(df_input) <- make.unique(colnames(df_input))

# 床边监护
setwd(Input_derived_data)
df_chartevents <- read.xlsx("mimic4_heart_sepsis_床边监护.xlsx") 
df_item <- read.xlsx("mimic4_d_items.xlsx")
length(unique(df_chartevents$hadm_id))
# [1] 1544

length(unique(subset(df_chartevents, hadm_id %in% df_cs_id$hadm_id)$hadm_id))
# [1] 119 
colnames(df_chartevents)

df_chartevents1 <- df_chartevents %>% 
  mutate(across(c(charttime, storetime), ~ as.POSIXct(as.numeric(.)*86400, origin = "1899-12-30", tz = "UTC"))) %>% 
  left_join(df_item,by = "itemid")

unique(df_chartevents1$label)

# 1. 获取两个表的唯一 ID 集合
ids_meds <- unique(df_input$hadm_id) 
ids_charts <- unique(df_chartevents$hadm_id)

# 2. 取交集 (重合的 ID)
common_ids <- intersect(ids_meds, ids_charts)
cat("重合的患者 ID 数量为：", length(common_ids), "\n")
# 重合的患者 ID 数量为： 1541 

# 3. 查看这 119 个心源性休克 (CS) 患者是否都在重合名单里
cs_ids <- unique(df_cs_id$hadm_id)
common_cs_ids <- intersect(common_ids, cs_ids)
cat("重合名单中包含的心源性休克患者数量为：", length(common_cs_ids), "\n")
# 重合名单中包含的心源性休克患者数量为： 119 

# 输出同时有床边监护和用药情况id的事件表
colnames(df_input)
df_admit <- final_df4 %>% distinct(hadm_id,admittime,edregtime)
df_input1 <- df_input %>% 
  filter(hadm_id %in% common_ids) %>% 
  # 关键步骤：将左表的 hadm_id 转换为字符型
  mutate(hadm_id = as.character(hadm_id)) %>% 
  mutate(across(c(starttime, endtime, storetime), 
                ~ as.POSIXct(as.numeric(.) * 86400, origin = "1899-12-30", tz = "UTC"))) %>% 
  # 同时也确保右表关联列也是字符型（虽然报错提示它已经是了，但这样写更稳健）
  left_join(df_admit %>% mutate(hadm_id = as.character(hadm_id)), 
            by = "hadm_id")

setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered") 
saveRDS(df_input1,"mimic4_sepsis_heart_filtered_用药信息.rds")

df_chartevents2 <- df_chartevents1 %>% filter(hadm_id %in% common_ids) %>% 
  mutate(hadm_id = as.character(hadm_id)) %>%    
  left_join(df_admit %>% mutate(hadm_id = as.character(hadm_id)), by = "hadm_id")

setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered") 
saveRDS(df_chartevents2,"mimic4_sepsis_heart_filtered_床边监护.rds")
unique(df_chartevents3$label)

setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered") 
df_input1 <- readRDS("mimic4_sepsis_heart_filtered_用药信息.rds") 
df_chartevents2 <- readRDS("mimic4_sepsis_heart_filtered_床边监护.rds")
df_input2 <- df_input1 %>% mutate(hosp_start_time = pmin(edregtime, admittime, na.rm = TRUE)) %>%
  mutate(STARTTIME = as.numeric(difftime(starttime, hosp_start_time, units = "days")),
         ENDTIME = as.numeric(difftime(endtime, hosp_start_time, units = "days")),
         STORETIME = as.numeric(difftime(endtime, hosp_start_time, units = "days"))
  )
df_chartevents3 <- df_chartevents2 %>% mutate(hosp_start_time = pmin(edregtime, admittime, na.rm = TRUE)) %>%
  mutate(TIME = as.numeric(difftime(charttime , admittime, units = "days")),STORETTIME = as.numeric(difftime(storetime , admittime, units = "days")))

patient_count_per_label <- df_chartevents3 %>%
  group_by(label) %>%
  summarise(unique_patients = n_distinct(hadm_id)) %>%
  arrange(desc(unique_patients))

patient_count_per_label1 <- df_input2 %>%
  group_by(label) %>%
  summarise(unique_patients = n_distinct(hadm_id)) %>%
  arrange(desc(unique_patients))


patient_count_per_label2 <- final_df4 %>%
  group_by(label) %>%
  summarise(unique_patients = n_distinct(hadm_id)) %>%
  arrange(desc(unique_patients))

{
  abp_labels <- c(
    "Arterial Blood Pressure mean",
    "Arterial Blood Pressure diastolic", 
    "Arterial Blood Pressure systolic"
  )
  
  # 直接计算有ABP的患者数
  abp_patients <- df_chartevents3 %>%
    filter(label %in% abp_labels) %>%
    pull(hadm_id) %>%
    unique()
  
  total_patients <- df_chartevents3 %>%
    pull(hadm_id) %>%
    unique()
  
  cat("总患者数:", length(total_patients))
  cat("至少有1种有创动脉血压的患者数:", length(abp_patients))
  cat("覆盖率:", round(length(abp_patients)/length(total_patients)*100, 2), "%")
  cat("无任何有创动脉血压的患者数:", length(total_patients) - length(abp_patients))
  }
# 总患者数: 1541至少有1种有创动脉血压的患者数: 911覆盖率: 59.12 %无任何有创动脉血压的患者数: 630

# 诊断
setwd(Input_derived_data)
df_diagnose <- read.xlsx("mimic4_spesis_诊断信息.xlsx") 
# 自动修复重复的列名
names(df_diagnose) <- make.names(names(df_diagnose), unique = TRUE)

# 查看修复后的列名
colnames(df_diagnose)
df_diagnose <- df_diagnose %>% select(,c("subject_id","hadm_id","seq_num","icd_code","icd_version","long_title") )
# 现在可以正常筛选了
df_diagnose1 <- df_diagnose %>% filter(hadm_id %in% common_ids) 

length(unique(df_diagnose1$hadm_id))
# [1] 1541

unique(df_diagnose1$long_title)

setwd(Output_derived_data)
write.xlsx(df_diagnose1,"mimic4_sepsis_heart_诊断.xlsx")
write.xlsx(df_diagnose7,"mimic4_sepsis_cs_nosurgery_诊断.xlsx")

setwd(Output_derived_data)
df_diagnose1 <- read.xlsx("mimic4_sepsis_heart_诊断.xlsx") 
df_diagnose6 <- df_diagnose1 %>% filter(hadm_id %in% df_cs_id$hadm_id)

# seq_num=1 是sepsis|septic患者
df_diagnose2 <- df_diagnose1 %>% filter(seq_num == 1) %>% filter(grepl("sepsis|septic", long_title, ignore.case = TRUE))
df_diagnose3 <- df_diagnose1 %>% filter(hadm_id %in% df_diagnose2$hadm_id)
length(unique(df_diagnose3$hadm_id))
# 869
# 
# # 基于你提供的 patients_with_surgery$long_title 原文
# # 筛选出所有与心脏手术相关的诊断（直接使用原文）
# {
#   cardiac_surgery_diagnoses <- c(
#     # 明确的心脏手术状态
#     "Aortocoronary bypass status",
#     "Presence of aortocoronary bypass graft",
#     "Presence of coronary angioplasty implant and graft",
#     "Heart valve replaced by transplant",
#     "Presence of other heart-valve replacement",
#     "Presence of cardiac pacemaker",
#     "Cardiac pacemaker in situ",
#     "Presence of automatic (implantable) cardiac defibrillator",
#     "Automatic implantable cardiac defibrillator in situ",
#     
#     # 心脏手术并发症
#     "Postoperative shock, cardiogenic",
#     "Postprocedural cardiogenic shock, initial encounter",
#     "Postprocedural heart failure following cardiac surgery",
#     "Postprocedural cardiac arrest following cardiac surgery",
#     "Functional disturbances following cardiac surgery",
#     "Atherosclerosis of autologous vein bypass graft",
#     "Atherosclerosis of artery bypass graft",
#     "Coronary atherosclerosis of autologous vein bypass graft",
#     "Coronary atherosclerosis of artery bypass graft",
#     "Mechanical complication of automatic implantable cardiac defibrillator",
#     "Mechanical complication of cardiac pacemaker",
#     "Infection and inflammatory reaction due to cardiac device, implant, and graft",
#     "Thrombosis due to cardiac prosthetic devices, implants and grafts, initial encounter",
#     "Other complications due to other cardiac device, implant, and graft",
#     "Other specified complication of cardiac prosthetic devices, implants and grafts, initial encounter",
#     "Unspecified complication of cardiac and vascular prosthetic device, implant and graft, initial encounter",
#     "Hemorrhage due to cardiac prosthetic devices, implants and grafts, initial encounter",
#     "Stenosis of other cardiac prosthetic devices, implants and grafts, initial encounter",
#     "Atherosclerosis of autologous vein coronary artery bypass graft(s) with unspecified angina pectoris",
#     "Atherosclerosis of autologous artery coronary artery bypass graft(s) with unspecified angina pectoris",
#     "Atherosclerosis of autologous vein coronary artery bypass graft(s) with unstable angina pectoris",
#     "Atherosclerosis of coronary artery bypass graft(s) without angina pectoris",
#     "Other postprocedural cardiac functional disturbances following other surgery",
#     "Postoperative shock, other",  # 需结合上下文判断是否为心脏术后
#     "Postoperative shock, unspecified",
#     "Postprocedural hypotension",
#     "Postprocedural cardiac arrest following other surgery",
#     "Other intraoperative cardiac functional disturbances during other surgery",
#     "Intraoperative cardiac arrest during other surgery"
#   )
#   
#   # 从诊断表中找出有心脏手术相关诊断的患者
#   patients_with_cardiac_surgery <- df_diagnose3 %>%
#     filter(long_title %in% cardiac_surgery_diagnoses) %>%
#     distinct(hadm_id)
# 
#   cat("找到有心脏手术史的患者数:", nrow(patients_with_cardiac_surgery), "\n")
#   
#   # 去除这些患者
#   df_diagnose2 <- df_diagnose1 %>%
#     filter(!hadm_id %in% patients_with_cardiac_surgery$hadm_id)
# 
#   # 查看结果
#   cat("原始休克患者数:", length(unique(df_diagnose1$hadm_id)), "\n")
#   cat("去除心脏手术患者后:", length(unique(df_diagnose2$hadm_id)), "\n")
#   cat("去除的患者数:", length(unique(df_diagnose1$hadm_id)) - length(unique(df_diagnose2$hadm_id)), "\n")
# } 
# # 找到有心脏手术史的患者数: 439 
# # 原始休克患者数: 1541 
# # 去除心脏手术患者后: 1102   
# # 去除的患者数: 439 

df_shock <- df_diagnose3 %>%
  filter(grepl("shock", long_title, ignore.case = TRUE)) %>%
  filter(long_title != "Severe sepsis without septic shock")

length(unique(df_shock$hadm_id)) 
# [1] 648
df_cs_only <- df_diagnose3 %>% filter(hadm_id %in% df_cs_id$hadm_id)
length(unique(df_cs_only$hadm_id)) 
# [1] 49
df_shock_only <- df_shock %>% filter(!hadm_id %in% df_cs_only$hadm_id)
length(unique(df_shock_only$hadm_id)) 
# [1] 599

unique(df_shock$long_title)
# [1] "Severe sepsis with septic shock"                    
# [2] "Cardiogenic shock"                                  
# [3] "Septic shock"                                       
# [4] "Other shock without mention of trauma"              
# [5] "Hypovolemic shock"                                  
# [6] "Other shock"                                        
# [7] "Postoperative shock, cardiogenic"                   
# [8] "Shock, unspecified"                                 
# [9] "Other postprocedural shock, initial encounter"      
# [10] "Postprocedural cardiogenic shock, initial encounter"

# 确定休克开始时间#-------------------------------------------------------------
{
  setwd(Input_derived_data)
  sofa_cohort_master <- read.xlsx("sofa_cohort_master.xlsx")
  df_filter_input <- df_input2 %>% filter(hadm_id %in% df_diagnose3$hadm_id)
  df_filter_chartevent <- df_chartevents3 %>% filter(hadm_id %in% df_diagnose3$hadm_id)
  df_filter_lab <- final_df4 %>% filter(hadm_id  %in% df_diagnose3$hadm_id)
  gc()
  df_s_input <- df_filter_input %>% filter(hadm_id %in% df_shock_only$hadm_id)
  df_s_chartevent <- df_filter_chartevent %>% filter(hadm_id %in% df_shock_only$hadm_id)
  length(unique(df_s_chartevent$hadm_id))
  df_s_lab <- df_filter_lab %>% filter(hadm_id %in% df_shock_only$hadm_id)
  
  # 定义液体 ItemID (晶体液/胶体液)
  fluid_itemids <- c(225158, 225828, 225943, 225159, 220864, 220862, 225161, 225823, 225825)
  
  df_fluid_600_time <- df_s_input %>%
    filter(itemid %in% fluid_itemids) %>%
    arrange(stay_id, starttime) %>%
    group_by(stay_id) %>%
    # 计算每个 stay_id 的累积入量
    mutate(cum_amount = cumsum(amount)) %>%
    # 找出第一次超过 600 mL 的那行记录
    filter(cum_amount >= 600) %>%
    summarise(fluid_600_time = min(endtime, na.rm = TRUE)) %>%
    ungroup()
  
  vasso_itemids <- c(221906, 222315, 221662, 221289, 221749)
  # 1. 221906 —— 去甲肾上腺素 (Norepinephrine)
  # 2. 222315 —— 血管加压素 (Vasopressin)
  # 3. 221662 —— 多巴胺 (Dopamine)
  # 4. 221289 —— 肾上腺素 (Epinephrine)
  # 5. 221749 —— 去氧肾上腺素 / 新福林 (Phenylephrine)
  df_vaso_start_time <- df_s_input %>%
    filter(itemid %in% vasso_itemids) %>%
    group_by(stay_id) %>%
    summarise(vaso_start_time = min(starttime, na.rm = TRUE)) %>%
    ungroup()
  
  df_lactate_2_time <- df_s_lab_with_stay %>%
    filter(grepl("Lactate", label, ignore.case = TRUE) & valuenum >= 2) %>%
    group_by(stay_id) %>%
    summarise(lactate_2_time = min(charttime, na.rm = TRUE)) %>%
    ungroup()
  
  # 1. 确保 df_shock 拥有 stay_id
  # 我们从之前的 sofa 队列中提取 ID 对应关系
  df_shock_with_id <- df_shock %>%
    # 必须关联，因为后续的 join 都依赖 stay_id
    inner_join(
      sofa_cohort_master %>% select(subject_id, hadm_id, stay_id) %>% distinct(),
      by = c("subject_id", "hadm_id")
    )
  
  # 2. 现在重新执行合并逻辑
  df_final_onset <- df_shock_with_id %>%
    # 此时可以 select stay_id 了
    select(subject_id, hadm_id, stay_id) %>%
    
    # 关联三个关键时间点
    left_join(df_fluid_600_time, by = "stay_id") %>%
    left_join(df_vaso_start_time, by = "stay_id") %>%
    left_join(df_lactate_2_time, by = "stay_id") %>%
    
    # 3. 过滤出三个条件都满足的人（严格遵循 Sepsis-3 诊断标准）
    filter(!is.na(fluid_600_time) & !is.na(vaso_start_time) & !is.na(lactate_2_time)) %>%
    
    # 4. 计算发生时刻
    mutate(
      # shock_onset_time 定义为：补液够了、药上了、乳酸高了，这三者最后达成的那一刻
      shock_onset_time = pmax(fluid_600_time, vaso_start_time, lactate_2_time, na.rm = TRUE)
    )
  
  # 打印最终确定的休克人数
  cat("--- 最终分析报告 ---\n")
  cat("原始休克诊断人数:", nrow(df_shock), "\n")
  # 原始休克诊断人数: 764 
  cat("符合[补液+升压药+高乳酸]严格判定的人数:", nrow(df_final_onset), "\n")
  # 符合[补液+升压药+高乳酸]严格判定的人数: 453 
  
  length(unique(df_final_onset$hadm_id))
  # [1] 371
}


# 确定cs开始时间#-------------------------------------------------------------
df_cs_input <- df_filter_input %>% filter(hadm_id %in% df_cs_only$hadm_id)
df_cs_chartevent <- df_filter_chartevent %>% filter(hadm_id %in% df_cs_only$hadm_id)
length(unique(df_cs_chartevent$hadm_id))
df_cs_lab <- df_filter_lab %>% filter(hadm_id %in% df_cs_only$hadm_id)
setwd(Input_derived_data)
sofa2_support_with_labels <- read.xlsx("sofa2_support_with_labels.xlsx")
head(sofa2_support_with_labels)
sofa2_support_with_labels <- sofa2_support_with_labels %>%    
  inner_join(
    sofa_cohort_master %>% select(hadm_id, stay_id) %>% distinct(),
    by = c("stay_id")
  )
df_cs_support <- sofa2_support_with_labels %>% filter(hadm_id %in% df_cs_only$hadm_id)

library(dplyr)
library(lubridate)
library(openxlsx)

# 提取这些精细指标的示例代码
df_cs_hemo <- df_cs_chartevent %>%
  filter(itemid %in% c(
    228368, 228177, # CI
    220060,         # PCWP (PAWP)
    220048, 228174  # CO (用于计算 CPO)
  )) 

# --- 0. 基础 ID 清理函数 ---
clean_ids <- function(df) {
  df %>% mutate(across(any_of(c("subject_id", "hadm_id", "stay_id")), as.numeric))
}

# --- 1. 数据准备与 ID 补全 ---
# 确保主队列干净，移除可能存在的旧 stay_id
df_cs_only <- clean_ids(df_cs_only) %>% select(-any_of("stay_id")) 
sofa_cohort_master_clean <- sofa_cohort_master %>% select(subject_id, hadm_id, stay_id, intime, outtime) %>% distinct() %>% clean_ids()

# 补全实验室检查表的 stay_id (因为 labevents 默认不带 stay_id)
df_cs_lab_with_id <- df_cs_lab %>%
  clean_ids() %>%
  inner_join(sofa_cohort_master_clean, by = c("subject_id", "hadm_id")) %>%
  filter(charttime >= (intime - hours(24)) & charttime <= outtime)

# 补全机械辅助表的 ID
sofa2_support_with_labels <- clean_ids(sofa2_support_with_labels) %>%
  inner_join(sofa_cohort_master_clean %>% select(hadm_id, stay_id), by = "stay_id")

# --- 2. 提取四个维度的核心时间点 ---

# A. 强心干预时间 (Dobutamine/Milrinone)
df_cs_intervention <- df_cs_input %>%
  clean_ids() %>%
  filter(itemid %in% c(221653, 221986)) %>% 
  group_by(stay_id) %>%
  summarise(time_intervention = min(starttime, na.rm = TRUE)) %>%
  ungroup()

# B. 循环受损时间 (低血压: SBP<90 或 MAP<65)
df_cs_low_bp <- df_cs_chartevent %>%
  clean_ids() %>%
  filter(
    (grepl("systolic", label, ignore.case = TRUE) & valuenum > 0 & valuenum < 90) |
      (grepl("mean", label, ignore.case = TRUE) & valuenum > 0 & valuenum < 65)
  ) %>%
  group_by(stay_id) %>%
  summarise(time_bp_low = min(charttime, na.rm = TRUE)) %>%
  ungroup()

# C. 组织灌注不足时间 (乳酸 >= 2)
df_cs_lactate <- df_cs_lab_with_id %>%
  filter(grepl("Lactate", label, ignore.case = TRUE) & valuenum >= 2) %>%
  group_by(stay_id) %>%
  summarise(time_lactate_high = min(charttime, na.rm = TRUE)) %>%
  ungroup()

# D. 机械辅助（IABP, Impella, VAD, ECMO）启动时间
df_mcs_onset <- sofa2_support_with_labels %>%
  filter(is_iabp == 1 | is_impella == 1 | is_vad == 1 | is_ecmo == 1) %>%
  group_by(stay_id) %>%
  mutate(across(c(charttime), ~ as.POSIXct(as.numeric(.)*86400, origin = "1899-12-30", tz = "UTC"))) %>% 
  summarise(time_mcs_start = min(charttime, na.rm = TRUE)) %>%
  ungroup()
length(unique(df_mcs_onset$hadm_id))
# --- 3. 汇总判定心源性休克开始时间 (Onset) ---

df_cs_final_onset <- df_cs_only %>%
  # 建立 stay_id 桥梁
  inner_join(
    sofa_cohort_master_clean %>% select(subject_id, hadm_id, stay_id), 
    by = c("subject_id", "hadm_id"),
    relationship = "many-to-many"
  ) %>%
  
  # 关联四个维度的时间轴
  left_join(df_cs_intervention, by = "stay_id") %>% 
  left_join(df_mcs_onset, by = "stay_id") %>%       
  left_join(df_cs_low_bp, by = "stay_id") %>%       
  left_join(df_cs_lactate, by = "stay_id") %>%      
  
  mutate(
    # 步骤 A: 确定“心脏泵衰竭”临床支持的最早时间 (强心药 或 机械辅助)
    time_cardiac_support = pmin(time_intervention, time_mcs_start, na.rm = TRUE),
    
    # 步骤 B: 确定“循环受损”的最早临床表现 (低血压 或 高乳酸)
    time_shock_indices = pmin(time_bp_low, time_lactate_high, na.rm = TRUE),
    
    # 步骤 C: 最终心源性休克开始时间 (当心脏支持已开始 且 循环指标已变差，取两者较晚者)
    cs_onset_time = pmax(time_cardiac_support, time_shock_indices, na.rm = TRUE)
  ) %>%
  
  # 过滤掉没有任何临床支持或休克证据的人
  filter(!is.na(cs_onset_time)) %>%
  
  # 针对同一个 hadm_id 进出多次 ICU 的情况，取第一次发生休克的时间
  group_by(hadm_id) %>%
  slice_min(cs_onset_time, n = 1, with_ties = FALSE) %>%
  ungroup()

# --- 4. 结果核查 ---
cat("心源性休克队列原始人数:", length(unique(df_cs_only$hadm_id)), "\n")
cat("成功判定 Onset 时间的人数:", length(unique(df_cs_final_onset$hadm_id)), "\n")

# 查看前几行确认时间格式
head(df_cs_final_onset %>% select(subject_id, hadm_id, stay_id, cs_onset_time))

df_diagnose3_cs <- df_diagnose3 %>% filter(hadm_id %in% df_cs_id$hadm_id)

setwd(Output_derived_data)
write.xlsx(df_diagnose3_cs,"df_diagnose3_cs.xlsx")


# 确定所有休克群体休克开始时间#-------------------------------------------------------------
setwd(Input_derived_data)
sofa_cohort_master <- read.xlsx("sofa_cohort_master.xlsx")
df_filter_input <- df_input2 %>% filter(hadm_id %in% df_diagnose3$hadm_id)
df_filter_chartevent <- df_chartevents3 %>% filter(hadm_id %in% df_diagnose3$hadm_id)
df_filter_lab <- final_df4 %>% filter(hadm_id %in% df_diagnose3$hadm_id)
gc()
df_s_all_input <- df_filter_input %>% filter(hadm_id %in% df_shock$hadm_id)
df_s_all_chartevent <- df_filter_chartevent %>% filter(hadm_id %in% df_shock$hadm_id)
length(unique(df_s_all_chartevent$hadm_id))
df_s_all_lab <- df_filter_lab %>% filter(hadm_id %in% df_shock$hadm_id)

# 补液+乳酸+血管升压素
{
  library(dplyr)
  library(lubridate)
  
  
  # 0. 环境准备与 ID 统一
  
  clean_ids <- function(df) {
    df %>% mutate(across(any_of(c("subject_id", "hadm_id", "stay_id")), as.numeric))
  }
  
  # 确保核心表 ID 干净
  df_shock_base <- clean_ids(df_shock) %>% select(-any_of("stay_id"))
  sofa_bridge <- sofa_cohort_master %>% 
    select(subject_id, hadm_id, stay_id, intime, outtime) %>% 
    distinct() %>% 
    clean_ids()
  
  
  # 1. 确定锚点：血管活性药物首次启动时间
  
  # 定义升压药 ItemID (去甲、多巴胺、血管加压素、肾上腺素、去氧肾)
  vasso_itemids <- c(221906, 221662, 222315, 221289, 221749)
  
  df_vaso_anchor <- df_s_all_input  %>%
    filter(itemid %in% vasso_itemids) %>%
    group_by(stay_id) %>%
    summarise(vaso_start_time = min(starttime, na.rm = TRUE)) %>%
    ungroup() %>%
    clean_ids()
  
  
  # 2. 验证条件1：药前 6 小时补液量 (基于 rate 精确计算)
  
  # 定义复苏液体 ItemID (晶体液/胶体液)
  fluid_itemids <- c(
    # 晶体液
    225158, 225828, 227073, 220861, 225823, 225825, 225159, 227533,
    # 胶体液
    220864, 220862, 220970,
    # 血制品
    225168, 227070, 220996, 225171, 225170
  )
  
  df_fluid_verified <- df_s_all_input  %>%
    filter(itemid %in% fluid_itemids) %>%
    inner_join(df_vaso_anchor, by = "stay_id") %>%
    # 定义 6 小时分析窗口：[药前6h, 药开始]
    mutate(
      win_start = vaso_start_time - hours(12),
      win_end   = vaso_start_time
    ) %>%
    # 筛选与窗口有时间交集的行
    filter(starttime <= win_end & endtime >= win_start) %>%
    # 计算相交部分的入量
    mutate(
      overlap_start = pmax(starttime, win_start),
      overlap_end   = pmin(endtime, win_end),
      overlap_hours = as.numeric(difftime(overlap_end, overlap_start, units = "hours")),
      # 计算总时长用于无流速时的分摊
      total_row_duration = as.numeric(difftime(endtime, starttime, units = "hours")),
      
      # 精确计算窗口内 mL 数
      amount_in_6h = case_when(
        !is.na(rate) & rate > 0 ~ rate * overlap_hours, # 有流速按流速算
        # total_row_duration > 0 ~ amount * (overlap_hours / total_row_duration), # 无流速按时间比例分摊(Bolus)
        total_row_duration > 0 ~ amount , # 无流速按时间比例分摊(Bolus)
        TRUE ~ amount # 瞬时记录直接取值
      )
    ) %>%
    group_by(stay_id, vaso_start_time) %>%
    summarise(total_fluid_6h = sum(amount_in_6h, na.rm = TRUE), .groups = "drop") %>%
    # 验证门槛：6 小时内补液 >= 600 mL
    filter(total_fluid_6h >= 600)
  length(unique(df_fluid_verified$stay_id))
  #184
  
  # 3. 验证条件2：药前后 6 小时高乳酸证据
  
  # 补全 Lab 表的 stay_id
  df_lab_with_stay <- df_s_all_lab %>%
    clean_ids() %>%
    inner_join(sofa_bridge, by = c("subject_id", "hadm_id"))
  # %>%
  #   filter(charttime >= (intime - hours(24)) & charttime <= outtime)
  
  df_lactate_verified <- df_lab_with_stay %>%
    filter(grepl("Lactate", label, ignore.case = TRUE) & valuenum >= 2) %>%
    inner_join(df_vaso_anchor, by = "stay_id") %>%
    # 验证窗口：药启动前后 6 小时 (共 12 小时)
    filter(charttime >= (vaso_start_time - hours(12)) & charttime <= (vaso_start_time + hours(12))) %>%
    group_by(stay_id, vaso_start_time) %>%
    # 记录该窗口内第一次乳酸达标的时间
    summarise(first_lac_time_in_6h = min(charttime, na.rm = TRUE), .groups = "drop")
  length(unique(df_lactate_verified$stay_id))
  # [1] 500
  
  # 4. 汇总判定休克发生时刻 (Shock Onset)
  
  df_final_shock_onset <- df_shock_base %>%
    # 关联 ID 桥梁
    inner_join(sofa_bridge %>% select(subject_id, hadm_id, stay_id), 
               by = c("subject_id", "hadm_id"), relationship = "many-to-many") %>%
    
    # 必须满足：有升压药启动
    inner_join(df_vaso_anchor, by = "stay_id") %>%
    
    # 必须满足：通过 6h 补液验证
    inner_join(df_fluid_verified, by = c("stay_id", "vaso_start_time")) %>%
    
    # 必须满足：通过 6h 乳酸验证
    inner_join(df_lactate_verified, by = c("stay_id", "vaso_start_time")) %>%
    
    mutate(
      # Time Zero 定义：在补液达标的前提下，药已上泵 且 乳酸已高 的那一刻
      # 取两者中最晚的一个 (pmax)
      shock_onset_time = pmax(vaso_start_time, first_lac_time_in_6h, na.rm = TRUE)
    ) %>% 
    
    # 针对同一住院多次进出 ICU，取最早的一次发病
    group_by(hadm_id) %>%
    slice_min(shock_onset_time, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  
  
  # 5. 输出最终报告
  
  cat("==========================================\n")
  cat("      脓毒性休克队列清洗报告 (Strict 12h)\n")
  cat("==========================================\n")
  cat("1. 原始诊断人数 (df_shock):      ", length(unique(df_shock_base$hadm_id)), "\n")
  cat("2. 启动升压药治疗人数:           ", nrow(df_vaso_anchor), "\n")
  cat("3. 通过补液+乳酸双重校验人数:    ", nrow(df_final_shock_onset), "\n")
  cat("------------------------------------------\n")
  cat("最终确定的 Time Zero 中位时间点（前5行）:\n")
  print(head(df_final_shock_onset %>% select(subject_id, stay_id, shock_onset_time)))
}
# 1. 原始诊断人数 (df_shock):       648 
# 2. 启动升压药治疗人数:            647 
# 3. 通过补液+乳酸双重校验人数:     113 

# 低血压MAP <= 65 + 前后12h乳酸Lactate >= 2
{
  
  # 1. 提取所有低血压候选时刻 (MAP <= 65)
  
  # 220052: Arterial Blood Pressure mean
  # 220181: Non Invasive Blood Pressure Mean
  df_hypotension_candidates <- df_s_all_chartevent %>%
    filter(itemid %in% c(220052, 220181) & valuenum > 0 & valuenum <= 65) %>%
    select(subject_id, hadm_id, stay_id, map_time = charttime, map_value = valuenum) %>%
    mutate(across(c(subject_id, hadm_id, stay_id), as.numeric))
  df_hypotension_candidates1 <- df_s_all_chartevent %>%
    filter(itemid %in% c(220052, 220181)) %>%
    select(subject_id, hadm_id, stay_id, map_time = charttime, map_value = valuenum) %>%
    mutate(across(c(subject_id, hadm_id, stay_id), as.numeric)) %>% 
    filter(hadm_id %in% df_cs_id$hadm_id)
  
  length(unique(df_hypotension_candidates$hadm_id))
  # [1] 368
  
  # 2. 提取所有高乳酸记录 (Lactate >= 2)
  
  # 这里的 df_lab_with_stay 是补全了 stay_id 的实验室表
  df_high_lactate <- df_lab_with_stay %>%
    filter(grepl("Lactate", label, ignore.case = TRUE) & valuenum >= 2) %>%
    select(stay_id, lactate_time = charttime, lactate_value = valuenum) %>%
    mutate(stay_id = as.numeric(stay_id))
  length(unique(df_high_lactate$stay_id))
  # 775
  
  # 3. 按照文章逻辑进行时空匹配
  
  df_shock_onset_tsai <- df_hypotension_candidates %>%
    # 将低血压时刻与乳酸记录进行关联 (根据 stay_id)
    inner_join(df_high_lactate, by = "stay_id", relationship = "many-to-many") %>%
    
    # 计算低血压时刻与乳酸检测时刻的时间差
    mutate(time_diff_hrs = as.numeric(difftime(map_time, lactate_time, units = "hours"))) %>%
    
    # 关键验证条件：乳酸必须在低血压时刻的前后 12 小时内
    filter(abs(time_diff_hrs) <= 12) %>%
    
    # 确定每个患者符合条件的第一个低血压时刻
    group_by(stay_id) %>%
    summarise(
      shock_onset_time = min(map_time, na.rm = TRUE),
      verified_lactate_val = first(lactate_value), # 仅作记录
      .groups = "drop"
    )
  
  
  # 4. 汇总最终队列
  
  df_final_shock_cohort <- df_shock_base %>%
    inner_join(sofa_bridge %>% select(subject_id, hadm_id, stay_id), 
               by = c("subject_id", "hadm_id")) %>%
    # 只保留符合上述 Tsai 逻辑（MAP低且有乳酸验证）的人
    inner_join(df_shock_onset_tsai, by = "stay_id") %>%
    
    # 针对同一个 hadm_id，保留最早的休克开始时间
    group_by(hadm_id) %>%
    slice_min(shock_onset_time, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # 输出结果报告
  cat("==========================================\n")
  cat("      Tsai (2025) 逻辑休克判定报告\n")
  cat("==========================================\n")
  cat("1. 原始诊断人数:            ", length(unique(df_shock_base$hadm_id)), "\n")
  cat("2. 符合 MAP<=65 且 12h内乳酸>=2 的人数: ", nrow(df_final_shock_cohort), "\n")
  cat("------------------------------------------\n")
  print(head(df_final_shock_cohort %>% select(subject_id, stay_id, shock_onset_time)))
}
# 1. 原始诊断人数:             648 
# 2. 符合 MAP<=65 且 12h内乳酸>=2 的人数:  328

# 血压+乳酸+血管升压素
{
  {
    library(dplyr)
    library(lubridate)
    
    
    # 0. 基础 ID 统一与清理
    
    clean_ids <- function(df) {
      df %>% mutate(across(any_of(c("subject_id", "hadm_id", "stay_id")), as.numeric))
    }
    
    df_shock_base <- clean_ids(df_shock_base) %>% select(-any_of("stay_id"))
    sofa_bridge <- sofa_cohort_master %>% 
      select(subject_id, hadm_id, stay_id, intime, outtime) %>% 
      distinct() %>% 
      clean_ids()
    
    
    # 1. 提取锚点：首次启动血管活性药物
    
    vasso_itemids <- c(221906, 221662, 222315, 221289, 221749)
    
    df_vaso_candidates <- df_s_all_input %>%
      filter(itemid %in% vasso_itemids) %>%
      group_by(stay_id) %>%
      summarise(vaso_start_time = min(starttime, na.rm = TRUE)) %>%
      ungroup() %>%
      clean_ids()
    
    # ------------------------------------------
    # 验证条件 3：排除入 ICU 后 2 小时内即启动药物的人
    # (为了保证有足够的休克前基线数据用于建模)
    # ------------------------------------------
    df_vaso_anchor <- df_vaso_candidates %>%
      inner_join(sofa_bridge, by = "stay_id") %>%
      mutate(hours_to_vaso = as.numeric(difftime(vaso_start_time, intime, units = "hours"))) %>%
      # filter(hours_to_vaso >= 2) %>%
      select(subject_id, hadm_id, stay_id, vaso_start_time, intime)
    
    
    # 2. 验证条件 1：药前后 12h 内有 MAP <= 65
    
    df_map_verified <- df_s_all_chartevent %>%
      filter(itemid %in% c(220052, 220181) & valuenum > 0 & valuenum <= 65) %>%
      clean_ids() %>%
      inner_join(df_vaso_anchor %>% select(stay_id, vaso_start_time), by = "stay_id") %>%
      filter(abs(difftime(charttime, vaso_start_time, units = "hours")) <= 12) %>%
      distinct(stay_id) # 只要存在即可
    
    
    # 3. 验证条件 2：药前后 12h 内有 Lactate >= 2
    
    df_lac_verified <- df_lab_with_stay %>%
      filter(grepl("Lactate", label, ignore.case = TRUE) & valuenum >= 2) %>%
      clean_ids() %>%
      inner_join(df_vaso_anchor %>% select(stay_id, vaso_start_time), by = "stay_id") %>%
      filter(abs(difftime(charttime, vaso_start_time, units = "hours")) <= 12) %>%
      distinct(stay_id) # 只要存在即可
    
    
    # 4. 汇总最终队列与休克开始时间
    
    df_final_shock_cohort <- df_shock_base %>%
      # 必须在血管活性药物的候选池中
      inner_join(df_vaso_anchor, by = c("subject_id", "hadm_id")) %>%
      # 必须通过 MAP 验证
      inner_join(df_map_verified, by = "stay_id") %>%
      # 必须通过乳酸验证
      inner_join(df_lac_verified, by = "stay_id") %>%
      
      mutate(
        # 最终定义：药物启动时刻即为休克起始时刻 (T0)
        shock_onset_time = vaso_start_time
      ) %>%
      
      # 针对同一个 hadm_id，保留最早的休克开始记录
      group_by(hadm_id) %>%
      slice_min(shock_onset_time, n = 1, with_ties = FALSE) %>%
      ungroup()
    
    # 输出结果报告
    cat("==========================================\n")
    cat("      新逻辑：药物锚点+双重验证报告\n")
    cat("==========================================\n")
    cat("1. 初始诊断人数:              ", length(unique(df_shock_base$hadm_id)), "\n")
    # cat("2. 排除早期用药后剩余人数:    ", nrow(df_vaso_anchor), "\n")
    cat("3. 最终符合药+低血压+高乳酸人数: ", nrow(df_final_shock_cohort), "\n")
    cat("------------------------------------------\n")
    # 预览对齐后的时间
    print(head(df_final_shock_cohort %>% select(hadm_id, stay_id, shock_onset_time)))
  }
}
# 1. 初始诊断人数:               648 
# 3. 最终符合药+低血压+高乳酸人数:  192 

# 合并可以确定休克开始时间的数据 #--------------------------------------------------------------------
shock_data <- final_df4 %>% filter(hadm_id %in% df_final_shock_cohort$hadm_id) %>% mutate(shock = 1)
length(unique(shock_data$hadm_id))
# 328
length(unique(subset(shock_data, Cardiogenic_shock == 1)$hadm_id))
noshock_data <- final_df4 %>% filter(hadm_id %in% df_diagnose3$hadm_id) %>% filter(!hadm_id %in% df_shock$hadm_id) %>% mutate(shock = 0,shock_onset_time = NA)
length(unique(noshock_data$hadm_id))
# 221
df_final_shock_cohort1 <- df_final_shock_cohort %>% select("hadm_id","shock_onset_time")
df_final_shock_cohort1 <- df_final_shock_cohort1 %>% mutate(hadm_id = as.character(hadm_id))

shock_data <- shock_data %>% left_join(df_final_shock_cohort1, by = "hadm_id")
all_data <- rbind(shock_data,noshock_data)      

colnames(all_data)  
# 定义需要转换的所有时间列名
time_cols <- c(
  "icu_intime_first", "charttime", "storetime", 
  "dischtime", "deathtime", "edregtime", "edouttime",  
  "icu_outtime_last", "shock_onset_time"
)

# 执行转换
all_data_rel <- all_data %>%
  mutate(
    hosp_start_time = as.POSIXct(hosp_start_time),
    
    # 批量将所有时间列转换为相对于hosp_start_time 的小时数
    across(
      all_of(time_cols), 
      ~ as.numeric(difftime(as.POSIXct(.), hosp_start_time, units = "days")),
      .names = "{.col}_rel" # 建议加上 _rel 后缀，保留原始列以备查验
    )
  )

# 3. 检查转换结果
# 现在你有了 shock_onset_time_rel 和 charttime_rel
colnames(all_data_rel)
# [1] "labevent_id"                "subject_id"                 "hadm_id"                   
# [4] "hosp_start_time"            "TIME"                       "LOS_DAYS"                  
# [7] "DEATH_DAYS"                 "icu_intime_first"           "specimen_id"               
# [10] "itemid"                     "charttime"                  "storetime"                 
# [13] "value"                      "valuenum"                   "valueuom"                  
# [16] "ref_range_lower"            "ref_range_upper"            "flag"                      
# [19] "priority"                   "comments"                   "label"                     
# [22] "fluid"                      "category"                   "Cardiogenic_shock"         
# [25] "admittime"                  "dischtime"                  "deathtime"                 
# [28] "admission_type"             "admission_location"         "discharge_location"        
# [31] "insurance"                  "language"                   "marital_status"            
# [34] "race"                       "edregtime"                  "edouttime"                 
# [37] "hospital_expire_flag"       "gender"                     "anchor_age"                
# [40] "anchor_year"                "anchor_year_group"          "dod"                       
# [43] "icu_stay_counts"            "icu_outtime_last"           "total_icu_los"             
# [46] "first_icu_unit"             "last_icu_unit"              "age"                       
# [49] "In_hospital_Mortality_Time" "Length_of_Stay"             "Gender"                    
# [52] "Outcome_Mortality"          "log_valuenum"               "shock"                     
# [55] "shock_onset_time"           "hosp_start_time_rel"        "icu_intime_first_rel"      
# [58] "charttime_rel"              "storetime_rel"              "dischtime_rel"             
# [61] "deathtime_rel"              "edregtime_rel"              "edouttime_rel"             
# [64] "icu_outtime_last_rel"       "shock_onset_time_rel"   
length(unique(all_data_rel$hadm_id[all_data_rel$Cardiogenic_shock == 1]))
# [1] 29
unique(all_data_rel$label)

setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
write.xlsx(all_data_rel,"lab_sepsis_seqnum1_filtered.xlsx")
# 判断icu结局 #------------------------------------------------------------
# 判断icu结局
{
  # 建立一个干净的 ID 映射表
  id_bridge <- sofa_cohort_master %>% 
    filter(hadm_id %in% all_data_rel$hadm_id) %>% 
    select(subject_id, hadm_id, stay_id, intime, outtime) %>% 
    distinct() %>%
    mutate(across(c(subject_id, hadm_id, stay_id), as.numeric))
  # 补全 stay_id
  nrow(id_bridge)
  # [1] 664
  
  # 患者住院期间进出icu次数
  hadm_stay_counts <- id_bridge %>% 
    group_by(hadm_id) %>%
    summarise(num_unique_stays = n_distinct(stay_id, na.rm = TRUE)) %>%
    arrange(desc(num_unique_stays))
  
  # 查看前几行
  print(head(hadm_stay_counts))
  length(unique(hadm_stay_counts$hadm_id))
  # [1] 549
  all_data_outcome  <- all_data_rel %>% distinct(hadm_id,hosp_start_time,shock_onset_time, deathtime,hospital_expire_flag)
  all_data_outcome$hadm_id <- as.numeric(all_data_outcome$hadm_id)
  all_data_outcome1 <- all_data_outcome %>% left_join(id_bridge,by = "hadm_id")
  # 判断icu结局
  icu_outcome <- all_data_outcome1 %>%
    distinct(hadm_id,stay_id,hosp_start_time,shock_onset_time,intime, outtime, deathtime,hospital_expire_flag) %>% 
    # 确保时间格式正确 
    mutate(across(c(intime, outtime), ~ as.POSIXct(as.numeric(.)*86400, origin = "1899-12-30", tz = "UTC"))) %>% 
    # --- 1. 判断休克是否发生在本次 ICU 期间 ---
    mutate(is_shock_onset_in_this_icu = case_when(
      is.na(shock_onset_time) ~ FALSE,  # 没休克的直接为 FALSE
      shock_onset_time >= intime & shock_onset_time <= outtime ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    
    # --- 2. 判断这行 stay_id 是第几次 ICU 记录 ---
    group_by(hadm_id) %>%
    # 按进入 ICU 时间排序
    arrange(intime, .by_group = TRUE) %>%
    # 使用 dense_rank 对 stay_id 进行排序编号 (1, 2, 3...)
    mutate(icu_stay_order = dense_rank(intime)) %>%
    ungroup() %>%
    
    # --- 3. 判断本次 ICU 结局 (死亡 vs 转出) ---
    mutate(this_icu_outcome = case_when(
      # 如果标记为死亡，且死亡时间在出 ICU 时间之前（加 1 小时容错防止系统记录延迟）
      hospital_expire_flag == 1 & !is.na(deathtime) & deathtime <= (outtime + hours(1)) ~ "Death in ICU",
      # 否则，视为活着转出 ICU（可能死在医院其他地方，或者康复出院）
      TRUE ~ "Discharged Alive from ICU"
    )) %>% 
    arrange(hadm_id,intime)
  
  # 统计检查结果
  # 1. 看看有多少休克确实发生在 ICU 记录窗口内
  cat("休克发生在对应 ICU 期间的人次数:\n")
  print(table(icu_outcome$is_shock_onset_in_this_icu))
  # FALSE  TRUE 
  # 338   326 
  # 2. 看看 ICU 结局分布
  cat("\n本次 ICU 住院结局分布:\n")
  print(table(icu_outcome$this_icu_outcome))
  # Death in ICU Discharged    Alive from ICU 
  # 160                          504 
  colnames(icu_outcome)
  
  # 患者在icu中休克住院号统计
  length(unique(subset(icu_outcome, is_shock_onset_in_this_icu)$hadm_id))
  # [1] 326
  colnames(icu_outcome)
  icu_outcome_noshock <- icu_outcome %>% filter(is.na(shock_onset_time))
}
# # A tibble: 6 × 2
# hadm_id num_unique_stays
# <dbl>            <int>
#   1 25640212                4
# 2 27477323                4
# 3 27546909                4
# 4 21170807                3
# 5 21422707                3
# 6 22811618                3
# 休克发生在对应 ICU 期间的人次数:
#   
#   FALSE  TRUE 
# 338   326 
# 
# 本次 ICU 住院结局分布:
#   
#   Death in ICU Discharged Alive from ICU 
# 160                       504

# 计算出室时间和死亡时间之间的绝对差值
{
  df_precise_check <- icu_outcome %>%
    filter(this_icu_outcome == "Death in ICU") %>%
    mutate(
      # 计算出室时间和死亡时间之间的绝对差值（分钟）
      time_diff_mins = as.numeric(abs(difftime(outtime, deathtime, units = "hours"))),
      
      # 如果差值在 1 分钟以内，可以认为它们在记录上是“同时”的
      is_exactly_same = (time_diff_mins <= 12)
    )
  
  # 查看有多少比例是“出室即死亡”
  summary(df_precise_check$time_diff_mins)
}
# 确定所有患者第一次进入icu的终点与结局
{
  df_first_icu <- icu_outcome %>%
    group_by(hadm_id) %>%
    # 按进入 ICU 时间排序，取最早的一次
    slice_min(intime, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  df_final_paths <- df_first_icu %>%
    mutate(
      # --- A. 确定终点时刻 (Endpoint Time) ---
      # 如果死在本次 ICU，终点是 deathtime；否则是离开 ICU 的 outtime
      endpoint_time = if_else(this_icu_outcome == "Death in ICU", deathtime, outtime),
      
      # --- B. 确定最终结局状态代码 ---
      # 状态定义：3 = 好转 (Recovery/Alive Discharge), 4 = 死亡 (Death in ICU)
      endpoint_state = if_else(this_icu_outcome == "Death in ICU", 4, 3),
      
      # --- C. 计算相对时间 (相对入室 intime 的小时数) ---
      # 终点时间（必须 > 0）
      t_end = as.numeric(difftime(endpoint_time, intime, units = "hours")),
      
      # 休克发生时间（仅对在本次 ICU 内休克的患者）
      # 如果休克发生在入室之前，我们将其标记为极小值(0.01)，表示初始即休克
      t_shock = case_when(
        is_shock_onset_in_this_icu == TRUE ~ as.numeric(difftime(shock_onset_time, intime, units = "hours")),
        TRUE ~ NA_real_
      ),
      t_shock = ifelse(t_shock <= 0, 0.01, t_shock) # 时间递增修正
    ) %>%
    # 逻辑清洗：剔除那些终点时间小于等于0的异常数据
    filter(t_end > 0)
  
  
  cat("--- 最终马尔可夫模型人群（仅限首诊ICU） ---\n")
  cat("总入组人数: ", nrow(df_final_paths), "\n")
  cat("期间发生休克人数: ", sum(df_final_paths$is_shock_onset_in_this_icu), "\n")
  cat("最终好转(存活出室)人数: ", sum(df_final_paths$endpoint_state == 3), "\n")
  cat("最终死亡(ICU内死亡)人数: ", sum(df_final_paths$endpoint_state == 4), "\n")
  
}
# --- 最终马尔可夫模型人群（仅限首诊ICU） ---
#   总入组人数:  549 
# 期间发生休克人数:  299 
# 最终好转(存活出室)人数:  420 
# 最终死亡(ICU内死亡)人数:  129 

# 输出首次icu临床终点
setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
write.xlsx(df_final_paths,"lab_sepsis_seqnum1_finalpath.xlsx")

# 二、住院全程数据检视 #-------------------------------------------------------- 
# 1.人口统计学数据和协变量分析-----------------------------------------------
setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
all_data_rel <- read.xlsx("lab_sepsis_seqnum1_filtered.xlsx")

length(unique(all_data_rel$hadm_id))
# [1] 549
cov1 <- all_data_rel[!duplicated(paste(all_data_rel$hadm_id)),]  #提取每个ID的第一行协变量值
colnames(cov1)
covname <- colnames(cov1[,c("Length_of_Stay","age","In_hospital_Mortality_Time")])
catcovname <- colnames(cov1[,c("Cardiogenic_shock","shock","Outcome_Mortality","Gender","race")])   #定义需要后续分析的连续协变量和分类协变量
colnames(cov1)
### 连续协变量描述性统计分析
summary(cov1[,covname])      #连续协变量描述性统计分析
# Min.   : 2.053   Min.   :20.00   Min.   : 2.053            
# 1st Qu.: 7.292   1st Qu.:61.00   1st Qu.: 5.089            
# Median :12.311   Median :71.00   Median :10.531            
# Mean   :14.231   Mean   :70.28   Mean   :11.363            
# 3rd Qu.:19.278   3rd Qu.:82.00   3rd Qu.:15.410            
# Max.   :39.625   Max.   :96.00   Max.   :37.211            
#                                  NA's   :367   

quantile(cov1$LOS_DAYS, probs = 0.95, na.rm = TRUE)
# 95% 
# 32.01986 
quantile(cov1$DEATH_DAYS, probs = 0.05,  na.rm = TRUE)
# 5% 
# 2.795208 
table(cov1$hospital_expire_flag[which(cov1$Cardiogenic_shock==1)]) # 心源性休克患者中死亡人数统计
# 0  1 
# 13 16 

### 分类协变量描述性统计分析
sapply(cov1[catcovname], table)  
sapply(cov1[,"gender"], table)  
# gender
# F    209
# M    340
sapply(cov1[,"Cardiogenic_shock"], table) 
# Cardiogenic_shock
# 0               520
# 1                29

# 输出表格
{
  library(tidyverse)
  library(psych)
  library(openxlsx)
  library(janitor)
  
  # 定义连续变量和分类变量
  concovname <- c("LOS_DAYS", "age", "DEATH_DAYS")
  catcovname <- c("hospital_expire_flag", "gender","race","shock","Cardiogenic_shock")
  
  # 定义分类变量的标签映射（根据你的数据描述进行调整）
  label_map <- list(
    hospital_expire_flag = c("0" = "Alive", "1" = "Dead"),
    gender = c("F" = "Female", "M" = "Male"),
    shock = c("0" = "Non-Shock", "1" = "Shock"),
    Cardiogenic_shock = c("0" = "Non-Cardiogenic_shock", "1" = "Cardiogenic_shock")
  )
  
  # 确保数据类型正确
  cov1 <- cov1 %>%
    mutate(across(all_of(concovname), as.numeric)) %>%
    mutate(across(all_of(catcovname), as.character)) %>%
    mutate(shock = as.character(shock))
  
  # --- 连续变量统计函数 ---
  cov_describe_custom <- function(df, vars, group_name) {
    desc_list <- lapply(vars, function(v) {
      vec <- df[[v]]
      if (is.null(vec) || all(is.na(vec))) {
        data.frame(n = 0, mean = NA, sd = NA, min = NA, median = NA, max = NA)
      } else {
        # 使用 psych::describe 并选择需要的列
        res <- psych::describe(vec, na.rm = TRUE)
        data.frame(n = res$n, mean = res$mean, sd = res$sd, 
                   min = res$min, median = res$median, max = res$max)
      }
    })
    
    desc <- do.call(rbind, desc_list)
    rownames(desc) <- vars
    
    desc %>%
      rownames_to_column(var = "Variable") %>%
      mutate(Group = as.character(group_name)) %>%
      select(Group, Variable, everything())
  }
  
  # --- 分类变量统计函数 ---
  cat_table_custom <- function(df, vars, group_name) {
    map_dfr(vars, ~ {
      var_name <- .x
      tab <- df %>%
        count(!!sym(var_name), name = "N") %>%
        filter(!is.na(!!sym(var_name))) %>% # 排除NA
        mutate(Percent = round(100 * N / sum(N), 2))
      
      # 添加 Total 行
      tab <- tab %>% adorn_totals("row")
      
      # 映射标签
      mapped_value <- as.character(tab[[var_name]])
      if (var_name %in% names(label_map)) {
        mapped_value <- recode(mapped_value, !!!label_map[[var_name]], .default = mapped_value)
      }
      
      tab %>%
        mutate(Value = mapped_value,
               Variable = var_name) %>%
        select(Variable, Value, N, Percent)
    }) %>%
      mutate(Group = as.character(group_name)) %>%
      select(Group, Variable, Value, N, Percent)
  }
  
  # --- A. 连续变量汇总 ---
  res_all_con <- cov_describe_custom(cov1, concovname, "ALL")
  
  # 按 shock 分组循环
  shock_levels <- sort(unique(cov1$shock))
  combined_con <- res_all_con
  
  for (s in shock_levels) {
    group_label <- recode(s, !!!label_map$shock)
    sub_df <- cov1 %>% filter(shock == s)
    temp <- cov_describe_custom(sub_df, concovname, group_label)
    combined_con <- bind_rows(combined_con, temp)
  }
  
  # --- B. 分类变量汇总 ---
  res_all_cat <- cat_table_custom(cov1, catcovname, "ALL")
  
  combined_cat <- res_all_cat
  for (s in shock_levels) {
    group_label <- recode(s, !!!label_map$shock)
    sub_df <- cov1 %>% filter(shock == s)
    temp <- cat_table_custom(sub_df, catcovname, group_label)
    combined_cat <- bind_rows(combined_cat, temp)
  }
  
  # 设置输出路径 (请根据你的电脑修改)
  output_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table\\filtered_Descriptive_Statistics_Shock_Groups.xlsx"
  
  wb <- createWorkbook()
  addWorksheet(wb, "Continuous_Vars")
  addWorksheet(wb, "Categorical_Vars")
  
  writeData(wb, "Continuous_Vars", combined_con)
  writeData(wb, "Categorical_Vars", combined_cat)
  
  saveWorkbook(wb, output_path, overwrite = TRUE)
  
  message("分析完成！结果已保存至: ", output_path)
  
}

{ 
  cov1 <- cov1 %>%
    mutate(race_group = case_when(
      grepl("WHITE|PORTUGUESE", race, ignore.case = TRUE) ~ "White",
      grepl("BLACK", race, ignore.case = TRUE) ~ "Black",
      grepl("ASIAN", race, ignore.case = TRUE) ~ "Asian",
      TRUE ~ "Other"  # 其余所有(Hispanic, Unknown, Native等)统一归为Other
    ))
  
  # 定义变量列表
  catcovname <- c("hospital_expire_flag", "gender", "race_group","Cardiogenic_shock")
  
  # 数值格式化
  fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "--", format(round(x, digits), nsmall = digits))
  }
  
  # 分类变量统计函数
  cat_table_calc <- function(df, vars, group_name) {
    df %>%
      select(all_of(vars)) %>%
      mutate(across(everything(), as.character)) %>%
      pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
      group_by(Variable, Value) %>%
      summarise(N = n(), .groups = "drop") %>%
      group_by(Variable) %>%
      mutate(Percent = N / sum(N) * 100, Group = group_name)
  }
  
  # 处理分类变量 
  res_all_cat <- cat_table_calc(cov1, catcovname, "ALL")
  res_s0_cat <- cat_table_calc(cov1 %>% filter(shock == 0), catcovname, "Non-Shock")
  res_s1_cat <- cat_table_calc(cov1 %>% filter(shock == 1), catcovname, "Shock")
  
  combined_cat <- bind_rows(res_all_cat, res_s0_cat, res_s1_cat)
  
  # 转换为宽表格式
  final_cat_wide <- combined_cat %>%
    mutate(`N (%)` = paste0(N, " (", fmt_num(Percent, 1), "%)")) %>%
    select(Group, Variable, Value, `N (%)`) %>%
    pivot_wider(names_from = Group, values_from = `N (%)`) %>%
    rename(Label = Value)
  
  # 处理连续变量 (假定 combined_con 已经存在)
  final_con_wide <- combined_con %>%
    mutate(
      `Mean (SD)` = paste0(fmt_num(mean), " (", fmt_num(sd), ")"),
      `Median [Min, Max]` = paste0(fmt_num(median), " [", fmt_num(min), ", ", fmt_num(max), "]")
    ) %>%
    select(Group, Variable, `Mean (SD)`, `Median [Min, Max]`) %>%
    pivot_longer(cols = c(`Mean (SD)`, `Median [Min, Max]`), names_to = "Statistic", values_to = "Value") %>%
    pivot_wider(names_from = Group, values_from = Value) %>%
    rename(Label = Statistic)
  
  # 组装
  res_list <- list()
  
  # --- 插入连续变量 ---
  for(v in unique(final_con_wide$Variable)) {
    sub <- final_con_wide %>% filter(Variable == v) %>% select(-Variable)
    header <- sub[1, ]; header[1, ] <- NA; header$Label <- v
    res_list[[paste0("con_", v)]] <- bind_rows(header, sub)
  }
  
  # --- 插入分类变量 (含种族特定排序) ---
  race_order <- c("White", "Black", "Asian", "Other")
  
  for(v in unique(final_cat_wide$Variable)) {
    sub <- final_cat_wide %>% filter(Variable == v) %>% select(-Variable)
    
    # 如果是种族，按指定顺序排序
    if(v == "race_group") {
      sub <- sub %>% 
        mutate(ord = match(Label, race_order)) %>% 
        arrange(ord) %>% 
        select(-ord)
    }
    
    header <- sub[1, ]; header[1, ] <- NA; header$Label <- v
    res_list[[paste0("cat_", v)]] <- bind_rows(header, sub)
  }
  
  # 垂直合并
  final_table_for_ppt <- bind_rows(res_list)
  
  output_dir <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  write.xlsx(final_table_for_ppt, file.path(output_dir, "filtered_PPT_Table_Race_Fixed.xlsx"), 
             na.string = "", overwrite = TRUE)
  
  message("表格已生成，种族已归类为 White, Black, Asian, Other。")
}


# 2. 实验室指标数据统计 #-----------------------------------------------------------------
# 统计每个 label 的测量人数
{
  label_counts <- all_data_rel %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  setwd(file.path(Output_derived_data, "data_filtered"))
  write.xlsx(label_counts,"filtred_mimic4_sepsis_实验室检查人数统计.xlsx")
  
  # 计算总人数及各指标覆盖率，提取缺失率 <= 30% 的标签
  total_n <- n_distinct(all_data_rel$hadm_id)
  
  keep_labels <- all_data_rel %>% group_by(label) %>% 
    summarise(n = n_distinct(hadm_id)) %>% filter(n / total_n >= 0.7) %>% pull(label)
  
  all_data_filtered <- all_data_rel %>% filter(label %in% keep_labels)
  
  label_stats <- all_data_rel %>% group_by(label) %>%
    summarise(missing_rate = 1 - (n_distinct(hadm_id) / total_n)) %>% arrange(missing_rate) 
}
df_lab <- all_data_filtered;
df_lab$valuenum <- as.numeric(as.character(df_lab$valuenum))
df_lab <- df_lab %>% 
  filter(!is.na(valuenum)) %>%  # 首先排除空值
  mutate(label = case_when(
    # 如果单位是%且数值在 20-60 之间，修正标签为 Hematocrit
    label == "Hemoglobin" & valueuom == "%" & valuenum > 20 ~ "Hematocrit",
    TRUE ~ as.character(label) )) %>% 
  filter(!(label == "Lactate" & (valuenum > 30 | valuenum < 0))) %>% 
  filter(!(label == "Glucose" & (valuenum > 1000 | valuenum < 0))) %>% 
  filter(!(label == "Base Excess" & (valuenum > 30 | valuenum < -30))) %>% 
  filter(!(label == "White Blood Cells" & (valuenum > 1000 | valuenum < 0))) %>% 
  filter(!(label == "Sodium" & (valuenum > 180 | valuenum < 90))) %>%
  filter(!(label == "Potassium" & (valuenum > 10 | valuenum < 1.0))) %>%
  filter(!(label == "Creatinine" & (valuenum > 20 | valuenum < 0.1))) %>% 
  # 2. 贫血与红细胞系列 (Hemoglobin, RBC, RDW)
  # 血红蛋白低于 3 就要输血，高于 25 极罕见 
  filter(!(label == "Red Blood Cells" & (valuenum > 12 | valuenum < 0.5))) %>%
  # RDW 是百分比，通常在 11-30 之间
  filter(!(label == "RDW" & (valuenum > 50 | valuenum < 5))) %>%
  
  # 3. 代谢与电解质系列 (Calcium, Total CO2)
  # 钙离子 8.5-10.5 为正常，低于 4 或高于 18 为生命极限
  filter(!(label == "Calcium, Total" & (valuenum > 20 | valuenum < 3))) %>%
  # Total CO2 反映碱储备，通常与 Bicarbonate 接近 (15-40)
  filter(!(label == "Calculated Total CO2" & (valuenum > 60 | valuenum < 2))) %>%
  
  # 4. 血气指标 (pCO2)
  # 二氧化碳分压 正常 35-45, 严重呼吸衰竭可能升至 100+, 但不会到 300
  filter(!(label == "pCO2" & (valuenum > 150 | valuenum < 5)))

log_needed_vars <- c(
  "NTproBNP", "proBNP, Pleural",
  "Asparate Aminotransferase (AST)", 
  "Alanine Aminotransferase (ALT)",
  "Creatine Kinase (CK)", 
  "Bilirubin, Total", "Bilirubin, Direct",
  "Troponin T","Lactate Dehydrogenase (LD)","Lactate" ,
  "Creatine Kinase, MB Isoenzyme"
)

df_lab_log <- df_lab %>%
  mutate(log_valuenum = if_else(label %in% log_needed_vars, log10(valuenum + 1), valuenum)) %>% 
  filter(!(label == "pH" & fluid != "Blood"))
setwd(file.path(Output_derived_data, "data_filtered"))
saveRDS(df_lab_log ,"filtered_mimic4_Sepsislog.rds") 

converted_labels <- df_lab_log %>%
  group_by(label) %>%
  summarise(is_log = any(log_valuenum != valuenum, na.rm = TRUE)) %>%
  filter(is_log == TRUE)
print(converted_labels$label)
# [1] "Alanine Aminotransferase (ALT)"  "Asparate Aminotransferase (AST)" "Bilirubin, Total"               
# [4] "Creatine Kinase (CK)"            "Creatine Kinase, MB Isoenzyme"   "Lactate"                        
# [7] "Lactate Dehydrogenase (LD)"      "Troponin T"  

# 3. 实验室指标统计表 #----------------------------------------------------
length(unique(df_lab_log$label))
unique(df_lab_log$label)


{
  library(dplyr)
  library(tidyr)
  library(broom) # 用于整洁地提取t检验结果
  
  # 定义生成包含 P 值的统计表函数
  generate_full_stat_table <- function(data) {
    
    # 1. 计算各组的基础统计量 (ALL, Non-Shock, Shock)
    # 计算分组统计 (0 和 1)
    group_stats <- data %>%
      group_by(label, shock) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = ifelse(shock == 1, "Shock", "Non-Shock"))
    
    # 计算全部统计 (ALL)
    all_stats <- data %>%
      group_by(label) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = "ALL")
    
    # 2. 核心步骤：对每个指标进行 t 检验计算 P 值
    p_values <- data %>%
      group_by(label) %>%
      filter(n_distinct(shock) == 2) %>% # 确保两组都有数据才能做检验
      do(tidy(t.test(valuenum ~ shock, data = .))) %>%
      select(label, p.value) %>%
      mutate(p_formatted = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ as.character(round(p.value, 3))
      ))
    
    # 3. 格式化统计字符串并合并 P 值
    final_stats <- bind_rows(group_stats, all_stats) %>%
      mutate(
        `Mean (SD)` = sprintf("%.2f (%.2f)", mean_val, sd_val),
        `Median [Min, Max]` = sprintf("%.2f [%.2f, %.2f]", med_val, min_val, max_val)
      ) %>%
      left_join(p_values, by = "label") %>%
      select(label, group_name, `Mean (SD)`, `Median [Min, Max]`, p_formatted)
    
    # 4. 转置表格为宽表格式
    table_wide <- final_stats %>%
      pivot_longer(cols = c(`Mean (SD)`, `Median [Min, Max]`), 
                   names_to = "Statistic", 
                   values_to = "Value") %>%
      pivot_wider(names_from = group_name, values_from = Value) %>%
      # 整理列顺序，并只在每个指标的第一行显示 P 值（美化）
      select(label, Statistic, ALL, `Non-Shock`, `Shock`, `P-value` = p_formatted) %>%
      mutate(`P-value` = ifelse(Statistic == "Median [Min, Max]", "", `P-value`)) # 第二行留白
    
    return(table_wide)
  }
  # --- 执行统计 ---
  biomarker_table <- generate_full_stat_table(df_lab_log)
  
  # 1. 创建分类与中英文对照映射表
  mapping_df <- data.frame(
    label = c(
      # --- 1. 心肌损伤与心脏功能指标 ---
      "Troponin T", "NTproBNP", "Creatine Kinase (CK)", "Creatine Kinase, MB Isoenzyme", "CK-MB Index", "proBNP, Pleural",
      # --- 2. 组织灌注与酸碱平衡指标 ---
      "Lactate", "Lactate Dehydrogenase (LD)", "Anion Gap", "Base Excess", "Bicarbonate", "Calculated Total CO2", 
      "pCO2", "pO2", "pH", "Oxygen Saturation", "Alveolar-arterial Gradient",
      # --- 3. 肝肾器官功能指标 ---
      "Alanine Aminotransferase (ALT)", "Asparate Aminotransferase (AST)", "Bilirubin, Direct", "Bilirubin, Indirect", 
      "Bilirubin, Total", "Albumin", "Creatinine", "Creatinine, Urine", "Urea Nitrogen", "Uric Acid",
      # --- 4. 血液学与凝血功能指标 ---
      "White Blood Cells", "Red Blood Cells", "Hemoglobin", "Hematocrit", "Platelet Count", "RDW", "RDW-SD", 
      "Neutrophils", "Lymphocytes", "Monocytes", "Immature Granulocytes", "Metamyelocytes", "Myelocytes", 
      "Bands", "Atypical Lymphocytes", "Nucleated Red Cells", "PT", "PTT", "INR(PT)", "Fibrinogen, Functional", "Antithrombin",
      # --- 5. 电解质与矿物质指标 ---
      "Potassium", "Potassium, Whole Blood", "Sodium", "Sodium, Whole Blood", "Chloride", "Chloride, Whole Blood", 
      "Calcium, Total", "Free Calcium", "Magnesium", "Phosphate",
      # --- 6. 炎症、代谢及其他指标 ---
      "C-Reactive Protein", "Glucose", "Osmolality, Measured", "Required O2", "Total Protein, Urine"
    ),
    label_cn = c(
      # 对应中文名称
      "Troponin T (肌钙蛋白T)", "NTproBNP (N端前B型利钠肽)", "Creatine Kinase (肌酸激酶)", "CK-MB Isoenzyme (肌酸激酶同工酶MB)", "CK-MB Index (CK-MB指数)", "proBNP, Pleural (胸水proBNP)",
      "Lactate (乳酸)", "Lactate Dehydrogenase (乳酸脱氢酶)", "Anion Gap (阴离子间隙)", "Base Excess (碱剩余)", "Bicarbonate (碳酸氢盐)", "Calculated Total CO2 (计算总二氧化碳)", 
      "pCO2 (二氧化碳分压)", "pO2 (氧分压)", "pH (酸碱度)", "Oxygen Saturation (氧饱和度)", "Alveolar-arterial Gradient (肺泡-动脉氧分压差)",
      "ALT (丙氨酸氨基转移酶)", "AST (天门冬氨酸氨基转移酶)", "Bilirubin, Direct (结合胆红素)", "Bilirubin, Indirect (非结合胆红素)", 
      "Bilirubin, Total (总胆红素)", "Albumin (白蛋白)", "Creatinine (肌酐)", "Creatinine, Urine (尿肌酐)", "Urea Nitrogen (尿素氮)", "Uric Acid (尿酸)",
      "White Blood Cells (白细胞计数)", "Red Blood Cells (红细胞计数)", "Hemoglobin (血红蛋白)", "Hematocrit (红细胞压积)", "Platelet Count (血小板计数)", "RDW (红细胞分布宽度)", "RDW-SD (红细胞分布宽度标准差)", 
      "Neutrophils (中性粒细胞)", "Lymphocytes (淋巴细胞)", "Monocytes (单核细胞)", "Immature Granulocytes (未成熟粒细胞)", "Metamyelocytes (晚幼粒细胞)", "Myelocytes (中幼粒细胞)", 
      "Bands (杆状核粒细胞)", "Atypical Lymphocytes (异型淋巴细胞)", "Nucleated Red Cells (有核红细胞)", "PT (凝血酶原时间)", "PTT (部分凝血活酶时间)", "INR(PT) (国际标准化比值)", "Fibrinogen (功能性纤维蛋白原)", "Antithrombin (抗凝血酶)",
      "Potassium (钾)", "Potassium, Whole Blood (全血钾)", "Sodium (钠)", "Sodium, Whole Blood (全血钠)", "Chloride (氯)", "Chloride, Whole Blood (全血氯)", 
      "Calcium, Total (总钙)", "Free Calcium (游离钙)", "Magnesium (镁)", "Phosphate (磷)",
      "C-Reactive Protein (C-反应蛋白)", "Glucose (血糖)", "Osmolality, Measured (实测渗透压)", "Required O2 (需氧量)", "Total Protein, Urine (尿总蛋白)"
    ),
    category = c(
      rep("1. Cardiac Injury & Function (心肌损伤与心脏功能指标)", 6),
      rep("2. Tissue Perfusion & Acid-Base Balance (组织灌注与酸碱平衡指标)", 11),
      rep("3. Hepatic & Renal Function (肝肾器官功能指标)", 10),
      rep("4. Hematology & Coagulation (血液学与凝血功能指标)", 21),
      rep("5. Electrolytes & Minerals (电解质与矿物质指标)", 10),
      rep("6. Inflammation, Metabolism & Others (炎症、代谢及其他指标)", 5)
    ),
    stringsAsFactors = FALSE
  )
  
  #  将映射合并到你生成的 biomarker_table 中
  final_ppt_table <- biomarker_table %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label) %>%
    select(category, label_cn, Statistic, ALL, `Non-Shock`, `Shock`, `P-value`)
  
  #  打印检查
  head(final_ppt_table)
  
  setwd("D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table")
  
  write.xlsx(final_ppt_table, "filtered_心肌标志物统计表_分类中文版.xlsx")
  
}

{
  library(dplyr)
  library(tidyr)
  library(broom)
  
  # 1. 计算各组的总人数 N (假设 subject_id 唯一，或直接按数据记录数计算)
  # 这里我们统计每组的观察数，用于表头显示
  n_all <- nrow(df_lab_log %>% distinct(hadm_id)) # 如果是按人统计
  n_no_shock <- nrow(df_lab_log %>% filter(shock == 0) %>% distinct(hadm_id))
  n_shock <- nrow(df_lab_log %>% filter(shock == 1) %>% distinct(hadm_id))
  
  # 如果数据是按记录算的，直接用 n_all <- nrow(df_lab_log) 等
  
  # 2. 定义格式化函数
  generate_ppt_table <- function(data) {
    
    # 获取所有指标列表
    labels <- unique(data$label)
    final_rows <- list()
    
    for (lab in labels) {
      # 提取单个指标数据
      sub_data <- data %>% filter(label == lab)
      
      # 计算 T 检验 P 值
      p_val <- t.test(valuenum ~ shock, data = sub_data)$p.value
      p_str <- if(p_val < 0.001) "P < 0.001" else sprintf("P = %.3f", p_val)
      
      # 计算三组统计量
      stats <- sub_data %>%
        group_by(shock) %>%
        summarise(
          mean_sd = sprintf("%.2f ( %.2f)", mean(valuenum, na.rm=T), sd(valuenum, na.rm=T)),
          med_range = sprintf("%.2f [ %.2f, %.2f]", median(valuenum, na.rm=T), min(valuenum, na.rm=T), max(valuenum, na.rm=T))
        )
      
      all_mean_sd <- sprintf("%.2f ( %.2f)", mean(sub_data$valuenum, na.rm=T), sd(sub_data$valuenum, na.rm=T))
      all_med_range <- sprintf("%.2f [ %.2f, %.2f]", median(sub_data$valuenum, na.rm=T), min(sub_data$valuenum, na.rm=T), max(sub_data$valuenum, na.rm=T))
      
      # 构建该指标的 4 行结构
      # 第一行：指标名称 (Label)
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      # 第四行：P 值
      
      df_block <- data.frame(
        Item = c(lab, "Mean (SD)", "Median [Min, Max]", p_str),
        All = c("", all_mean_sd, all_med_range, ""),
        No_shock = c("", stats$mean_sd[stats$shock == 0], stats$med_range[stats$shock == 0], ""),
        shock = c("", stats$mean_sd[stats$shock == 1], stats$med_range[stats$shock == 1], ""),
        stringsAsFactors = FALSE
      )
      
      final_rows[[lab]] <- df_block
    }
    
    # 合并所有指标块
    res_table <- bind_rows(final_rows)
    
    # 修改表头名称，加入 N 数值
    colnames(res_table) <- c(
      "", 
      paste0("All\n(N=", n_all, ")"), 
      paste0("No Cardiogenic Shock\n(N=", n_no_shock, ")"), 
      paste0("Cardiogenic Shock\n(N=", n_shock, ")")
    )
    
    return(res_table)
  }
  
  # --- 执行生成 ---
  ppt_stat_table <- generate_ppt_table(df_lab_log)
  
  # 查看结果
  print(ppt_stat_table)
}

{
  library(dplyr)
  library(openxlsx)
  
  # --- 1. 获取分组人数 N ---
  n_all <- 520
  n_shock <- 299
  n_no_shock <- n_all - n_shock
  
  # --- 2. 处理函数：将扁平表转换为块状 PPT 格式 ---
  # 请确保 final_ppt_table 包含列：category, label_cn, Statistic, ALL, Non-Shock, Shock, P-value
  df_input <- final_ppt_table 
  
  categories <- unique(df_input$category)
  final_list <- list()
  
  for (cat in categories) {
    
    # A. 插入分类大标题行 (例如：人口统计学特征)
    cat_header <- data.frame(
      Item = cat, 
      All = "", 
      No_Shock = "", 
      Shock = "", 
      stringsAsFactors = FALSE
    )
    final_list[[paste0(cat, "_header")]] <- cat_header
    
    # 获取该分类下的所有指标
    cat_data <- df_input %>% filter(category == cat)
    labels <- unique(cat_data$label_cn)
    
    for (lab in labels) {
      # 提取该指标的数据子集
      sub <- cat_data %>% filter(label_cn == lab)
      
      # 获取并格式化 P 值
      p_val_raw <- sub$`P-value`[1]
      p_suffix <- if(is.na(p_val_raw) || p_val_raw == "") "" else paste0(" (P ", p_val_raw, ")")
      
      # 组合 Label 和 P 值作为块的第一行
      display_label <- paste0(lab, p_suffix)
      
      # 构建 3 行块 (不再有单独的 P 值行)
      # 第一行：指标名称 + P值
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      block <- data.frame(
        Item = c(display_label, "  Mean (SD)", "  Median [Min, Max]"), # 前面加空格方便PPT对齐
        All = c("", 
                sub$ALL[sub$Statistic == "Mean (SD)"], 
                sub$ALL[sub$Statistic == "Median [Min, Max]"]),
        No_Shock = c("", 
                     sub$`Non-Shock`[sub$Statistic == "Mean (SD)"], 
                     sub$`Non-Shock`[sub$Statistic == "Median [Min, Max]"]),
        Shock = c("", 
                  sub$`Shock`[sub$Statistic == "Mean (SD)"], 
                  sub$`Shock`[sub$Statistic == "Median [Min, Max]"]),
        stringsAsFactors = FALSE
      )
      
      final_list[[paste0(cat, "_", lab)]] <- block
    }
  }
  
  # --- 3. 合并所有行并设置表头 ---
  final_table_ppt = bind_rows(final_list)
  
  # 设置表头，包含换行符 N 值
  colnames(final_table_ppt) <- c(
    "Characteristic", 
    paste0("All\n(N=", n_all, ")"), 
    paste0("No Shock\n(N=", n_no_shock, ")"), 
    paste0("Shock\n(N=", n_shock, ")")
  )
  
  # --- 4. 保存结果 ---
  output_file <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table\\PPT格式统计表_P值合并版.xlsx"
  write.xlsx(final_table_ppt, output_file)
  
  print(paste("转换完成！文件保存在：", output_file))
}

{
  label_counts <- df_lab_log %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  
  #  将映射合并到你生成的 biomarker_table 中
  final_label_counts <- label_counts %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label)  
  
  setwd("D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table")
  write.xlsx(final_label_counts,"中文分类filtered_mimic4_sepsis_实验室检查人数统计.xlsx")
}

# 休克组和非休克组实验室指标对比 #----------------------------------------------
unique(df_lab_log$label)
target_labelx1 <- c(
  # --- 1. 心肌损伤与心脏功能指标 ---
  "Troponin T", "NTproBNP", "Creatine Kinase (CK)", "Creatine Kinase, MB Isoenzyme", "CK-MB Index", "proBNP, Pleural"
)

target_labelx2 <- c(
  # --- 2. 组织灌注与酸碱平衡指标 ---
  "Lactate", "Lactate Dehydrogenase (LD)", "Anion Gap", "Base Excess", "Bicarbonate", "Calculated Total CO2", 
  "pCO2", "pO2", "pH", "Oxygen Saturation", "Alveolar-arterial Gradient"
)
target_labelx3 <- c( 
  # --- 3. 肝肾器官功能指标 ---
  "Alanine Aminotransferase (ALT)", "Asparate Aminotransferase (AST)", "Bilirubin, Direct", "Bilirubin, Indirect", 
  "Bilirubin, Total", "Albumin", "Creatinine", "Creatinine, Urine", "Urea Nitrogen", "Uric Acid" )
target_labelx4 <- c(
  # --- 4. 血液学与凝血功能指标 ---
  "White Blood Cells", "Red Blood Cells", "Hemoglobin", "Hematocrit", "Platelet Count", "RDW", "RDW-SD", 
  "Neutrophils", "Lymphocytes", "Monocytes", "Immature Granulocytes", "Metamyelocytes", "Myelocytes", 
  "Bands", "Atypical Lymphocytes", "Nucleated Red Cells", "PT", "PTT", "INR(PT)", "Fibrinogen, Functional", "Antithrombin"
)
target_labelx5<- c( 
  # --- 5. 电解质与矿物质指标 ---
  "Potassium", "Potassium, Whole Blood", "Sodium", "Sodium, Whole Blood", "Chloride", "Chloride, Whole Blood", 
  "Calcium, Total", "Free Calcium", "Magnesium", "Phosphate")
target_labelx6 <- c( 
  # --- 6. 炎症、代谢及其他指标 ---
  "C-Reactive Protein", "Glucose", "Osmolality, Measured", "Required O2", "Total Protein, Urine"
)

# 1. 配置参数与路径
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)

base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots"

# 1. 定义一个统一的 P 值格式化函数
format_p_value <- function(p) {
  if (is.na(p)) return("P = N/A")
  if (p < 0.001) return("P < 0.001")
  return(paste0("P = ", round(p, 3)))
}

# 2. 开始大循环
for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    
    message("正在绘制并计算双重显著性: ", target_biomarker)
    
    # --- A. 数据预处理 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        Group = factor(ifelse(shock == 1, "Shock Group", "Non-Shock Group"), 
                       levels = c("Non-Shock Group", "Shock Group")),
        # 为 Log 统计准备数据 (加极小值防止 log(0) 报错)
        val_log = log10(valuenum + 0.001),
        Phase = case_when(
          shock == 0 ~ "Non-Shock",
          charttime_rel < shock_onset_time_rel ~ "Pre-Shock",
          charttime_rel >= shock_onset_time_rel ~ "Post-Shock"
        )
      )
    
    if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) next
    
    #计算显著性 ---
    # 1. 原始数值的显著性 (Linear Scale P)
    p_raw_text <- "P = N/A"
    try({
      fit_raw <- lmer(valuenum ~ Group + charttime_rel + (1|subject_id), data = plot_df)
      p_raw <- summary(fit_raw)$coefficients["GroupShock Group", "Pr(>|t|)"]
      p_raw_text <- format_p_value(p_raw)
    }, silent = TRUE)
    
    # 2. 对数转换后的显著性 (Log-scale P)
    p_log_text <- "P = N/A"
    try({
      fit_log <- lmer(val_log ~ Group + charttime_rel + (1|subject_id), data = plot_df)
      p_log <- summary(fit_log)$coefficients["GroupShock Group", "Pr(>|t|)"]
      p_log_text <- format_p_value(p_log)
    }, silent = TRUE)
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # --- C. 绘图：线性坐标图 (标注原始 P) ---
    p1_linear <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(
        subtitle = paste0("Linear-scale diff: ", p_raw_text),
        x = "Days from Admission", y = paste0(target_biomarker, " (", unit, ")")
      ) +
      theme_cor + ggtitle("Linear Scale")
    
    # --- D. 绘图：对数坐标图 (标注对数 P) ---
    p1_log <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_y_log10() + 
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(
        subtitle = paste0("Log-scale diff: ", p_log_text), # 这里显示 Log P
        x = "Days from Admission", y = paste0("Log10 [ ", target_biomarker, " ]")
      ) +
      theme_cor + ggtitle("Log10 Scale")
    
    # --- E. 合并与保存 ---
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_Dual_P.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}
message("所有带显著性检验的轨迹图已保存。")

# 死亡组和非死亡组实验室指标对比 #----------------------------------------------

# 1. 辅助函数：格式化 P 值显示
format_p <- function(p) {
  if (is.na(p)) return("P = N/A")
  if (p < 0.001) return("P < 0.001")
  return(paste0("P = ", round(p, 3)))
}

# 2. 配置参数与新路径
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_Death_vs_Survival"

# 3. 开始大循环
for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    message("正在绘制死亡 vs 存活显著性: ", target_biomarker)
    
    # --- A. 数据预处理 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        # 【修改关键点：使用 hospital_expire_flag 分组】
        # 0 为存活 (Survivors), 1 为死亡 (Non-Survivors)
        Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                       levels = c("Survivors", "Non-Survivors")),
        
        # 这里的 Phase 逻辑保留，反映数据点在休克发生的前后关系（如果有的话）
        Phase = case_when(
          shock == 0 ~ "Non-Shock",
          charttime_rel < shock_onset_time_rel ~ "Pre-Shock",
          charttime_rel >= shock_onset_time_rel ~ "Post-Shock"
        )
      )
    
    # 健壮性检查
    if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
      message("跳过 ", target_biomarker, "：有效对比数据不足")
      next
    }
    
    # --- B. 计算显著性 (LMM) ---
    p_linear <- NA
    p_log_val <- NA
    
    # 线性模型显著性：Group 对数值的影响
    try({
      model_linear <- lmer(valuenum ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
      # 【提取关键：行名变为 GroupNon-Survivors】
      p_linear <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 对数模型显著性
    try({
      model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
      p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # --- C. 绘图部分 ---
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # 定义颜色：存活组蓝色，死亡组红色
    outcome_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
    
    # 线性坐标图
    p1_linear <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = outcome_colors) +
      scale_fill_manual(values = outcome_colors) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(x = "Days from Admission", y = paste0(target_biomarker, " (", unit, ")"),
           subtitle = paste0("Outcome Difference: ", format_p(p_linear))) +
      theme_cor + ggtitle("Linear Scale")
    
    # 对数坐标图
    p1_log <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = outcome_colors) +
      scale_fill_manual(values = outcome_colors) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      scale_y_log10() + 
      labs(x = "Days from Admission", y = paste0("Log10 [ ", target_biomarker, " ]"),
           subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
      ) +
      theme_cor + ggtitle("Log10 Scale")
    
    # 合并与保存
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_Death_vs_Survivors.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有死亡 vs 存活对比图已完成保存。") 
# 对齐休克为0时刻 休克组实验室指标ct图#---------------------------------

# 定义组名和对应的 label 列表
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)

base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_对齐休克0时刻"

for (group_name in names(label_groups)) {
  # group_name = "心肌损伤与心脏功能指标"
  # 创建子文件夹
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    # target_biomarker = "Troponin T" 
    message("正在绘制对齐图 [", group_name, "]: ", target_biomarker)
    
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        Group = ifelse(shock == 1, "Shock Group", "Non-Shock Group"),
        # 核心：时间对齐逻辑 (T=0 为休克发生点或入ICU点)
        time_aligned = ifelse(shock == 1, 
                              charttime_rel - shock_onset_time_rel, 
                              charttime_rel - icu_intime_first_rel),
        # 定义阶段
        Phase = case_when(
          shock == 0 ~ "Non-Shock",
          time_aligned < 0 ~ "Pre-Shock",
          time_aligned >= 0 ~ "Post-Shock"
        ),
        # 用于分段画趋势线
        LineGroup = interaction(Group, Phase)
      )  %>% 
      filter(shock == 1)
    
    # 健壮性检查
    if (nrow(plot_df) < 10) {
      message("跳过 ", target_biomarker, "：有效数据点太少")
      next
    }
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.3, size = 1, stroke = 0) + 
      # geom_line(aes(group = hadm_id),size = 0.8, alpha = 0.2) +
      # 关键：group = LineGroup 实现 0 点处的断开拟合
      geom_smooth(aes(group = LineGroup), method = "loess", size = 0.8, alpha = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 16)) +
      labs(x = "Days relative to Shock Onset (T=0)", y = paste0(target_biomarker, " (", unit, ")")) +
      theme_cor + ggtitle("Linear Scale Aligned")
    
    p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.3, size = 1, stroke = 0) + 
      # geom_line(aes(group = hadm_id),size = 0.8, alpha = 0.2) +
      geom_smooth(aes(group = LineGroup), method = "loess", size = 0.8, alpha = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      scale_y_log10() + 
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 16)) +
      labs(x = "Days relative to Shock Onset (T=0)", y = paste0("Log10 [ ", target_biomarker, " ]")) +
      theme_cor + ggtitle("Log10 Scale Aligned")
    
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "none")
    
    # 文件名清洗
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    ggsave(file.path(output_dir, paste0(safe_name, "_Aligned_T0.png")), 
           plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有 6 组对齐图片已处理完毕并按文件夹分类保存。")
colnames(plot_df)
# 对齐休克为0时刻 心源性休克组 vs 非心源性休克实验室指标ct图#---------------------------------
{
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_对齐休克为0时刻_CS_vs_SS"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>%
        mutate(
          # 重新定义对比组，并设置为 Factor（固定基准组为 Non-Cardio）
          Group = ifelse(Cardiogenic_shock == 1, "Cardiogenic Shock", "Non-Cardio Shock"),
          Group = factor(Group, levels = c("Non-Cardio Shock", "Cardiogenic Shock")),
          
          # 统一以休克时刻为 T=0
          time_aligned = charttime_rel - shock_onset_time_rel,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        )  
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足")
        next
      }
      
      # --- B. 计算显著性 (LMM) ---
      p_linear <- NA
      p_log_val <- NA
      
      # 线性模型显著性：Group 对数值的影响
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear <- summary(model_linear)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置
      group_colors <- c("Cardiogenic Shock" = "#377EB8", "Non-Cardio Shock" = "#E41A1C")
      
      # --- 4. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Outcome Difference: ", format_p(p_log_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Shock Type", fill = "Shock Type"
        ) +
        theme_cor +
        theme(legend.position = "none") 
      
      # --- 5. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]"),
          color = "Shock Type", fill = "Shock Type",
          subtitle = paste0("Outcome Difference: ", format_p(p_linear))
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 6. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_vs_NonCS.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("所有心源性 vs 非心源性对比图已保存。")
} 
# 对齐休克为0时刻：死亡组 vs 存活组 实验室指标ct图对比 #---------------------------------
{
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 修改输出根目录为 Outcome 比较
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_对齐休克为0时刻_Death_vs_Survival"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制结局对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        # 仅针对休克患者观察其结局差异
        filter(shock == 1) %>%
        mutate(
          
          Group = ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"),
          Group = factor(Group, levels = c("Survivors", "Non-Survivors")),
          
          # 统一以休克时刻为 T=0
          time_aligned = charttime_rel - shock_onset_time_rel,
          
          # 定义阶段
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          # 用于分段画趋势线
          LineGroup = interaction(Group, Phase)
        )  
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足")
        next
      }
      
      p_linear <- NA
      p_log_val <- NA
      
      # 线性模型显著性：Group 对数值的影响
      try({
        model_linear <- lmer(valuenum ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
        # 【提取关键：行名变为 GroupNon-Survivors】
        p_linear <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # --- 颜色设置：死亡红色，存活蓝色 ---
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- 1. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Outcome Difference: ", format_p(p_linear)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Hospital Outcome", fill = "Hospital Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 2. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]"),
          color = "Hospital Outcome", fill = "Hospital Outcome",
          subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 3. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_Outcome_Comparison.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("所有死亡 vs 存活对比图已完成。")
}
# 对齐休克为0时刻：非心源性休克：死亡组 vs 存活组 实验室指标ct图对比 #---------------------------------
{
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  
  # 辅助函数：格式化P值
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) return("N/A")
    if (p < 0.001) return("< 0.001")
    return(as.character(round(p, 3)))
  }
  
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 修改输出根目录
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_对齐休克为0时刻心源性休克_Death_vs_Survival"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制心源性休克结局对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 1. 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1 & Cardiogenic_shock == 1) %>% # 严格限定为心源性休克人群
        mutate(
          # 定义结局组，Survivors作为基准
          Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                         levels = c("Survivors", "Non-Survivors")),
          # 时间对齐 T=0
          time_aligned = charttime_rel - shock_onset_time_rel,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        ) 
      # filter(time_aligned >= -7 & time_aligned <= 7) # 观察核心窗口
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：有效数据不足或单一结局")
        next
      }
      
      # --- 2. 计算显著性 (LMM) ---
      p_linear_val <- NA
      p_log_val <- NA
      
      # 线性模型显著性
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear_val <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置：死亡红色，存活蓝色
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- 3. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Overall Diff P: ", format_p(p_linear_val)), # 修正这里
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Hospital Outcome", fill = "Hospital Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 4. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          subtitle = paste0("Overall Diff (Log) P: ", format_p(p_log_val)), # 修正这里
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]")
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 5. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      # 文件名处理
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_Outcome_T0.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("心源性休克死亡 vs 存活所有轨迹图已绘制完毕。")
}
# 休克发生前两组患者的实验室指标对比 #------------------------------------------
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_PreShock"

for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    
    message("正在绘制休克前对比图: ", target_biomarker)
    
    # --- 核心逻辑：剔除休克后的测量点 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      filter(
        # 条件1：如果是休克组，只取入院后、发病前的点
        (shock == 1 & charttime_rel >= 0 & charttime_rel < shock_onset_time_rel) |
          # 条件2：如果是非休克组，只取入院后的点
          (shock == 0 & charttime_rel >= 0)
      ) %>% 
      # # 限制观察时间范围，比如入院后前 7 天（防止长住院病人的尾部噪音）
      # filter(charttime_rel <= 7) %>%
      mutate(Group = factor(ifelse(shock == 1, "Future Shock", "Non-Shock"),
                            levels = c("Non-Shock", "Future Shock")))
    
    # 健壮性检查
    if (nrow(plot_df) < 10 | length(unique(plot_df$shock)) < 2) {
      message("跳过 ", target_biomarker, "：有效对比数据太少")
      next
    }
    
    # 线性模型显著性：Group 对数值的影响
    try({
      model_linear <- lmer(valuenum ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
      
      p_linear <- summary(model_linear)$coefficients["GroupFuture Shock", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 对数模型显著性
    try({
      model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel + (1 | subject_id), data = plot_df)
      p_log_val <- summary(model_log)$coefficients["GroupFuture Shock", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # --- 2. 线性坐标图 (仅含休克前点) ---
    p_linear <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(
        x = "Days from Admission",
        y = paste0(target_biomarker, " (", unit, ")"),
        subtitle = paste0("Outcome Difference: ", format_p(p_linear))
      ) +
      theme_cor + ggtitle("Linear Scale (Pre-Shock Only)")
    
    # --- 3. 对数坐标图 (仅含休克前点) ---
    p_log <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_y_log10() +
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(x = "Days from Admission", y = paste0("Log10[ ",target_biomarker, " (", unit, ") ]"),
           subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
      ) +
      theme_cor + ggtitle("Log10 Scale (Pre-Shock Only)")
    
    # --- 4. 合并与保存 ---
    combined_plot <- (p_linear | p_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_PreShock_Comparison.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有休克前对比图绘制完毕。")

# 休克发生后两组患者的实验室指标对比 #------------------------------------------
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)

# 建议存放在专门的 "Pre-Shock_Comparison" 文件夹
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\Biomarker_Plots_post_Shock"

for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    
    message("正在绘制休克后对比图: ", target_biomarker)
    
    # --- 核心逻辑：剔除休克后的测量点 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      filter(
        # 条件1：如果是休克组，只取入院后、发病后的点
        (shock == 1 & charttime_rel >= 0 & charttime_rel > shock_onset_time_rel) |
          # 条件2：如果是非休克组，只取入院后的点
          (shock == 0 & charttime_rel >= 0)
      ) %>% 
      # # 限制观察时间范围，比如入院后前 7 天（防止长住院病人的尾部噪音）
      # filter(charttime_rel <= 7) %>%
      mutate(Group = factor(ifelse(shock == 1, "Future Shock", "Non-Shock"),
                            levels = c("Non-Shock", "Future Shock")))
    
    # 健壮性检查
    if (nrow(plot_df) < 10 | length(unique(plot_df$shock)) < 2) {
      message("跳过 ", target_biomarker, "：有效对比数据太少")
      next
    }
    
    
    
    p_linear <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(
        x = "Days from Admission",
        y = paste0(target_biomarker, " (", unit, ")")
      ) +
      theme_cor + ggtitle("Linear Scale (Pre-Shock Only)")
    
    
    p_log <- ggplot(plot_df, aes(x = charttime_rel, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_y_log10() +
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(x = "Days from Admission", y = paste0("Log10[ ",target_biomarker, " (", unit, ") ]")) +
      theme_cor + ggtitle("Log10 Scale (Pre-Shock Only)")
    
    # --- 4. 合并与保存 ---
    combined_plot <- (p_linear | p_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_PreShock_Comparison.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有休克后对比图绘制完毕。")

# 三、第一次进入icu数据检视 #----------------------------------------------------------------
# 1. 提取首诊 ICU 期间的 labevents 记录
{
  df_final_paths <- df_final_paths %>% rename(first_intime = intime, first_outtime = outtime)
  colnames(df_final_paths)
  colnames(all_data_rel)
  df_first_icu_course <- all_data_rel %>%
    # 先确保实验室数据的 ID 和时间格式正确
    mutate(across(c(subject_id, hadm_id), as.numeric),
           charttime = as.POSIXct(charttime)) %>%
    
    # 关联第一步确定的“首诊窗口” (只根据 subject_id 关联)
    inner_join(df_final_paths %>% select(hadm_id, first_intime, first_outtime, is_shock_onset_in_this_icu,this_icu_outcome,endpoint_time,endpoint_state), 
               by = "hadm_id") %>%
    
    # 核心过滤：只保留化验时间落在首诊 ICU 进出时间范围内的记录
    # 我们通常会给 intime 往前加一个小的 buffer（如 1 小时），捕捉刚入室时的化验
    filter(charttime >= (first_intime - hours(1)) & charttime <= first_outtime)  
  
  cat("--- 提取完成 ---\n")
  cat("提取前化验记录总数: ", nrow(all_data_rel), "\n")
  cat(" 患者人数: ", n_distinct(df_first_icu_course$hadm_id), "\n")
  cat("首诊ICU内的化验记录总数: ", nrow(df_first_icu_course), "\n")
  
}
# 提取前化验记录总数:  464906 
# 患者人数:  548 
# 首诊ICU内的化验记录总数:  284621 

# 提取前化验记录总数:  355599 
# 患者人数:  548 
# 首诊ICU内的化验记录总数:  246365 
{
  
  first_icu_lab <- df_first_icu_course %>%
    mutate(
      # 1. 首先确保所有涉及计算的列均为 POSIXct 时间格式
      across(c(first_intime, charttime, endpoint_time, shock_onset_time, first_outtime), as.POSIXct),
      
      # 2. 计算相对时间（单位：天）
      # charttime_rel_icu: 实验室检查发生的时刻
      charttime_rel_icu = as.numeric(difftime(charttime, first_intime, units = "days")),
      
      # endpoint_rel_icu: 临床终点（死亡或出室）时刻
      endpoint_rel_icu = as.numeric(difftime(endpoint_time, first_intime, units = "days")),
      
      # shock_onset_rel_icu: 休克发生的时刻
      shock_onset_rel_icu = case_when(
        as.character(is_shock_onset_in_this_icu)  == "TRUE" ~ as.numeric(difftime(shock_onset_time, first_intime, units = "days")),
        TRUE ~ NA_real_ ),
      
      # outtime_rel_icu: 离开此次 ICU 的时刻
      outtime_rel_icu = as.numeric(difftime(first_outtime, first_intime, units = "days")),
      # deathtime_rel_icu: 此次 ICU内死亡 的时刻
      deathtime_rel_icu = case_when(
        this_icu_outcome == "Death in ICU" ~ as.numeric(difftime(deathtime, first_intime, units = "days")),
        TRUE ~ NA_real_ ) # 这里使用 NA_real_ 确保类型匹配（数值型）
    )
  unique(first_icu_lab$this_icu_outcome)
  # 检查是否存在负值（除了 charttime 允许有少量 buffer 负值外，其余通常应 >= 0）
  cat("相对时间计算完成。摘要统计：\n")
  summary(first_icu_lab %>% select(charttime_rel_icu, shock_onset_rel_icu, endpoint_rel_icu))
}
# charttime_rel_icu  shock_onset_rel_icu endpoint_rel_icu  
# Min.   :-0.04133   Min.   : 0.001      Min.   : 0.08198  
# 1st Qu.: 1.09514   1st Qu.: 0.245      1st Qu.: 6.04413  
# Median : 3.25347   Median : 0.505      Median :11.39175  
# Mean   : 5.13653   Mean   : 1.384      Mean   :12.53225  
# 3rd Qu.: 7.45347   3rd Qu.: 1.452      3rd Qu.:16.95755  
# Max.   :36.39053   Max.   :13.026      Max.   :36.96642  
#                    NA's   :68510  
# charttime_rel_icu  shock_onset_rel_icu endpoint_rel_icu  
# Min.   :-0.04133   Min.   : 0.001      Min.   : 0.08198  
# 1st Qu.: 1.14375   1st Qu.: 0.245      1st Qu.: 6.31554  
# Median : 3.36569   Median : 0.520      Median :11.61806  
# Mean   : 5.32466   Mean   : 1.396      Mean   :12.96411  
# 3rd Qu.: 7.74583   3rd Qu.: 1.527      3rd Qu.:17.22223  
# Max.   :36.39053   Max.   :13.026      Max.   :36.96642  
#                    NA's   :56328   
# 1. 提取首诊 ICU 期间的 Chartevents 记录
{
  df_first_icu_chartevent_course <- df_chartevents3 %>%
    filter(hadm_id %in% first_icu_lab$hadm_id) %>% 
    # 1. 确保 ID 为数值型，时间为 POSIXct
    mutate(across(c(subject_id, hadm_id), as.numeric),
           charttime = as.POSIXct(charttime)) %>%
    
    # 2. 关联首诊路径表   
    inner_join(df_final_paths %>% 
                 select(hadm_id, first_intime, first_outtime, 
                        is_shock_onset_in_this_icu, this_icu_outcome, 
                        endpoint_time, deathtime, shock_onset_time), 
               by = "hadm_id") %>%
    
    # 3. 核心过滤：只保留在首诊 ICU 时间窗内的监护记录
    # 同样给予 1 小时的 buffer 以捕捉入室瞬间的数据
    filter(charttime >= (first_intime - hours(1)) & charttime <= first_outtime)
  
  cat("--- Chartevents 提取完成 ---\n")
  cat("首诊患者人数: ", n_distinct(df_first_icu_chartevent_course$hadm_id), "\n")
  cat("提取出的监护记录总数: ", nrow(df_first_icu_chartevent_course), "\n")
  
  first_icu_chartevent <- df_first_icu_chartevent_course %>%
    mutate(
      # 确保涉及计算的列均为时间格式
      across(c(first_intime, charttime, endpoint_time, shock_onset_time, deathtime), as.POSIXct),
      
      # 计算相对时间（单位：days）
      # charttime_rel_icu: 监护记录发生的时刻
      charttime_rel_icu = as.numeric(difftime(charttime, first_intime, units = "days")),
      
      # endpoint_rel_icu: 临床终点时刻
      endpoint_rel_icu = as.numeric(difftime(endpoint_time, first_intime, units = "days")),
      
      # shock_onset_rel_icu: 休克发生的时刻
      shock_onset_rel_icu = case_when(
        as.character(is_shock_onset_in_this_icu) == "TRUE" ~ 
          as.numeric(difftime(shock_onset_time, first_intime, units = "days")),
        TRUE ~ NA_real_ 
      ),
      
      # outtime_rel_icu: 离开 ICU 时刻
      outtime_rel_icu = as.numeric(difftime(first_outtime, first_intime, units = "days")),
      
      # deathtime_rel_icu: ICU 内死亡时刻
      deathtime_rel_icu = case_when(
        this_icu_outcome == "Death in ICU" ~ 
          as.numeric(difftime(deathtime, first_intime, units = "days")),
        TRUE ~ NA_real_
      )
    )
  
  # 打印摘要检查
  summary(first_icu_chartevent %>% select(charttime_rel_icu, shock_onset_rel_icu, endpoint_rel_icu))
}
# --- Chartevents 提取完成 ---
#   首诊患者人数:  548 
# 提取出的监护记录总数:  357413 
# charttime_rel_icu  shock_onset_rel_icu endpoint_rel_icu  
# Min.   :-0.04167   Min.   : 0.001      Min.   : 0.08198  
# 1st Qu.: 2.11319   1st Qu.: 0.235      1st Qu.: 7.91190  
# Median : 4.77024   Median : 0.505      Median :12.60041  
# Mean   : 6.39413   Mean   : 1.370      Mean   :14.28008  
# 3rd Qu.: 9.15924   3rd Qu.: 1.527      3rd Qu.:18.41304  
# Max.   :36.72248   Max.   :13.026      Max.   :36.96642  
#                    NA's   :64461   
colnames(first_icu_chartevent)
colnames(first_icu_lab)


# 1.人口统计学数据和协变量分析-----------------------------------------------
colnames(first_icu_lab)
length(unique(first_icu_lab$hadm_id))
# [1] 548
cov1 <- first_icu_lab[!duplicated(paste(first_icu_lab$hadm_id)),]  #提取每个ID的第一行协变量值
colnames(cov1)
covname <- colnames(cov1[,c("outtime_rel_icu","age","endpoint_rel_icu","deathtime_rel_icu","shock_onset_rel_icu")])
catcovname <- colnames(cov1[,c("Cardiogenic_shock","this_icu_outcome","is_shock_onset_in_this_icu","Gender","race")])   #定义需要后续分析的连续协变量和分类协变量
colnames(cov1)
### 连续协变量描述性统计分析
summary(cov1[,covname])      #连续协变量描述性统计分析
# outtime_rel_icu         age        endpoint_rel_icu   deathtime_rel_icu shock_onset_rel_icu
# Min.   : 0.08198   Min.   :20.00   Min.   : 0.08198   Min.   : 0.2522   Min.   : 0.000926  
# 1st Qu.: 2.56578   1st Qu.:61.75   1st Qu.: 2.51765   1st Qu.: 3.5299   1st Qu.: 0.218559  
# Median : 5.23405   Median :71.00   Median : 5.21656   Median : 7.8069   Median : 0.532477  
# Mean   : 7.51317   Mean   :70.30   Mean   : 7.47846   Mean   : 9.4519   Mean   : 1.213554  
# 3rd Qu.:10.93084   3rd Qu.:82.00   3rd Qu.:10.88819   3rd Qu.:13.1566   3rd Qu.: 1.283194  
# Max.   :36.96642   Max.   :96.00   Max.   :36.96642   Max.   :29.1854   Max.   :13.026111  
#                                                       NA's   :419       NA's   :249        

quantile(cov1$outtime_rel_icu , probs = 0.95, na.rm = TRUE)
# 95% 
# 20.56503 
quantile(cov1$deathtime_rel_icu, probs = 0.05,  na.rm = TRUE)
# 5% 
# 2.162134 
table(cov1$this_icu_outcome) 
# Death in ICU Discharged Alive from ICU 
# 129                       419 
table(cov1$this_icu_outcome[which(cov1$Cardiogenic_shock==1)]) # 心源性休克患者中死亡人数统计
# Death in ICU Discharged Alive from ICU 
# 9                        20 

### 分类协变量描述性统计分析
sapply(cov1[catcovname], table)  
sapply(cov1[,"gender"], table)  
# gender
# F    209
# M    339
sapply(cov1[,"Cardiogenic_shock"], table) 
# Cardiogenic_shock
# 0               519
# 1                29
sapply(cov1[,"is_shock_onset_in_this_icu"], table) 
# is_shock_onset_in_this_icu
# FALSE                        249
# TRUE                         299

# 输出表格
{
  library(tidyverse)
  library(psych)
  library(openxlsx)
  library(janitor)
  
  # --- 1. 变量定义 ---
  # 连续变量
  concovname <- c("age", "outtime_rel_icu", "endpoint_rel_icu", "deathtime_rel_icu")
  # 分类变量 (注意：race_group 是稍后生成的)
  catcovname <- c("this_icu_outcome", "gender", "race_group", "Cardiogenic_shock")
  
  # --- 2. 数据清洗与 0/1 转换 ---
  cov1_clean <- cov1 %>%
    # A. 种族四分类
    mutate(race_group = case_when(
      grepl("WHITE|PORTUGUESE", race, ignore.case = TRUE) ~ "White",
      grepl("BLACK", race, ignore.case = TRUE) ~ "Black",
      grepl("ASIAN", race, ignore.case = TRUE) ~ "Asian",
      TRUE ~ "Other"
    )) %>%
    # B. 将 is_shock_onset_in_this_icu 转换为数值 0 和 1
    mutate(is_shock_onset_in_this_icu = case_when(
      as.character(is_shock_onset_in_this_icu) %in% c("TRUE", "1") ~ 1,
      as.character(is_shock_onset_in_this_icu) %in% c("FALSE", "0") ~ 0,
      TRUE ~ NA_real_
    )) %>%
    # C. 统一数据类型
    mutate(across(all_of(concovname), as.numeric)) %>%
    mutate(across(all_of(catcovname), as.character))
  
  # --- 3. 定义标签映射与辅助函数 ---
  label_map <- list(
    this_icu_outcome = c("Death in ICU" = "Dead", "Discharged Alive from ICU" = "Alive"),
    gender = c("F" = "Female", "M" = "Male"),
    Cardiogenic_shock = c("0" = "No", "1" = "Yes")
  )
  
  fmt_num <- function(x, digits = 2) {
    ifelse(is.na(x), "--", format(round(as.numeric(x), digits), nsmall = digits))
  }
  
  # --- 4. 统计函数定义 ---
  
  # 连续变量：输出 Mean (SD) 和 Median [Range]
  cov_describe_custom <- function(df, vars, group_name) {
    map_dfr(vars, function(v) {
      res <- psych::describe(df[[v]], na.rm = TRUE)
      data.frame(
        Group = group_name,
        Variable = v,
        "Mean (SD)" = paste0(fmt_num(res$mean), " (", fmt_num(res$sd), ")"),
        "Median [Min, Max]" = paste0(fmt_num(res$median), " [", fmt_num(res$min), ", ", fmt_num(res$max), "]"),
        check.names = FALSE
      )
    })
  }
  
  # 分类变量：输出 N (%)
  cat_table_custom <- function(df, vars, group_name) {
    map_dfr(vars, function(v) {
      tab <- df %>%
        count(!!sym(v), name = "N") %>%
        filter(!is.na(!!sym(v))) %>%
        mutate(Percent = N / sum(N) * 100)
      
      val_labels <- as.character(tab[[v]])
      if (v %in% names(label_map)) {
        val_labels <- recode(val_labels, !!!label_map[[v]], .default = val_labels)
      }
      
      data.frame(
        Group = group_name,
        Variable = v,
        Label = val_labels,
        "N (%)" = paste0(tab$N, " (", fmt_num(tab$Percent, 1), "%)"),
        check.names = FALSE
      )
    })
  }
  
  # --- 5. 执行分组统计 (基于 0/1 分组) ---
  
  # A. 连续变量汇总
  combined_con <- bind_rows(
    cov_describe_custom(cov1_clean, concovname, "ALL"),
    cov_describe_custom(cov1_clean %>% filter(is_shock_onset_in_this_icu == 0), concovname, "Non-Shock (0)"),
    cov_describe_custom(cov1_clean %>% filter(is_shock_onset_in_this_icu == 1), concovname, "Shock (1)")
  ) %>%
    pivot_longer(cols = c("Mean (SD)", "Median [Min, Max]"), names_to = "Statistic", values_to = "Value") %>%
    pivot_wider(names_from = Group, values_from = Value) %>%
    rename(Label = Statistic)
  
  # B. 分类变量汇总
  combined_cat <- bind_rows(
    cat_table_custom(cov1_clean, catcovname, "ALL"),
    cat_table_custom(cov1_clean %>% filter(is_shock_onset_in_this_icu == 0), catcovname, "Non-Shock (0)"),
    cat_table_custom(cov1_clean %>% filter(is_shock_onset_in_this_icu == 1), catcovname, "Shock (1)")
  ) %>%
    pivot_wider(names_from = Group, values_from = "N (%)")
  
  # --- 6. 组装 PPT 格式最终表 (插入标题行 & 种族排序) ---
  res_list <- list()
  
  # 插入连续变量块
  for(v in unique(combined_con$Variable)) {
    sub <- combined_con %>% filter(Variable == v) %>% select(-Variable)
    header <- sub[1, ]; header[1, ] <- NA; header$Label <- v
    res_list[[paste0("con_", v)]] <- bind_rows(header, sub)
  }
  
  # 插入分类变量块 (含种族特定排序)
  race_order <- c("White", "Black", "Asian", "Other")
  for(v in unique(combined_cat$Variable)) {
    sub <- combined_cat %>% filter(Variable == v) %>% select(-Variable)
    
    if(v == "race_group") {
      sub <- sub %>% mutate(ord = match(Label, race_order)) %>% arrange(ord) %>% select(-ord)
    }
    
    header <- sub[1, ]; header[1, ] <- NA; header$Label <- v
    res_list[[paste0("cat_", v)]] <- bind_rows(header, sub)
  }
  
  final_ppt_table <- bind_rows(res_list)
  
  # --- 7. 导出 Excel ---
  output_dir <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table"
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  file_name <- "Table1_ICU_Shock_Comparison_POP.xlsx"
  write.xlsx(final_ppt_table, file.path(output_dir, file_name), na.string = "", overwrite = TRUE)
  
  cat("\n==========================================\n")
  cat("✅ 处理完成！\n")
  cat("is_shock_onset_in_this_icu 已转换为 0/1 \n")
  cat("表格已保存至: ", file.path(output_dir, file_name), "\n")
  cat("==========================================\n")
}


# 2. 实验室指标数据统计 #-----------------------------------------------------------------
# 统计每个 label 的测量人数
{
  first_icu_lab1 <- first_icu_lab %>%  mutate(
    orignal_shock = shock,
    shock = is_shock_onset_in_this_icu,
    original_hospital_expire_flag = hospital_expire_flag,
    hospital_expire_flag = case_when(this_icu_outcome == "Death in ICU" ~ 1,
                                     TRUE ~ 0 
    )
  )
  
  label_counts <- first_icu_lab1 %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  
  output_dir <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table"
  file_name <- "Table2.ICU实验室检查人数统计.xlsx"
  write.xlsx(label_counts, file.path(output_dir, file_name))
  
  
  # 计算总人数及各指标覆盖率，提取缺失率 <= 30% 的标签
  total_n <- n_distinct(first_icu_lab1$hadm_id)
  
  keep_labels <- first_icu_lab1 %>% group_by(label) %>% 
    summarise(n = n_distinct(hadm_id)) %>% filter(n / total_n >= 0.7) %>% pull(label)
  
  all_data_filtered <- first_icu_lab1 %>% filter(label %in% keep_labels)
  
  label_stats <- first_icu_lab1 %>% group_by(label) %>%
    summarise(missing_rate = 1 - (n_distinct(hadm_id) / total_n)) %>% arrange(missing_rate) 
}
df_lab <- all_data_filtered
df_lab$valuenum <- as.numeric(as.character(df_lab$valuenum))
df_lab <- df_lab %>% 
  filter(!is.na(valuenum)) %>%  # 首先排除空值
  mutate(label = case_when(
    # 如果单位是%且数值在 20-60 之间，修正标签为 Hematocrit
    label == "Hemoglobin" & valueuom == "%" & valuenum > 20 ~ "Hematocrit",
    TRUE ~ as.character(label) )) %>% 
  filter(!(label == "Lactate" & (valuenum > 30 | valuenum < 0))) %>% 
  filter(!(label == "Glucose" & (valuenum > 1000 | valuenum < 0))) %>% 
  filter(!(label == "Base Excess" & (valuenum > 30 | valuenum < -30))) %>% 
  filter(!(label == "White Blood Cells" & (valuenum > 1000 | valuenum < 0))) %>% 
  filter(!(label == "Sodium" & (valuenum > 180 | valuenum < 90))) %>%
  filter(!(label == "Potassium" & (valuenum > 10 | valuenum < 1.0))) %>%
  filter(!(label == "Creatinine" & (valuenum > 20 | valuenum < 0.1))) %>% 
  # 2. 贫血与红细胞系列 (Hemoglobin, RBC, RDW)
  # 血红蛋白低于 3 就要输血，高于 25 极罕见 
  filter(!(label == "Red Blood Cells" & (valuenum > 12 | valuenum < 0.5))) %>%
  # RDW 是百分比，通常在 11-30 之间
  filter(!(label == "RDW" & (valuenum > 50 | valuenum < 5))) %>%
  
  # 3. 代谢与电解质系列 (Calcium, Total CO2)
  # 钙离子 8.5-10.5 为正常，低于 4 或高于 18 为生命极限
  filter(!(label == "Calcium, Total" & (valuenum > 20 | valuenum < 3))) %>%
  # Total CO2 反映碱储备，通常与 Bicarbonate 接近 (15-40)
  filter(!(label == "Calculated Total CO2" & (valuenum > 60 | valuenum < 2))) %>%
  
  # 4. 血气指标 (pCO2)
  # 二氧化碳分压 正常 35-45, 严重呼吸衰竭可能升至 100+, 但不会到 300
  filter(!(label == "pCO2" & (valuenum > 150 | valuenum < 5)))

log_needed_vars <- c(
  "NTproBNP", "proBNP, Pleural",
  "Asparate Aminotransferase (AST)", 
  "Alanine Aminotransferase (ALT)",
  "Creatine Kinase (CK)", 
  "Bilirubin, Total", "Bilirubin, Direct",
  "Troponin T","Lactate Dehydrogenase (LD)","Lactate" ,
  "Creatine Kinase, MB Isoenzyme"
)

df_lab_log <- df_lab %>%
  mutate(log_valuenum = if_else(label %in% log_needed_vars, log10(valuenum + 1), valuenum)) %>% 
  filter(!(label == "pH" & fluid != "Blood"))
setwd(file.path(Output_derived_data, "data_filtered"))
saveRDS(df_lab_log ,"icu_filtered_mimic4_Sepsislog.rds") 

converted_labels <- df_lab_log %>%
  group_by(label) %>%
  summarise(is_log = any(log_valuenum != valuenum, na.rm = TRUE)) %>%
  filter(is_log == TRUE)
print(converted_labels$label)
# [1] "Alanine Aminotransferase (ALT)"  "Asparate Aminotransferase (AST)" "Bilirubin, Total"               
# [4] "Creatine Kinase (CK)"            "Creatine Kinase, MB Isoenzyme"   "Lactate"                        
# [7] "Lactate Dehydrogenase (LD)"      "Troponin T"  

# 3. 实验室指标统计表 #----------------------------------------------------
setwd(file.path(Output_derived_data, "data_filtered"))
df_lab_log <- readRDS("icu_filtered_mimic4_Sepsislog.rds") 

length(unique(df_lab_log$label))
unique(df_lab_log$label)
length(unique(df_lab_log$hadm_id))
length(unique(df_lab_log$hadm_id[df_lab_log$shock == 1]))
{
  library(dplyr)
  library(tidyr)
  library(broom) # 用于整洁地提取t检验结果
  
  # 定义生成包含 P 值的统计表函数
  generate_full_stat_table <- function(data) {
    
    # 1. 计算各组的基础统计量 (ALL, Non-Shock, Shock)
    # 计算分组统计 (0 和 1)
    group_stats <- data %>%
      group_by(label, shock) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = ifelse(shock == 1, "Shock", "Non-Shock"))
    
    # 计算全部统计 (ALL)
    all_stats <- data %>%
      group_by(label) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = "ALL")
    
    # 2. 核心步骤：对每个指标进行 t 检验计算 P 值
    p_values <- data %>%
      group_by(label) %>%
      filter(n_distinct(shock) == 2) %>% # 确保两组都有数据才能做检验
      do(tidy(t.test(valuenum ~ shock, data = .))) %>%
      select(label, p.value) %>%
      mutate(p_formatted = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ as.character(round(p.value, 3))
      ))
    
    # 3. 格式化统计字符串并合并 P 值
    final_stats <- bind_rows(group_stats, all_stats) %>%
      mutate(
        `Mean (SD)` = sprintf("%.2f (%.2f)", mean_val, sd_val),
        `Median [Min, Max]` = sprintf("%.2f [%.2f, %.2f]", med_val, min_val, max_val)
      ) %>%
      left_join(p_values, by = "label") %>%
      select(label, group_name, `Mean (SD)`, `Median [Min, Max]`, p_formatted)
    
    # 4. 转置表格为宽表格式
    table_wide <- final_stats %>%
      pivot_longer(cols = c(`Mean (SD)`, `Median [Min, Max]`), 
                   names_to = "Statistic", 
                   values_to = "Value") %>%
      pivot_wider(names_from = group_name, values_from = Value) %>%
      # 整理列顺序，并只在每个指标的第一行显示 P 值（美化）
      select(label, Statistic, ALL, `Non-Shock`, `Shock`, `P-value` = p_formatted) %>%
      mutate(`P-value` = ifelse(Statistic == "Median [Min, Max]", "", `P-value`)) # 第二行留白
    
    return(table_wide)
  }
  # --- 执行统计 ---
  biomarker_table <- generate_full_stat_table(df_lab_log)
  
  # 1. 创建分类与中英文对照映射表
  mapping_df <- data.frame(
    label = c(
      # --- 1. 心肌损伤与心脏功能指标 ---
      "Troponin T", "NTproBNP", "Creatine Kinase (CK)", "Creatine Kinase, MB Isoenzyme", "CK-MB Index", "proBNP, Pleural",
      # --- 2. 组织灌注与酸碱平衡指标 ---
      "Lactate", "Lactate Dehydrogenase (LD)", "Anion Gap", "Base Excess", "Bicarbonate", "Calculated Total CO2", 
      "pCO2", "pO2", "pH", "Oxygen Saturation", "Alveolar-arterial Gradient",
      # --- 3. 肝肾器官功能指标 ---
      "Alanine Aminotransferase (ALT)", "Asparate Aminotransferase (AST)", "Bilirubin, Direct", "Bilirubin, Indirect", 
      "Bilirubin, Total", "Albumin", "Creatinine", "Creatinine, Urine", "Urea Nitrogen", "Uric Acid",
      # --- 4. 血液学与凝血功能指标 ---
      "White Blood Cells", "Red Blood Cells", "Hemoglobin", "Hematocrit", "Platelet Count", "RDW", "RDW-SD", 
      "Neutrophils", "Lymphocytes", "Monocytes", "Immature Granulocytes", "Metamyelocytes", "Myelocytes", 
      "Bands", "Atypical Lymphocytes", "Nucleated Red Cells", "PT", "PTT", "INR(PT)", "Fibrinogen, Functional", "Antithrombin",
      # --- 5. 电解质与矿物质指标 ---
      "Potassium", "Potassium, Whole Blood", "Sodium", "Sodium, Whole Blood", "Chloride", "Chloride, Whole Blood", 
      "Calcium, Total", "Free Calcium", "Magnesium", "Phosphate",
      # --- 6. 炎症、代谢及其他指标 ---
      "C-Reactive Protein", "Glucose", "Osmolality, Measured", "Required O2", "Total Protein, Urine"
    ),
    label_cn = c(
      # 对应中文名称
      "Troponin T (肌钙蛋白T)", "NTproBNP (N端前B型利钠肽)", "Creatine Kinase (肌酸激酶)", "CK-MB Isoenzyme (肌酸激酶同工酶MB)", "CK-MB Index (CK-MB指数)", "proBNP, Pleural (胸水proBNP)",
      "Lactate (乳酸)", "Lactate Dehydrogenase (乳酸脱氢酶)", "Anion Gap (阴离子间隙)", "Base Excess (碱剩余)", "Bicarbonate (碳酸氢盐)", "Calculated Total CO2 (计算总二氧化碳)", 
      "pCO2 (二氧化碳分压)", "pO2 (氧分压)", "pH (酸碱度)", "Oxygen Saturation (氧饱和度)", "Alveolar-arterial Gradient (肺泡-动脉氧分压差)",
      "ALT (丙氨酸氨基转移酶)", "AST (天门冬氨酸氨基转移酶)", "Bilirubin, Direct (结合胆红素)", "Bilirubin, Indirect (非结合胆红素)", 
      "Bilirubin, Total (总胆红素)", "Albumin (白蛋白)", "Creatinine (肌酐)", "Creatinine, Urine (尿肌酐)", "Urea Nitrogen (尿素氮)", "Uric Acid (尿酸)",
      "White Blood Cells (白细胞计数)", "Red Blood Cells (红细胞计数)", "Hemoglobin (血红蛋白)", "Hematocrit (红细胞压积)", "Platelet Count (血小板计数)", "RDW (红细胞分布宽度)", "RDW-SD (红细胞分布宽度标准差)", 
      "Neutrophils (中性粒细胞)", "Lymphocytes (淋巴细胞)", "Monocytes (单核细胞)", "Immature Granulocytes (未成熟粒细胞)", "Metamyelocytes (晚幼粒细胞)", "Myelocytes (中幼粒细胞)", 
      "Bands (杆状核粒细胞)", "Atypical Lymphocytes (异型淋巴细胞)", "Nucleated Red Cells (有核红细胞)", "PT (凝血酶原时间)", "PTT (部分凝血活酶时间)", "INR(PT) (国际标准化比值)", "Fibrinogen (功能性纤维蛋白原)", "Antithrombin (抗凝血酶)",
      "Potassium (钾)", "Potassium, Whole Blood (全血钾)", "Sodium (钠)", "Sodium, Whole Blood (全血钠)", "Chloride (氯)", "Chloride, Whole Blood (全血氯)", 
      "Calcium, Total (总钙)", "Free Calcium (游离钙)", "Magnesium (镁)", "Phosphate (磷)",
      "C-Reactive Protein (C-反应蛋白)", "Glucose (血糖)", "Osmolality, Measured (实测渗透压)", "Required O2 (需氧量)", "Total Protein, Urine (尿总蛋白)"
    ),
    category = c(
      rep("1. Cardiac Injury & Function (心肌损伤与心脏功能指标)", 6),
      rep("2. Tissue Perfusion & Acid-Base Balance (组织灌注与酸碱平衡指标)", 11),
      rep("3. Hepatic & Renal Function (肝肾器官功能指标)", 10),
      rep("4. Hematology & Coagulation (血液学与凝血功能指标)", 21),
      rep("5. Electrolytes & Minerals (电解质与矿物质指标)", 10),
      rep("6. Inflammation, Metabolism & Others (炎症、代谢及其他指标)", 5)
    ),
    stringsAsFactors = FALSE
  )
  
  #  将映射合并到你生成的 biomarker_table 中
  final_ppt_table <- biomarker_table %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label) %>%
    select(category, label_cn, Statistic, ALL, `Non-Shock`, `Shock`, `P-value`)
  
  #  打印检查
  head(final_ppt_table)
  
  setwd("D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table")
  
  write.xlsx(final_ppt_table, "icu_filtered_心肌标志物统计表_分类中文版.xlsx")
  
}

{
  library(dplyr)
  library(tidyr)
  library(broom)
  
  # 1. 计算各组的总人数 N (假设 subject_id 唯一，或直接按数据记录数计算)
  # 这里我们统计每组的观察数，用于表头显示
  n_all <- nrow(df_lab_log %>% distinct(hadm_id)) # 如果是按人统计
  n_no_shock <- nrow(df_lab_log %>% filter(shock == 0) %>% distinct(hadm_id))
  n_shock <- nrow(df_lab_log %>% filter(shock == 1) %>% distinct(hadm_id))
  
  # 如果数据是按记录算的，直接用 n_all <- nrow(df_lab_log) 等
  
  # 2. 定义格式化函数
  generate_ppt_table <- function(data) {
    
    # 获取所有指标列表
    labels <- unique(data$label)
    final_rows <- list()
    
    for (lab in labels) {
      # 提取单个指标数据
      sub_data <- data %>% filter(label == lab)
      
      # 计算 T 检验 P 值
      p_val <- t.test(valuenum ~ shock, data = sub_data)$p.value
      p_str <- if(p_val < 0.001) "P < 0.001" else sprintf("P = %.3f", p_val)
      
      # 计算三组统计量
      stats <- sub_data %>%
        group_by(shock) %>%
        summarise(
          mean_sd = sprintf("%.2f ( %.2f)", mean(valuenum, na.rm=T), sd(valuenum, na.rm=T)),
          med_range = sprintf("%.2f [ %.2f, %.2f]", median(valuenum, na.rm=T), min(valuenum, na.rm=T), max(valuenum, na.rm=T))
        )
      
      all_mean_sd <- sprintf("%.2f ( %.2f)", mean(sub_data$valuenum, na.rm=T), sd(sub_data$valuenum, na.rm=T))
      all_med_range <- sprintf("%.2f [ %.2f, %.2f]", median(sub_data$valuenum, na.rm=T), min(sub_data$valuenum, na.rm=T), max(sub_data$valuenum, na.rm=T))
      
      # 构建该指标的 4 行结构
      # 第一行：指标名称 (Label)
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      # 第四行：P 值
      
      df_block <- data.frame(
        Item = c(lab, "Mean (SD)", "Median [Min, Max]", p_str),
        All = c("", all_mean_sd, all_med_range, ""),
        No_shock = c("", stats$mean_sd[stats$shock == 0], stats$med_range[stats$shock == 0], ""),
        shock = c("", stats$mean_sd[stats$shock == 1], stats$med_range[stats$shock == 1], ""),
        stringsAsFactors = FALSE
      )
      
      final_rows[[lab]] <- df_block
    }
    
    # 合并所有指标块
    res_table <- bind_rows(final_rows)
    
    # 修改表头名称，加入 N 数值
    colnames(res_table) <- c(
      "", 
      paste0("All\n(N=", n_all, ")"), 
      paste0("No Cardiogenic Shock\n(N=", n_no_shock, ")"), 
      paste0("Cardiogenic Shock\n(N=", n_shock, ")")
    )
    
    return(res_table)
  }
  
  # --- 执行生成 ---
  ppt_stat_table <- generate_ppt_table(df_lab_log)
  
  # 查看结果
  print(ppt_stat_table)
}

{
  library(dplyr)
  library(openxlsx)
  
  # --- 1. 获取分组人数 N ---
  n_all <- 548
  n_shock <- 299
  n_no_shock <- n_all - n_shock
  
  # --- 2. 处理函数：将扁平表转换为块状 PPT 格式 ---
  # 请确保 final_ppt_table 包含列：category, label_cn, Statistic, ALL, Non-Shock, Shock, P-value
  df_input <- final_ppt_table 
  
  categories <- unique(df_input$category)
  final_list <- list()
  
  for (cat in categories) {
    
    # A. 插入分类大标题行 (例如：人口统计学特征)
    cat_header <- data.frame(
      Item = cat, 
      All = "", 
      No_Shock = "", 
      Shock = "", 
      stringsAsFactors = FALSE
    )
    final_list[[paste0(cat, "_header")]] <- cat_header
    
    # 获取该分类下的所有指标
    cat_data <- df_input %>% filter(category == cat)
    labels <- unique(cat_data$label_cn)
    
    for (lab in labels) {
      # 提取该指标的数据子集
      sub <- cat_data %>% filter(label_cn == lab)
      
      # 获取并格式化 P 值
      p_val_raw <- sub$`P-value`[1]
      p_suffix <- if(is.na(p_val_raw) || p_val_raw == "") "" else paste0(" (P ", p_val_raw, ")")
      
      # 组合 Label 和 P 值作为块的第一行
      display_label <- paste0(lab, p_suffix)
      
      # 构建 3 行块 (不再有单独的 P 值行)
      # 第一行：指标名称 + P值
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      block <- data.frame(
        Item = c(display_label, "  Mean (SD)", "  Median [Min, Max]"), # 前面加空格方便PPT对齐
        All = c("", 
                sub$ALL[sub$Statistic == "Mean (SD)"], 
                sub$ALL[sub$Statistic == "Median [Min, Max]"]),
        No_Shock = c("", 
                     sub$`Non-Shock`[sub$Statistic == "Mean (SD)"], 
                     sub$`Non-Shock`[sub$Statistic == "Median [Min, Max]"]),
        Shock = c("", 
                  sub$`Shock`[sub$Statistic == "Mean (SD)"], 
                  sub$`Shock`[sub$Statistic == "Median [Min, Max]"]),
        stringsAsFactors = FALSE
      )
      
      final_list[[paste0(cat, "_", lab)]] <- block
    }
  }
  
  # --- 3. 合并所有行并设置表头 ---
  final_table_ppt = bind_rows(final_list)
  
  # 设置表头，包含换行符 N 值
  colnames(final_table_ppt) <- c(
    "Characteristic", 
    paste0("All\n(N=", n_all, ")"), 
    paste0("No Shock\n(N=", n_no_shock, ")"), 
    paste0("Shock\n(N=", n_shock, ")")
  )
  
  # --- 4. 保存结果 ---
  output_file <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table\\Table3_icu_PPT格式统计表_P值合并版.xlsx"
  write.xlsx(final_table_ppt, output_file)
  
  print(paste("转换完成！文件保存在：", output_file))
}

{
  label_counts <- df_lab_log %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  
  #  将映射合并到你生成的 biomarker_table 中
  final_label_counts <- label_counts %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label)  
  
  setwd("D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table")
  write.xlsx(final_label_counts,"Table4_icu_中文分类filtered_mimic4_sepsis_实验室检查人数统计.xlsx")
}

# 2. 床边监护指标数据统计 #-----------------------------------------------------------------
first_icu_chartevent1 <- first_icu_chartevent %>%  mutate( 
  shock = is_shock_onset_in_this_icu, 
  hospital_expire_flag = case_when(this_icu_outcome == "Death in ICU" ~ 1,
                                   TRUE ~ 0 
  )
)
# 统计每个 label 的测量人数
{
  label_chartevent_counts <- first_icu_chartevent1 %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  
  output_dir <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table"
  file_name <- "Table5.ICU床边监护检查人数统计.xlsx"
  write.xlsx(label_counts, file.path(output_dir, file_name))
  
  setwd(Output_derived_data)
  saveRDS(first_icu_chartevent1,"ICU床边监护_filtered_mimic4_Sepsis.rds")
  
  chart_labels <- c(
    # --- 核心循环指标 (用于判定休克/血压) ---
    "Arterial Blood Pressure mean",
    "Arterial Blood Pressure systolic",
    "Arterial Blood Pressure diastolic",
    "PAEDP",                          # 肺动脉舒张末压 (反映左心前负荷)
    "Unassisted Systole",             # 反映自身心室功能
    
    # --- 神经系统指标 (用于 SOFA 评分) ---
    "GCS - Motor Response",
    "GCS - Verbal Response",
    
    # --- 呼吸系统指标 (用于氧合评估) ---
    "FiO2ApacheIIValue",              # 吸氧浓度
    "Total PEEP Level",               # PEEP水平
    "PSV Level",                      # 压力支持模式
    "PCV Level",                      # 压力控制模式
    
    # --- 肾脏支持与 CRRT (反映器官衰竭程度) ---
    "Filter Pressure",
    "Effluent Pressure",
    "Return Pressure",
    "Replacement Rate",
    "Dialysate Rate",
    
    # --- 实验室指标补充 (存在于 ChartEvents 中的) ---
    "Direct Bilirubin",
    "Chloride (serum)",
    "Total Protein",
    "CreatinineApacheIIValue",
    "CreatinineApacheIIScore"
  )
  
  first_icu_chartevent2 <- first_icu_chartevent1 %>%
    # 筛选指定的标签
    filter(label %in% chart_labels) %>%
    # 剔除 valuenum 为空或明显错误的数据 (例如血压为0或负数)
    filter(!is.na(valuenum)) %>%
    filter(!(label %in% c("Arterial Blood Pressure mean", "Arterial Blood Pressure systolic") & valuenum <= 0))  
  
  cat("筛选后的指标种类:", n_distinct( first_icu_chartevent2$label), "\n")
  cat("筛选后的总记录条数:", nrow( first_icu_chartevent2), "\n")
  
  # 查看各指标的数据密度，确认哪些指标最完整
  label_density <- first_icu_chartevent2 %>%
    group_by(label) %>%
    summarise(n_total = n(), n_patients = n_distinct(stay_id)) %>%
    arrange(desc(n_patients)) 
  
  print(label_density)
}

# 3. 床边监护指标统计表 #----------------------------------------------------
length(unique(first_icu_chartevent2$label))
unique(first_icu_chartevent2$label)
length(unique(first_icu_chartevent2$hadm_id))
length(unique(first_icu_chartevent2$hadm_id[first_icu_chartevent2$shock == 1]))
first_icu_chartevent3 <- first_icu_chartevent2 %>% filter(label %in% c("GCS - Motor Response", "GCS - Verbal Response","Arterial Blood Pressure mean"))
{
  library(dplyr)
  library(tidyr)
  library(broom) # 用于整洁地提取t检验结果
  
  # 定义生成包含 P 值的统计表函数
  generate_full_stat_table <- function(data) {
    
    # 1. 计算各组的基础统计量 (ALL, Non-Shock, Shock)
    # 计算分组统计 (0 和 1)
    group_stats <- data %>%
      group_by(label, shock) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = ifelse(shock == 1, "Shock", "Non-Shock"))
    
    # 计算全部统计 (ALL)
    all_stats <- data %>%
      group_by(label) %>%
      summarise(
        mean_val = mean(valuenum, na.rm = TRUE),
        sd_val   = sd(valuenum, na.rm = TRUE),
        med_val  = median(valuenum, na.rm = TRUE),
        min_val  = min(valuenum, na.rm = TRUE),
        max_val  = max(valuenum, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(group_name = "ALL")
    
    # 2. 核心步骤：对每个指标进行 t 检验计算 P 值
    p_values <- data %>%
      group_by(label) %>%
      filter(n_distinct(shock) == 2) %>% # 确保两组都有数据才能做检验
      do(tidy(t.test(valuenum ~ shock, data = .))) %>%
      select(label, p.value) %>%
      mutate(p_formatted = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ as.character(round(p.value, 3))
      ))
    
    # 3. 格式化统计字符串并合并 P 值
    final_stats <- bind_rows(group_stats, all_stats) %>%
      mutate(
        `Mean (SD)` = sprintf("%.2f (%.2f)", mean_val, sd_val),
        `Median [Min, Max]` = sprintf("%.2f [%.2f, %.2f]", med_val, min_val, max_val)
      ) %>%
      left_join(p_values, by = "label") %>%
      select(label, group_name, `Mean (SD)`, `Median [Min, Max]`, p_formatted)
    
    # 4. 转置表格为宽表格式
    table_wide <- final_stats %>%
      pivot_longer(cols = c(`Mean (SD)`, `Median [Min, Max]`), 
                   names_to = "Statistic", 
                   values_to = "Value") %>%
      pivot_wider(names_from = group_name, values_from = Value) %>%
      # 整理列顺序，并只在每个指标的第一行显示 P 值（美化）
      select(label, Statistic, ALL, `Non-Shock`, `Shock`, `P-value` = p_formatted) %>%
      mutate(`P-value` = ifelse(Statistic == "Median [Min, Max]", "", `P-value`)) # 第二行留白
    
    return(table_wide)
  }
  # --- 执行统计 ---
  biomarker_table <- generate_full_stat_table(first_icu_chartevent3)
  
  # 1. 创建分类与中英文对照映射表
  mapping_df <- data.frame(
    label = c(
      # --- 1. 心肌损伤与心脏功能指标 ---
      "GCS - Motor Response", "GCS - Verbal Response","Arterial Blood Pressure mean" 
    ),
    label_cn = c( "GCS - Motor Response(格拉斯哥昏迷评分-运动反应)", "GCS - Verbal Response(格拉斯哥昏迷评分-语言反应)", 
                  "Arterial Blood Pressure mean(有创动脉血压平均压)" 
    ),
    category = c(
      rep("床边监护", 3) 
    ),
    stringsAsFactors = FALSE
  )
  
  #  将映射合并到你生成的 biomarker_table 中
  final_ppt_table <- biomarker_table %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label) %>%
    select(category, label_cn, Statistic, ALL, `Non-Shock`, `Shock`, `P-value`)
  
  #  打印检查
  head(final_ppt_table)
}

{
  library(dplyr)
  library(tidyr)
  library(broom)
  
  # 1. 计算各组的总人数 N (假设 subject_id 唯一，或直接按数据记录数计算)
  # 这里我们统计每组的观察数，用于表头显示
  n_all <- nrow(first_icu_chartevent3 %>% distinct(hadm_id)) # 如果是按人统计
  n_no_shock <- nrow(first_icu_chartevent3 %>% filter(shock == 0) %>% distinct(hadm_id))
  n_shock <- nrow(first_icu_chartevent3 %>% filter(shock == 1) %>% distinct(hadm_id))
  
  # 如果数据是按记录算的，直接用 n_all <- nrow(df_lab_log) 等
  
  # 2. 定义格式化函数
  generate_ppt_table <- function(data) {
    
    # 获取所有指标列表
    labels <- unique(data$label)
    final_rows <- list()
    
    for (lab in labels) {
      # 提取单个指标数据
      sub_data <- data %>% filter(label == lab)
      
      # 计算 T 检验 P 值
      p_val <- t.test(valuenum ~ shock, data = sub_data)$p.value
      p_str <- if(p_val < 0.001) "P < 0.001" else sprintf("P = %.3f", p_val)
      
      # 计算三组统计量
      stats <- sub_data %>%
        group_by(shock) %>%
        summarise(
          mean_sd = sprintf("%.2f ( %.2f)", mean(valuenum, na.rm=T), sd(valuenum, na.rm=T)),
          med_range = sprintf("%.2f [ %.2f, %.2f]", median(valuenum, na.rm=T), min(valuenum, na.rm=T), max(valuenum, na.rm=T))
        )
      
      all_mean_sd <- sprintf("%.2f ( %.2f)", mean(sub_data$valuenum, na.rm=T), sd(sub_data$valuenum, na.rm=T))
      all_med_range <- sprintf("%.2f [ %.2f, %.2f]", median(sub_data$valuenum, na.rm=T), min(sub_data$valuenum, na.rm=T), max(sub_data$valuenum, na.rm=T))
      
      # 构建该指标的 4 行结构
      # 第一行：指标名称 (Label)
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      # 第四行：P 值
      
      df_block <- data.frame(
        Item = c(lab, "Mean (SD)", "Median [Min, Max]", p_str),
        All = c("", all_mean_sd, all_med_range, ""),
        No_shock = c("", stats$mean_sd[stats$shock == 0], stats$med_range[stats$shock == 0], ""),
        shock = c("", stats$mean_sd[stats$shock == 1], stats$med_range[stats$shock == 1], ""),
        stringsAsFactors = FALSE
      )
      
      final_rows[[lab]] <- df_block
    }
    
    # 合并所有指标块
    res_table <- bind_rows(final_rows)
    
    # 修改表头名称，加入 N 数值
    colnames(res_table) <- c(
      "", 
      paste0("All\n(N=", n_all, ")"), 
      paste0("No Cardiogenic Shock\n(N=", n_no_shock, ")"), 
      paste0("Cardiogenic Shock\n(N=", n_shock, ")")
    )
    
    return(res_table)
  }
  
  # --- 执行生成 ---
  ppt_stat_table <- generate_ppt_table(first_icu_chartevent3)
  
  # 查看结果
  print(ppt_stat_table)
}

{
  library(dplyr)
  library(openxlsx)
  
  # --- 1. 获取分组人数 N ---
  n_all <- 548
  n_shock <- 299
  n_no_shock <- n_all - n_shock
  
  # --- 2. 处理函数：将扁平表转换为块状 PPT 格式 ---
  # 请确保 final_ppt_table 包含列：category, label_cn, Statistic, ALL, Non-Shock, Shock, P-value
  df_input <- final_ppt_table 
  
  categories <- unique(df_input$category)
  final_list <- list()
  
  for (cat in categories) {
    
    # A. 插入分类大标题行 (例如：人口统计学特征)
    cat_header <- data.frame(
      Item = cat, 
      All = "", 
      No_Shock = "", 
      Shock = "", 
      stringsAsFactors = FALSE
    )
    final_list[[paste0(cat, "_header")]] <- cat_header
    
    # 获取该分类下的所有指标
    cat_data <- df_input %>% filter(category == cat)
    labels <- unique(cat_data$label_cn)
    
    for (lab in labels) {
      # 提取该指标的数据子集
      sub <- cat_data %>% filter(label_cn == lab)
      
      # 获取并格式化 P 值
      p_val_raw <- sub$`P-value`[1]
      p_suffix <- if(is.na(p_val_raw) || p_val_raw == "") "" else paste0(" (P ", p_val_raw, ")")
      
      # 组合 Label 和 P 值作为块的第一行
      display_label <- paste0(lab, p_suffix)
      
      # 构建 3 行块 (不再有单独的 P 值行)
      # 第一行：指标名称 + P值
      # 第二行：Mean (SD)
      # 第三行：Median [Min, Max]
      block <- data.frame(
        Item = c(display_label, "  Mean (SD)", "  Median [Min, Max]"), # 前面加空格方便PPT对齐
        All = c("", 
                sub$ALL[sub$Statistic == "Mean (SD)"], 
                sub$ALL[sub$Statistic == "Median [Min, Max]"]),
        No_Shock = c("", 
                     sub$`Non-Shock`[sub$Statistic == "Mean (SD)"], 
                     sub$`Non-Shock`[sub$Statistic == "Median [Min, Max]"]),
        Shock = c("", 
                  sub$`Shock`[sub$Statistic == "Mean (SD)"], 
                  sub$`Shock`[sub$Statistic == "Median [Min, Max]"]),
        stringsAsFactors = FALSE
      )
      
      final_list[[paste0(cat, "_", lab)]] <- block
    }
  }
  
  # --- 3. 合并所有行并设置表头 ---
  final_table_ppt = bind_rows(final_list)
  
  # 设置表头，包含换行符 N 值
  colnames(final_table_ppt) <- c(
    "Characteristic", 
    paste0("All\n(N=", n_all, ")"), 
    paste0("No Shock\n(N=", n_no_shock, ")"), 
    paste0("Shock\n(N=", n_shock, ")")
  )
  
  # --- 4. 保存结果 ---
  output_file <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table\\Table6_icu_chartevents_PPT格式统计表_P值合并版.xlsx"
  write.xlsx(final_table_ppt, output_file)
  
  print(paste("转换完成！文件保存在：", output_file))
}

{
  label_counts <- first_icu_chartevent3 %>%
    group_by(label) %>%
    summarise(
      n_measurements = n(),  # 总测量次数
      n_patients = n_distinct(hadm_id),  # 不同患者数量（添加逗号）
      measurements_per_patient = n_measurements / n_patients  # 平均每人测量次数
    ) %>%
    arrange(desc(measurements_per_patient ))  # 按患者数降序排列 
  
  #  将映射合并到你生成的 biomarker_table 中
  final_label_counts <- label_counts %>%
    left_join(mapping_df, by = "label") %>%
    # 使用中文名替换原始 label，并根据类别排序
    arrange(category, label)  
  
  setwd("D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\EDA_table")
  write.xlsx(final_label_counts,"Table7_icu_床边监护_中文分类filtered_mimic4_sepsis_实验室检查人数统计.xlsx")
}

# 休克组和非休克组实验室指标对比 #----------------------------------------------
setwd(file.path(Output_derived_data, "data_filtered"))
df_lab_log <- readRDS("icu_filtered_mimic4_Sepsislog.rds") 

unique(df_lab_log$label)
target_labelx1 <- c(
  # --- 1. 心肌损伤与心脏功能指标 ---
  "Troponin T", "NTproBNP", "Creatine Kinase (CK)", "Creatine Kinase, MB Isoenzyme", "CK-MB Index", "proBNP, Pleural"
)

target_labelx2 <- c(
  # --- 2. 组织灌注与酸碱平衡指标 ---
  "Lactate", "Lactate Dehydrogenase (LD)", "Anion Gap", "Base Excess", "Bicarbonate", "Calculated Total CO2", 
  "pCO2", "pO2", "pH", "Oxygen Saturation", "Alveolar-arterial Gradient"
)
target_labelx3 <- c( 
  # --- 3. 肝肾器官功能指标 ---
  "Alanine Aminotransferase (ALT)", "Asparate Aminotransferase (AST)", "Bilirubin, Direct", "Bilirubin, Indirect", 
  "Bilirubin, Total", "Albumin", "Creatinine", "Creatinine, Urine", "Urea Nitrogen", "Uric Acid" )
target_labelx4 <- c(
  # --- 4. 血液学与凝血功能指标 ---
  "White Blood Cells", "Red Blood Cells", "Hemoglobin", "Hematocrit", "Platelet Count", "RDW", "RDW-SD", 
  "Neutrophils", "Lymphocytes", "Monocytes", "Immature Granulocytes", "Metamyelocytes", "Myelocytes", 
  "Bands", "Atypical Lymphocytes", "Nucleated Red Cells", "PT", "PTT", "INR(PT)", "Fibrinogen, Functional", "Antithrombin"
)
target_labelx5<- c( 
  # --- 5. 电解质与矿物质指标 ---
  "Potassium", "Potassium, Whole Blood", "Sodium", "Sodium, Whole Blood", "Chloride", "Chloride, Whole Blood", 
  "Calcium, Total", "Free Calcium", "Magnesium", "Phosphate")
target_labelx6 <- c( 
  # --- 6. 炎症、代谢及其他指标 ---
  "C-Reactive Protein", "Glucose", "Osmolality, Measured", "Required O2", "Total Protein, Urine"
)

# 1. 配置参数与路径
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)

base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots"

# 1. 定义一个统一的 P 值格式化函数
format_p_value <- function(p) {
  if (is.na(p)) return("P = N/A")
  if (p < 0.001) return("P < 0.001")
  return(paste0("P = ", round(p, 3)))
}

# 2. 开始大循环
for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    
    message("正在绘制并计算双重显著性: ", target_biomarker)
    colnames(df_lab_log)
    # --- A. 数据预处理 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        Group = factor(ifelse(shock == 1, "Shock Group", "Non-Shock Group"), 
                       levels = c("Non-Shock Group", "Shock Group")),
        # 为 Log 统计准备数据 (加极小值防止 log(0) 报错)
        val_log = log10(valuenum + 0.001),
        Phase = case_when(
          shock == 0 ~ "Non-Shock",
          charttime_rel_icu < shock_onset_rel_icu ~ "Pre-Shock",
          charttime_rel_icu >= shock_onset_rel_icu ~ "Post-Shock"
        )
      )
    
    if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) next
    
    #计算显著性 ---
    # 1. 原始数值的显著性 (Linear Scale P)
    p_raw_text <- "P = N/A"
    try({
      fit_raw <- lmer(valuenum ~ Group + charttime_rel_icu + (1|subject_id), data = plot_df)
      p_raw <- summary(fit_raw)$coefficients["GroupShock Group", "Pr(>|t|)"]
      p_raw_text <- format_p_value(p_raw)
    }, silent = TRUE)
    
    # 2. 对数转换后的显著性 (Log-scale P)
    p_log_text <- "P = N/A"
    try({
      fit_log <- lmer(val_log ~ Group + charttime_rel_icu + (1|subject_id), data = plot_df)
      p_log <- summary(fit_log)$coefficients["GroupShock Group", "Pr(>|t|)"]
      p_log_text <- format_p_value(p_log)
    }, silent = TRUE)
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # --- C. 绘图：线性坐标图 (标注原始 P) ---
    p1_linear <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(
        subtitle = paste0("Linear-scale diff: ", p_raw_text),
        x = "Days from Admission", y = paste0(target_biomarker, " (", unit, ")")
      ) +
      theme_cor + ggtitle("Linear Scale")
    
    # --- D. 绘图：对数坐标图 (标注对数 P) ---
    p1_log <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_y_log10() + 
      scale_color_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_fill_manual(values = c("Shock Group" = "#E41A1C", "Non-Shock Group" = "#377EB8")) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(
        subtitle = paste0("Log-scale diff: ", p_log_text), # 这里显示 Log P
        x = "Days from Admission", y = paste0("Log10 [ ", target_biomarker, " ]")
      ) +
      theme_cor + ggtitle("Log10 Scale")
    
    # --- E. 合并与保存 ---
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_Dual_P.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}
message("所有带显著性检验的轨迹图已保存。")

# 死亡组和非死亡组实验室指标对比 #----------------------------------------------

# 1. 辅助函数：格式化 P 值显示
format_p <- function(p) {
  if (is.na(p)) return("P = N/A")
  if (p < 0.001) return("P < 0.001")
  return(paste0("P = ", round(p, 3)))
}

# 2. 配置参数与新路径
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_Death_vs_Survival"

# 3. 开始大循环
for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    message("正在绘制死亡 vs 存活显著性: ", target_biomarker)
    
    # --- A. 数据预处理 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        # 【修改关键点：使用 hospital_expire_flag 分组】
        # 0 为存活 (Survivors), 1 为死亡 (Non-Survivors)
        Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                       levels = c("Survivors", "Non-Survivors")),
        
        # 这里的 Phase 逻辑保留，反映数据点在休克发生的前后关系（如果有的话）
        Phase = case_when(
          shock == 0 ~ "Non-Shock",
          charttime_rel_icu < shock_onset_rel_icu ~ "Pre-Shock",
          charttime_rel_icu >= shock_onset_rel_icu ~ "Post-Shock"
        )
      )
    
    # 健壮性检查
    if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
      message("跳过 ", target_biomarker, "：有效对比数据不足")
      next
    }
    
    # --- B. 计算显著性 (LMM) ---
    p_linear <- NA
    p_log_val <- NA
    
    # 线性模型显著性：Group 对数值的影响
    try({
      model_linear <- lmer(valuenum ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
      # 【提取关键：行名变为 GroupNon-Survivors】
      p_linear <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 对数模型显著性
    try({
      model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
      p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # --- C. 绘图部分 ---
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # 定义颜色：存活组蓝色，死亡组红色
    outcome_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
    
    # 线性坐标图
    p1_linear <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = outcome_colors) +
      scale_fill_manual(values = outcome_colors) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      labs(x = "Days from Admission", y = paste0(target_biomarker, " (", unit, ")"),
           subtitle = paste0("Outcome Difference: ", format_p(p_linear))) +
      theme_cor + ggtitle("Linear Scale")
    
    # 对数坐标图
    p1_log <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(aes(shape = Phase), alpha = 0.2, size = 1) + 
      geom_smooth(method = "loess", size = 0.8, alpha = 0.2) + 
      scale_color_manual(values = outcome_colors) +
      scale_fill_manual(values = outcome_colors) +
      scale_shape_manual(values = c("Non-Shock" = 16, "Pre-Shock" = 1, "Post-Shock" = 17)) +
      scale_y_log10() + 
      labs(x = "Days from Admission", y = paste0("Log10 [ ", target_biomarker, " ]"),
           subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
      ) +
      theme_cor + ggtitle("Log10 Scale")
    
    # 合并与保存
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_Death_vs_Survivors.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有死亡 vs 存活对比图已完成保存。") 
# 对齐休克为0时刻 休克组实验室指标ct图#---------------------------------

# 定义组名和对应的 label 列表
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)

base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_对齐休克0时刻"

for (group_name in names(label_groups)) {
  # group_name = "心肌损伤与心脏功能指标"
  # 创建子文件夹
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    # target_biomarker = "Troponin T"
    message("正在绘制对齐图 [", group_name, "]: ", target_biomarker)
    
    plot_df <- df_lab_log  %>% 
      filter(shock == 1)  %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      mutate(
        
        # 核心：时间对齐逻辑 (T=0 为休克发生点或入ICU点)
        time_aligned = ifelse(shock == 1, 
                              charttime_rel_icu - shock_onset_rel_icu, 
                              charttime_rel_icu - icu_intime_first_rel),
        # 定义阶段
        Phase = case_when(
          time_aligned < 0 ~ "Pre-Shock",
          time_aligned >= 0 ~ "Post-Shock"
        ),
        # 用于分段画趋势线
        LineGroup = interaction(Phase)
      ) 
    
    # 健壮性检查
    if (nrow(plot_df) < 10) {
      message("跳过 ", target_biomarker, "：有效数据点太少")
      next
    }
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color ="#E41A1C" )) +
      geom_point(aes(shape = Phase), alpha = 0.3, size = 1, stroke = 0) + 
      # geom_line(aes(group = hadm_id),size = 0.8, alpha = 0.2) +
      # 关键：group = LineGroup 实现 0 点处的断开拟合
      geom_smooth(aes(group = LineGroup), method = "loess", size = 0.8, alpha = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) + 
      scale_shape_manual(values = c("Pre-Shock" = 16, "Post-Shock" = 16)) +
      labs(x = "Days relative to Shock Onset (T=0)", y = paste0(target_biomarker, " (", unit, ")")) +
      theme_cor + ggtitle("Linear Scale Aligned")
    
    p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = "#E41A1C" )) +
      geom_point(aes(shape = Phase), alpha = 0.3, size = 1, stroke = 0) + 
      # geom_line(aes(group = hadm_id),size = 0.8, alpha = 0.2) +
      geom_smooth(aes(group = LineGroup), method = "loess", size = 0.8, alpha = 0.2) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      scale_y_log10() +  
      scale_shape_manual(values = c("Pre-Shock" = 16, "Post-Shock" = 16)) +
      labs(x = "Days relative to Shock Onset (T=0)", y = paste0("Log10 [ ", target_biomarker, " ]")) +
      theme_cor + ggtitle("Log10 Scale Aligned")
    
    combined_plot <- (p1_linear | p1_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "none")
    
    # 文件名清洗
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    ggsave(file.path(output_dir, paste0(safe_name, "_Aligned_T0.png")), 
           plot = combined_plot, width = 12, height = 6, dpi = 300)
    
    
  }
}

message("所有 6 组对齐图片已处理完毕并按文件夹分类保存。")
colnames(plot_df)
# 对齐休克为0时刻 心源性休克组 vs 非心源性休克实验室指标ct图#---------------------------------
{
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_对齐休克为0时刻_CS_vs_SS"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>%
        mutate(
          # 重新定义对比组，并设置为 Factor（固定基准组为 Non-Cardio）
          Group = ifelse(Cardiogenic_shock == 1, "Cardiogenic Shock", "Non-Cardio Shock"),
          Group = factor(Group, levels = c("Non-Cardio Shock", "Cardiogenic Shock")),
          
          # 统一以休克时刻为 T=0
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        )  
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足")
        next
      }
      
      # --- B. 计算显著性 (LMM) ---
      p_linear <- NA
      p_log_val <- NA
      
      # 线性模型显著性：Group 对数值的影响
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear <- summary(model_linear)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置
      group_colors <- c("Cardiogenic Shock" = "#377EB8", "Non-Cardio Shock" = "#E41A1C")
      
      # --- 4. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Outcome Difference: ", format_p(p_log_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Shock Type", fill = "Shock Type"
        ) +
        theme_cor +
        theme(legend.position = "none") 
      
      # --- 5. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]"),
          color = "Shock Type", fill = "Shock Type",
          subtitle = paste0("Outcome Difference: ", format_p(p_linear))
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 6. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_vs_NonCS.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("所有心源性 vs 非心源性对比图已保存。")
} 
# 对齐休克为0时刻：死亡组 vs 存活组 实验室指标ct图对比 #---------------------------------
{
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 修改输出根目录为 Outcome 比较
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_对齐休克为0时刻_Death_vs_Survival"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制结局对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        # 仅针对休克患者观察其结局差异
        filter(shock == 1) %>%
        mutate(
          
          Group = ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"),
          Group = factor(Group, levels = c("Survivors", "Non-Survivors")),
          
          # 统一以休克时刻为 T=0
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          
          # 定义阶段
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          # 用于分段画趋势线
          LineGroup = interaction(Group, Phase)
        )  
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足")
        next
      }
      
      p_linear <- NA
      p_log_val <- NA
      
      # 线性模型显著性：Group 对数值的影响
      try({
        model_linear <- lmer(valuenum ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
        # 【提取关键：行名变为 GroupNon-Survivors】
        p_linear <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # --- 颜色设置：死亡红色，存活蓝色 ---
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- 1. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Outcome Difference: ", format_p(p_linear)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Hospital Outcome", fill = "Hospital Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 2. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]"),
          color = "Hospital Outcome", fill = "Hospital Outcome",
          subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 3. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_Outcome_Comparison.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("所有死亡 vs 存活对比图已完成。")
}
# 对齐休克为0时刻：非心源性休克：死亡组 vs 存活组 实验室指标ct图对比 #---------------------------------
{
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  
  # 辅助函数：格式化P值
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) return("N/A")
    if (p < 0.001) return("< 0.001")
    return(as.character(round(p, 3)))
  }
  
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 修改输出根目录
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_对齐休克为0时刻心源性休克_Death_vs_Survival"
  
  for (group_name in names(label_groups)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制心源性休克结局对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 1. 数据预处理 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1 & Cardiogenic_shock == 1) %>% # 严格限定为心源性休克人群
        mutate(
          # 定义结局组，Survivors作为基准
          Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                         levels = c("Survivors", "Non-Survivors")),
          # 时间对齐 T=0
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        ) 
      # filter(time_aligned >= -7 & time_aligned <= 7) # 观察核心窗口
      
      # 健壮性检查
      if (nrow(plot_df) < 10 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：有效数据不足或单一结局")
        next
      }
      
      # --- 2. 计算显著性 (LMM) ---
      p_linear_val <- NA
      p_log_val <- NA
      
      # 线性模型显著性
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear_val <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置：死亡红色，存活蓝色
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- 3. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Linear Scale",
          subtitle = paste0("Overall Diff P: ", format_p(p_linear_val)), # 修正这里
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Hospital Outcome", fill = "Hospital Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 4. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.15, size = 0.8) + 
        geom_smooth(aes(group = LineGroup), method = "loess", size = 1, alpha = 0.2) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = "Log10 Scale",
          subtitle = paste0("Overall Diff (Log) P: ", format_p(p_log_val)), # 修正这里
          x = "Days relative to Shock Onset (T=0)",
          y = paste0("Log10 [ ", target_biomarker, " ]")
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 5. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      # 文件名处理
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_Outcome_T0.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("心源性休克死亡 vs 存活所有轨迹图已绘制完毕。")
}
# 休克发生前两组患者的实验室指标对比 #------------------------------------------
label_groups <- list(
  "心肌损伤与心脏功能指标" = target_labelx1,
  "组织灌注与酸碱平衡指标" = target_labelx2,
  "肝肾器官功能指标"       = target_labelx3,
  "血液学与凝血功能指标"   = target_labelx4,
  "电解质与矿物质指标"     = target_labelx5,
  "炎症、代谢及其他指标"   = target_labelx6
)
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_PreShock_13days"
format_p <- function(p) {
  if (is.na(p)) return("P = N/A")
  if (p < 0.001) return("P < 0.001")
  return(paste0("P = ", round(p, 3)))
}
for (group_name in names(label_groups)) {
  
  output_dir <- file.path(base_path, group_name)
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  target_list <- label_groups[[group_name]]
  
  for (target_biomarker in target_list) {
    
    message("正在绘制休克前对比图: ", target_biomarker)
    
    # --- 核心逻辑：剔除休克后的测量点 ---
    plot_df <- df_lab_log %>%
      filter(label == target_biomarker & !is.na(valuenum)) %>%
      filter(
        # 条件1：如果是休克组，只取入院后、发病前的点
        (shock == 1 & charttime_rel_icu >= 0 & charttime_rel_icu < shock_onset_rel_icu) |
          # 条件2：如果是非休克组，只取入院后的点
          (shock == 0 & charttime_rel_icu >= 0)
      ) %>% 
      # # 限制观察时间范围，比如入院后前 7 天（防止长住院病人的尾部噪音）
      filter(charttime_rel_icu <= 13) %>%
      mutate(Group = factor(ifelse(shock == 1, "Future Shock", "Non-Shock"),
                            levels = c("Non-Shock", "Future Shock")))
    
    # 健壮性检查
    if (nrow(plot_df) < 10 | length(unique(plot_df$shock)) < 2) {
      message("跳过 ", target_biomarker, "：有效对比数据太少")
      next
    }
    
    # 线性模型显著性：Group 对数值的影响
    try({
      model_linear <- lmer(valuenum ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
      
      p_linear <- summary(model_linear)$coefficients["GroupFuture Shock", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 对数模型显著性
    try({
      model_log <- lmer(log10(valuenum + 0.001) ~ Group + charttime_rel_icu + (1 | subject_id), data = plot_df)
      p_log_val <- summary(model_log)$coefficients["GroupFuture Shock", "Pr(>|t|)"]
    }, silent = TRUE)
    
    # 获取单位
    unit <- plot_df$valueuom[1]
    if (is.na(unit) || unit == "") unit <- "Value"
    
    # --- 2. 线性坐标图 (仅含休克前点) ---
    p_linear <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(
        x = "Days from Admission",
        y = paste0(target_biomarker, " (", unit, ")"),
        subtitle = paste0("Outcome Difference: ", format_p(p_linear))
      ) +
      theme_cor + ggtitle("Linear Scale (Pre-Shock Only)")
    
    # --- 3. 对数坐标图 (仅含休克前点) ---
    p_log <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
      geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
      geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
      scale_y_log10() +
      scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
      labs(x = "Days from Admission", y = paste0("Log10[ ",target_biomarker, " (", unit, ") ]"),
           subtitle = paste0("Outcome Difference: ", format_p(p_log_val))
      ) +
      theme_cor + ggtitle("Log10 Scale (Pre-Shock Only)")
    
    # --- 4. 合并与保存 ---
    combined_plot <- (p_linear | p_log) + 
      plot_layout(guides = 'collect') & 
      theme(legend.position = "bottom")
    
    safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
    file_path <- file.path(output_dir, paste0(safe_name, "_PreShock_Comparison.png"))
    ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
  }
}

message("所有休克前对比图绘制完毕。")

# 休克发生后两组患者的实验室指标对比 #------------------------------------------
{
  # 辅助函数：格式化P值
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) return("N/A")
    if (p < 0.001) return("< 0.001")
    return(as.character(round(p, 3)))
  }
  
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 建议存放在专门的 "Pre-Shock_Comparison" 文件夹
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_post_Shock"
  
  for (group_name in names(label_groups)) {
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      message("正在处理: ", target_biomarker)
      
      # 1. 数据准备
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(
          (shock == 1 & charttime_rel_icu >= 0 & charttime_rel_icu < shock_onset_rel_icu) | # 注意：Pre-Shock应该是小于onset
            (shock == 0 & charttime_rel_icu >= 0)
        ) %>% 
        mutate(Group = factor(ifelse(shock == 1, "Future Shock", "Non-Shock"),
                              levels = c("Non-Shock", "Future Shock")))
      
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据量不足")
        next
      }
      
      # --- 2. 统计建模 (移入循环内部) ---
      p_linear_val <- "N/A"
      p_log_val <- "N/A"
      
      # 线性模型：加入交互项 Group * charttime_rel_icu
      try({
        model_linear <- lmer(valuenum ~ Group * charttime_rel_icu + (1 | subject_id), data = plot_df)
        coef_summary <- summary(model_linear)$coefficients
        # 提取 Group 的主效应 P 值 (整体差异)
        p_linear_val <- format_p(coef_summary["GroupFuture Shock", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # 对数模型：处理偏态数据（如肌钙蛋白、乳酸）
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group * charttime_rel_icu + (1 | subject_id), data = plot_df)
        coef_log_summary <- summary(model_log)$coefficients
        p_log_val <- format_p(coef_log_summary["GroupFuture Shock", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # --- 3. 绘图 ---
      # 线性图
      p_linear_plot <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.1, size = 0.5) + 
        geom_smooth(method = "loess", size = 1.2) + 
        scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
        scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
        labs(x = "Days", y = target_biomarker, 
             title = paste0("Linear (P: ", p_linear_val, ")")) +
        theme_cor
      
      # 对数图
      p_log_plot <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.1, size = 0.5) + 
        geom_smooth(method = "loess", size = 1.2) + 
        scale_y_log10() +
        scale_color_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
        scale_fill_manual(values = c("Future Shock" = "#E41A1C", "Non-Shock" = "#377EB8")) +
        labs(x = "Days", y = paste0("Log10 ", target_biomarker),
             title = paste0("Log-scale (P: ", p_log_val, ")")) +
        theme_cor
      
      # 合并保存
      combined_plot <- (p_linear_plot | p_log_plot) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      ggsave(file.path(output_dir, paste0(safe_name, "_PreShock.png")), combined_plot, width = 12, height = 6)
    }
  }
}
# 休克发生后休克患者死亡存活组的实验室指标对比 #------------------------------------------
# 休克发生后：对比休克患者中的 死亡组 (4) vs 存活组 (3)
{
  # 1. 辅助函数：格式化P值
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) return("N/A")
    if (p < 0.001) return("< 0.001")
    return(as.character(round(p, 3)))
  }
  
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 设置输出路径
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_post_Shock_survival_death"
  
  for (group_name in names(label_groups)) {
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      message("正在处理: ", target_biomarker)
      
      # --- 2. 数据准备 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>%                      # 仅保留休克患者
        filter(charttime_rel_icu >= shock_onset_rel_icu) %>% # 仅保留休克发生后的数据点
        # 时间对齐：以休克发生为 0 点
        mutate(Days_Post_Shock = charttime_rel_icu - shock_onset_rel_icu) %>%
        # 分组逻辑修改：4为死亡，3为存活
        mutate(Group = factor(ifelse(endpoint_state == 4, "Deceased", "Survivor"),
                              levels = c("Survivor", "Deceased")))
      
      # 健壮性检查
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足以对比存活与死亡")
        next
      }
      
      # --- 3. 统计建模 ---
      p_linear_val <- "N/A"
      p_log_val <- "N/A"
      
      # 线性混合模型：检验两组之间是否存在显著性差异
      try({
        model_linear <- lmer(valuenum ~ Group * Days_Post_Shock + (1 | subject_id), data = plot_df)
        coef_summary <- summary(model_linear)$coefficients
        # 提取 GroupDeceased 的 P 值
        p_linear_val <- format_p(coef_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # 对数混合模型
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group * Days_Post_Shock + (1 | subject_id), data = plot_df)
        coef_log_summary <- summary(model_log)$coefficients
        p_log_val <- format_p(coef_log_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # --- 4. 绘图 ---
      # 存活组 Survivor(3) 用绿色，死亡组 Deceased(4) 用红色
      color_values <- c("Deceased" = "#E41A1C", "Survivor" = "#377EB8")
      
      p_linear_plot <- ggplot(plot_df, aes(x = Days_Post_Shock, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Days post shock onset", y = target_biomarker, 
             title = paste0("Linear (Group P: ", p_linear_val, ")")) +
        theme_cor
      
      p_log_plot <- ggplot(plot_df, aes(x = Days_Post_Shock, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_y_log10() +
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Days post shock onset", y = paste0("Log10 ", target_biomarker),
             title = paste0("Log-scale (Group P: ", p_log_val, ")")) +
        theme_cor
      
      # 合并保存
      combined_plot <- (p_linear_plot | p_log_plot) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      ggsave(file.path(output_dir, paste0(safe_name, "_PostShock_Survivor_vs_Death.png")), 
             combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
}
# 休克发生后休克患者死亡存活组的实验室指标对比30days #------------------------------------------
# 休克发生后：对比休克患者中的 死亡组 (4) vs 存活组 (3)
{
  # 1. 辅助函数：格式化P值
  format_p <- function(p) {
    if (is.null(p) || is.na(p)) return("N/A")
    if (p < 0.001) return("< 0.001")
    return(as.character(round(p, 3)))
  }
  
  label_groups <- list(
    "心肌损伤与心脏功能指标" = target_labelx1,
    "组织灌注与酸碱平衡指标" = target_labelx2,
    "肝肾器官功能指标"       = target_labelx3,
    "血液学与凝血功能指标"   = target_labelx4,
    "电解质与矿物质指标"     = target_labelx5,
    "炎症、代谢及其他指标"   = target_labelx6
  )
  
  # 设置输出路径
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_post_Shock_survival_death_30days"
  
  for (group_name in names(label_groups)) {
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      message("正在处理: ", target_biomarker)
      
      # --- 2. 数据准备 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>%                      # 仅保留休克患者
        filter(charttime_rel_icu >= shock_onset_rel_icu) %>% # 仅保留休克发生后的数据点
        # 时间对齐：以休克发生为 0 点
        mutate(Days_Post_Shock = charttime_rel_icu - shock_onset_rel_icu) %>%
        # 分组逻辑修改：4为死亡，3为存活
        mutate(Group = factor(ifelse(endpoint_state == 4, "Deceased", "Survivor"),
                              levels = c("Survivor", "Deceased"))) %>% 
        filter(charttime_rel_icu <= 30)
      
      # 健壮性检查
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：数据不足以对比存活与死亡")
        next
      }
      
      # --- 3. 统计建模 ---
      p_linear_val <- "N/A"
      p_log_val <- "N/A"
      
      # 线性混合模型：检验两组之间是否存在显著性差异
      try({
        model_linear <- lmer(valuenum ~ Group * Days_Post_Shock + (1 | subject_id), data = plot_df)
        coef_summary <- summary(model_linear)$coefficients
        # 提取 GroupDeceased 的 P 值
        p_linear_val <- format_p(coef_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # 对数混合模型
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group * Days_Post_Shock + (1 | subject_id), data = plot_df)
        coef_log_summary <- summary(model_log)$coefficients
        p_log_val <- format_p(coef_log_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # --- 4. 绘图 ---
      # 存活组 Survivor(3) 用绿色，死亡组 Deceased(4) 用红色
      color_values <- c("Deceased" = "#E41A1C", "Survivor" = "#377EB8")
      
      p_linear_plot <- ggplot(plot_df, aes(x = Days_Post_Shock, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Days post shock onset", y = target_biomarker, 
             title = paste0("Linear (Group P: ", p_linear_val, ")")) +
        theme_cor
      
      p_log_plot <- ggplot(plot_df, aes(x = Days_Post_Shock, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_y_log10() +
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Days post shock onset", y = paste0("Log10 ", target_biomarker),
             title = paste0("Log-scale (Group P: ", p_log_val, ")")) +
        theme_cor
      
      # 合并保存
      combined_plot <- (p_linear_plot | p_log_plot) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      ggsave(file.path(output_dir, paste0(safe_name, "_PostShock_Survivor_vs_Death.png")), 
             combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
}
# 4. 床边监护ct图 #-------------------------------------------------------------
# 对齐休克为0时刻 心源性休克组 vs 非心源性休克床边监护指标ct图#---------------------------------
# 1. 提取 hadm_id 与 Cardiogenic_shock 的对应关系
df_cs_mapping <- first_icu_lab %>% select(hadm_id, Cardiogenic_shock) %>% distinct() %>%
  mutate(hadm_id = as.numeric(hadm_id))

# 2. 将标签关联到监护数据表中
first_icu_chartevent3 <- first_icu_chartevent3 %>% mutate(hadm_id = as.numeric(hadm_id)) %>%
  left_join(df_cs_mapping, by = "hadm_id") %>% 
  filter(!(label == "Arterial Blood Pressure mean" & valuenum > 300))

{
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(mgcv) # gam平滑需要
  
  # 辅助函数：格式化 P 值
  format_p <- function(p) {
    if (is.na(p)) return("P = N/A")
    if (p < 0.001) return("P < 0.001")
    return(paste0("P = ", round(p, 3)))
  }
  
  hemo_labels <- list( "循环与意识监护指标" = c("GCS - Verbal Response", "GCS - Motor Response", "Arterial Blood Pressure mean"))
  
  # 修改输出路径
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Hemo_Plots_Aligned_CS_vs_SS"
  
  for (group_name in names(hemo_labels)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- hemo_labels[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制对比图 [", group_name, "]: ", target_biomarker)
      
      # --- 1. 数据预处理 ---
      plot_df <- first_icu_chartevent3 %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>%
        mutate(
          # 定义组别并锁定因子水平
          Group = factor(ifelse(Cardiogenic_shock == 1, "Cardiogenic Shock", "Non-Cardio Shock"),
                         levels = c("Non-Cardio Shock", "Cardiogenic Shock")),
          # 时间对齐 (T=0 为休克发生点)
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        ) %>%
        # 限制观察窗口：休克前后 7 天
        filter(time_aligned >= -7 & time_aligned <= 7)
      
      # 健壮性检查
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：有效数据不足")
        next
      }
      
      # --- 2. 计算显著性 (LMM) ---
      p_linear_text <- "P = N/A"
      p_log_text <- "P = N/A"
      
      # 线性模型显著性
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_lin <- summary(model_linear)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
        p_linear_text <- format_p(p_lin)
      }, silent = TRUE)
      
      # 对数模型显著性
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log <- summary(model_log)$coefficients["GroupCardiogenic Shock", "Pr(>|t|)"]
        p_log_text <- format_p(p_log)
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$unitname[1] # 注意监护表可能叫 unitname 或 valueuom
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置
      group_colors <- c("Cardiogenic Shock" = "#377EB8", "Non-Cardio Shock" = "#E41A1C")
      
      # --- 3. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + # 监护数据点多，降低alpha和size
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6)) +
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + # 使用 gam 拟合大数据
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Linear)"),
          subtitle = paste0("Group Diff: ", p_linear_text),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Shock Type", fill = "Shock Type"
        ) +
        theme_bw() + 
        theme(legend.position = "none") 
      
      # --- 4. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + 
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + 
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Log10)"),
          subtitle = paste0("Group Diff (Log): ", p_log_text),
          x = "Days relative to Shock Onset (T=0)",
          y = "Log10 Scale"
        ) +
        theme_bw() +
        theme(legend.position = "none")
      
      # --- 5. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_vs_NonCS.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("监护指标对齐图已全部生成。")
  }
# 对齐休克为0时刻：死亡组 vs 存活组 监护指标ct图对比 #---------------------------------
{
  
  # 修改输出根目录
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Hemo_Plots_Death_vs_Survival"
  
  for (group_name in names(hemo_labels)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- hemo_labels[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制结局对比图: ", target_biomarker)
      
      # --- A. 数据预处理 ---
      plot_df <- first_icu_chartevent3 %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1) %>% # 仅针对休克患者观察结局
        # 血压异常值清洗（防止 7559 等噪声拉歪曲线）
        filter(!(label == "Arterial Blood Pressure mean" & (valuenum <= 0 | valuenum > 300))) %>%
        mutate(
          # 定义结局组，Survivors为基准
          Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                         levels = c("Survivors", "Non-Survivors")),
          # 时间轴对齐 T=0
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        ) %>%
        # 观察窗口：休克前后 7 天
        filter(time_aligned >= -7 & time_aligned <= 7)
      
      # 健壮性检查
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：有效数据不足")
        next
      }
      
      # --- B. 计算整体显著性 (LMM) ---
      p_linear_val <- NA
      p_log_val <- NA
      
      # 线性空间 LMM
      try({
        fit_lin <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear_val <- summary(fit_lin)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数空间 LMM (处理 Sepsis 指标常见的偏态)
      try({
        fit_log <- lmer(log10(valuenum + 0.1) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(fit_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$unitname[1]
      if (is.na(unit) || unit == "") unit <- "Value"
      
      # 颜色设置：死亡红色，存活蓝色
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- C. 绘图：线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + # 降低 alpha 处理海量点
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + 
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Linear)"),
          subtitle = paste0("Outcome Diff P: ", format_p(p_linear_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Outcome", fill = "Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- D. 绘图：对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + 
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + 
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Log10)"),
          subtitle = paste0("Outcome Diff (Log) P: ", format_p(p_log_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = "Log10 Scale"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- E. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_Outcome_Trajectory.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("✅ 监护指标结局对比图绘制完成，保存在：", base_path)
} 
# 对齐休克为0时刻：心源性休克：死亡组 vs 存活组 监护指标ct图对比 #---------------------------------
{
  # 修改输出根目录
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Hemo_Plots_CS_Death_vs_Survival"
  
  for (group_name in names(hemo_labels)) {
    
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    # 使用监护指标列表
    target_list <- hemo_labels[[group_name]]
    
    for (target_biomarker in target_list) {
      
      message("正在绘制心源性休克结局对比图 [监护]: ", target_biomarker)
      
      # --- 1. 数据预处理 ---
      plot_df <- first_icu_chartevent3 %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 1 & Cardiogenic_shock == 1) %>% 
        # 针对血压指标进行初步清洗，防止异常值
        filter(!(grepl("Pressure", label) & (valuenum <= 0 | valuenum > 300))) %>%
        mutate(
          Group = factor(ifelse(hospital_expire_flag == 1, "Non-Survivors", "Survivors"), 
                         levels = c("Survivors", "Non-Survivors")),
          time_aligned = charttime_rel_icu - shock_onset_rel_icu,
          Phase = ifelse(time_aligned < 0, "Pre-Shock", "Post-Shock"),
          LineGroup = interaction(Group, Phase)
        ) %>%
        # 限制观察时间：休克前后 7 天
        filter(time_aligned >= -7 & time_aligned <= 7)
      
      # 健壮性检查
      if (nrow(plot_df) < 20 || length(unique(plot_df$Group)) < 2) {
        message("跳过 ", target_biomarker, "：有效数据不足")
        next
      }
      
      # --- 2. 计算显著性 (LMM) ---
      p_linear_val <- NA
      p_log_val <- NA
      
      # 线性模型显著性
      try({
        model_linear <- lmer(valuenum ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_linear_val <- summary(model_linear)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 对数模型显著性 (处理偏态，GCS等指标不建议用对数，但逻辑保留)
      try({
        model_log <- lmer(log10(valuenum + 0.1) ~ Group + time_aligned + (1 | subject_id), data = plot_df)
        p_log_val <- summary(model_log)$coefficients["GroupNon-Survivors", "Pr(>|t|)"]
      }, silent = TRUE)
      
      # 获取单位
      unit <- plot_df$valueuom[1]
      if (is.na(unit) || unit == "") unit <- "unit"
      
      # 颜色设置：死亡红色，存活蓝色
      group_colors <- c("Non-Survivors" = "#E41A1C", "Survivors" = "#377EB8")
      
      # --- 3. 线性坐标图 ---
      p1_linear <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + # 极低透明度
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + 
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Linear)"),
          subtitle = paste0("Overall Diff P: ", format_p(p_linear_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = paste0(target_biomarker, " (", unit, ")"),
          color = "Outcome", fill = "Outcome"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 4. 对数坐标图 ---
      p1_log <- ggplot(plot_df, aes(x = time_aligned, y = valuenum, color = Group, fill = Group)) +
        geom_point(aes(shape = Phase), alpha = 0.2, size = 0.6) + 
        geom_smooth(aes(group = LineGroup), method = "gam", size = 1, alpha = 0.3) + 
        geom_line(aes(aes(group = hadm_id ), alpha = 0.2, size = 0.6))+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_y_log10() + 
        scale_color_manual(values = group_colors) +
        scale_fill_manual(values = group_colors) +
        scale_shape_manual(values = c("Pre-Shock" = 1, "Post-Shock" = 16)) +
        labs(
          title = paste(target_biomarker, "(Log10)"),
          subtitle = paste0("Overall Diff (Log) P: ", format_p(p_log_val)),
          x = "Days relative to Shock Onset (T=0)",
          y = "Log10 Scale"
        ) +
        theme_cor +
        theme(legend.position = "none")
      
      # --- 5. 合并与输出 ---
      combined_plot <- (p1_linear | p1_log) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      file_path <- file.path(output_dir, paste0(safe_name, "_CS_Outcome_T0.png"))
      
      ggsave(file_path, plot = combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
  
  message("✅ 心源性休克结局对比图（监护指标）已绘制完毕。")
} 
# 6.实验室指标间相关性图 #--------------------------------------------------------------
# 6.1 住院期间平均值计算相关性 #--------------------------------------------------------------
# 第一次icu病程数据 #---------------------------------------------------
setwd(file.path(Output_derived_data,"data_filtered"))
df_lab_log <- readRDS("icu_filtered_mimic4_Sepsislog.rds")
folder_name <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\correlation\\548第一次住icu的病程记录\\住院期间平均值"
if (!dir.exists(folder_name)) { dir.create(folder_name, recursive = TRUE) }

unique(df_lab_log$label)
colnames(df_lab_log)
df_lab_log1 <- df_lab_log %>% select(hadm_id,charttime_rel_icu,valuenum,label,this_icu_outcome,outtime_rel_icu,deathtime_rel_icu ,is_shock_onset_in_this_icu)


lab_wider <- df_lab_log1 %>%
  group_by(hadm_id, label) %>%
  summarise(avg_value = mean(as.numeric(valuenum) , na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = avg_value)

# 设定阈值
p_threshold <- 0.05      # 显著性水平
r_threshold <- 0.3       # 相关系数绝对值阈值（可选，0.3以上通常认为有相关性）

target_labels <- unique(df_lab_log$label)
 
# 循环开始
for (c in 1:length(target_labels)) {
  start_e <- c + 1 
  if (start_e > length(target_labels)) next 
  
  for (e in start_e:length(target_labels)) {
    var_x <- target_labels[c]
    var_y <- target_labels[e]
    
    # 1. 【核心：计算相关性】
    # 提取两列数据并去掉缺失值
    x_data <- lab_wider[[var_x]]
    y_data <- lab_wider[[var_y]]
    
    # 至少需要一些非空数据点才能计算
    if (sum(!is.na(x_data) & !is.na(y_data)) < 5) {
      message(paste("跳过:", var_x, "vs", var_y, "(有效数据点太少)"))
      next
    }
    
    # 使用 tryCatch 防止因数据问题（如标准差为0）导致报错中断循环
    test_result <- tryCatch({
      cor.test(x_data, y_data, method = "pearson", use = "complete.obs")
    }, error = function(err) return(NULL))
    
    if (is.null(test_result)) next
    
    r_val <- test_result$estimate
    p_val <- test_result$p.value
    
    
    # 如果满足显著性条件，在控制台打印醒目的提示
    if (r_val > r_threshold) {
      # 使用 repeat 符号让提示更显眼
      cat("\n", rep("!", 20), "\n")
      message(paste0("发现显著相关性!!! [ ", var_x, " ] vs [ ", var_y, " ]"))
      message(paste0("R = ", round(r_val, 3), " | P-value = ", format.pval(p_val, digits = 3)))
      cat(rep("!", 20), "\n\n")
      
      # 可选：如果相关，给文件名加个前缀，方便在文件夹里一眼看到
      file_prefix <- "SIG_"
      # 3. 【绘图与保存】
      p3 <- ggplot(data = lab_wider, mapping = aes(x = .data[[var_x]], y = .data[[var_y]])) +
        geom_point(color = "blue", alpha = 0.3, size = 1.5) +
        geom_smooth(color = "red", method = "lm", size = 0.9, formula = y ~ x) +
        stat_cor(method = "pearson", color = "black", digits = 3) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10, face = "bold", color = "black"),
              axis.ticks = element_line(color = "black", linewidth = 1),
              axis.ticks.length = unit(0.15, "cm"),
              legend.position = "none") + 
        xlab(var_x) +
        ylab(var_y)
      
      # 保存图片
      safe_name_x <- gsub("[^a-zA-Z0-9]", "_", var_x)
      safe_name_y <- gsub("[^a-zA-Z0-9]", "_", var_y)
      file_name <- paste0(file_prefix, safe_name_x, "_vs_", safe_name_y, "_corr.png")
      
      ggsave(filename = file_name, plot = p3, path = folder_name, 
             width = 6, height = 6, units = "in", dpi = 300)
      
    } else {
      file_prefix <- ""
    }
    
    
  }
}

# 热图
{
  {
    library(corrplot)
    library(dplyr)
    
    # 1. 准备数据
    cor_data <- lab_wider %>% select(-hadm_id)
    M <- cor(cor_data, method = "pearson", use = "pairwise.complete.obs")
    
    # 2. 聚类排序
    ord <- corrMatOrder(M, order = "hclust")
    M_ord <- M[ord, ord]
    n <- ncol(M_ord)
    
    # 3. 颜色与保存设置
    col_palette <- colorRampPalette(c("#377EB8", "white", "#E41A1C"))(200)
    label_color <- "#FF8C00" # 橙色
    
    file_path <- file.path(folder_name, "corrplot_Staircase_Final_Fixed.png")
    # 增加图片宽度到 3200，给左侧留位置
    png(file_path, width = 3200, height = 3000, res = 300)
    
    # 【修正点 1】：加大左边距 (第二个参数) 和顶边距 (第三个参数)
    # c(bottom, left, top, right)
    par(mar = c(2, 15, 15, 2)) 
    
    # A. 绘制基础热图
    # 【修正点 2】：使用 xlim 参数！
    # 默认 x 轴是 0.5 到 n+0.5。我们把它改成 -12 到 n+0.5，
    # 这样左边就有了一大块可以写字的“画布区域”。
    corrplot(M_ord, 
             method = "color", 
             type = "upper",            
             col = col_palette,
             addCoef.col = "black",     
             number.cex = 0.45,         
             tl.pos = "n",              
             diag = TRUE,               
             outline = "grey90",        
             cl.pos = "r",
             xlim = c(-12, n + 0.5)     # 核心：撑开左侧坐标轴
    )
    
    # 【关键】：开启裁剪忽略模式，允许在绘图区外写字
    par(xpd = TRUE)
    
    var_names <- colnames(M_ord)
    
    # B. 手动添加顶部标签 (X轴)
    text(x = 1:n, 
         y = n + 1.5,                  # 向上偏移一点，防止压线
         labels = var_names, 
         srt = 90, 
         adj = 0, 
         cex = 0.55, 
         col = label_color)            
    
    # C. 【修正点 3】：手动添加 Y 轴台阶标签
    for (i in 1:n) {
      # 这里的 x 坐标现在可以在 -12 到 1 之间自由浮动了
      # 调大偏移量（如 -0.8），确保文字不会碰到方块
      text(x = i - 0.8, 
           y = n - i + 1, 
           labels = var_names[i], 
           col = label_color, 
           adj = 1,                    # 右对齐：文字末尾对齐方块左沿
           cex = 0.55, 
           font = 2)
    }
    
    dev.off()
    cat("✅ 最终版热图已生成。如果左边还缺，请继续调大 xlim 的负值：\n", file_path)
  }
}
 
# 住院天数大于2小于40的患者log #---------------------------------------------------
folder_name <- "D:\\Lab_project\\2026work\\Markfu\\PLOT\\sepsis\\correlation\\住院天数大于2小于40的患者脓毒症相关指标相关性\\住院期间平均值log"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  print("文件夹创建成功")
} else {
  print("文件夹已存在")
}

# 找出至少测过这些生物标志物之一的患者
{
  target_heart_markers <- c(
    "NTproBNP", 
    # "proBNP, Pleural", 
    "Troponin T", 
    # "CK-MB Index", 
    "Creatine Kinase, MB Isoenzyme" 
    # "Creatine Kinase (CK)"
  )
  
  valid_ids <- df_lab %>%
    filter(label %in% target_heart_markers) %>%
    pull(hadm_id) %>%
    unique()
  
  # 3. 统计一下符合条件的患者数量
  cat("至少测过其中一项心脏标志物的患者(hadm_id)总数为：", length(valid_ids), "\n")
  
  # 4. 提取这些患者的所有化验记录 (为了后续分析)
  df_heart_patients <- df_lab %>%
    filter(hadm_id %in% valid_ids)
  
  unique(df_heart_patients[which(df_heart_patients$label == "Creatine Kinase, MB Isoenzyme"), "valueuom"])
  # ng/mL
  unique(df_heart_patients$valueuom[df_heart_patients$label == "Creatine Kinase, MB Isoenzyme"])
  
  df_heart_patients %>% group_by(shock) %>% summarise(n_units = n_distinct(hadm_id))
  # shock n_units
  # <dbl>   <int>
  #   1                 0    6076
  # 2                 1     255
  
  df_heart_patients1 <- df_heart_patients %>% select(hadm_id,TIME,valuenum,label,Outcome_Mortality,Length_of_Stay,In_hospital_Mortality_Time,shock)
  heart_patients_wider <- df_heart_patients1 %>%
    group_by(hadm_id, label) %>%
    summarise(avg_value = mean(as.numeric(valuenum) , na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = label, values_from = avg_value)
}
colnames(df_lab_log)
df_lab1 <- df_lab_log %>% select(hadm_id,TIME,log_valuenum,label,Outcome_Mortality,Length_of_Stay,In_hospital_Mortality_Time,shock)

lab_wider <- df_lab1 %>%
  group_by(hadm_id, label) %>%
  summarise(avg_value = mean(as.numeric(log_valuenum) , na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = avg_value)


# 设定阈值
p_threshold <- 0.05      # 显著性水平
r_threshold <- 0.3       # 相关系数绝对值阈值（可选，0.3以上通常认为有相关性）

# 循环开始
for (c in 1:length(target_labels)) {
  start_e <- c + 1 
  if (start_e > length(target_labels)) next 
  
  for (e in start_e:length(target_labels)) {
    var_x <- target_labels[c]
    var_y <- target_labels[e]
    
    # 1. 【核心：计算相关性】
    # 提取两列数据并去掉缺失值
    x_data <- lab_wider[[var_x]]
    y_data <- lab_wider[[var_y]]
    
    # 至少需要一些非空数据点才能计算
    if (sum(!is.na(x_data) & !is.na(y_data)) < 5) {
      message(paste("跳过:", var_x, "vs", var_y, "(有效数据点太少)"))
      next
    }
    # “使用线性混合效应模型（LMM）分析指标的组间差异，将‘组别’和‘进入ICU后的相对时间’作为固定效应，将‘患者ID’作为随机截距以修正重复测量数据的内相关性。通过对固定效应 Group 的 $t$ 检验来评估两组间的显著性差异。”
    # 使用 tryCatch 防止因数据问题（如标准差为0）导致报错中断循环
    test_result <- tryCatch({
      cor.test(x_data, y_data, method = "pearson", use = "complete.obs")
    }, error = function(err) return(NULL))
    
    if (is.null(test_result)) next
    
    r_val <- test_result$estimate
    p_val <- test_result$p.value
    
    
    # 如果满足显著性条件，在控制台打印醒目的提示
    if (r_val > r_threshold) {
      # 使用 repeat 符号让提示更显眼
      cat("\n", rep("!", 20), "\n")
      message(paste0("发现显著相关性!!! [ ", var_x, " ] vs [ ", var_y, " ]"))
      message(paste0("R = ", round(r_val, 3), " | P-value = ", format.pval(p_val, digits = 3)))
      cat(rep("!", 20), "\n\n")
      
      # 可选：如果相关，给文件名加个前缀，方便在文件夹里一眼看到
      file_prefix <- "SIG_"
      # 3. 【绘图与保存】
      p3 <- ggplot(data = lab_wider, mapping = aes(x = .data[[var_x]], y = .data[[var_y]])) +
        geom_point(color = "blue", alpha = 0.3, size = 1.5) +
        geom_smooth(color = "red", method = "lm", size = 0.9, formula = y ~ x) +
        stat_cor(method = "pearson", color = "black", digits = 3) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10, face = "bold", color = "black"),
              axis.ticks = element_line(color = "black", linewidth = 1),
              axis.ticks.length = unit(0.15, "cm"),
              legend.position = "none") + 
        xlab(var_x) +
        ylab(var_y)
      
      # 保存图片
      safe_name_x <- gsub("[^a-zA-Z0-9]", "_", var_x)
      safe_name_y <- gsub("[^a-zA-Z0-9]", "_", var_y)
      file_name <- paste0(file_prefix, safe_name_x, "_vs_", safe_name_y, "_corr.png")
      
      ggsave(filename = file_name, plot = p3, path = folder_name, 
             width = 6, height = 6, units = "in", dpi = 300)
      
    } else {
      file_prefix <- ""
    }
    
    
  }
}


# 住院天数大于2小于40的中心源性休克患者 #---------------------------------------------------
folder_name <- "D:\\Lab_project\\2026work\\Markfu\\PLOT\\sepsis\\correlation\\住院天数大于2小于40的患者脓毒症相关指标相关性\\心源性休克患者住院期间平均值"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  print("文件夹创建成功")
} else {
  print("文件夹已存在")
}


lab_wider1 <- df_lab1 %>%
  filter(shock == 1) %>% 
  group_by(hadm_id, label) %>%
  summarise(avg_value = mean(as.numeric(valuenum) , na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = avg_value)


# 设定阈值
p_threshold <- 0.05      # 显著性水平
r_threshold <- 0.3       # 相关系数绝对值阈值（可选，0.3以上通常认为有相关性）

# 循环开始
for (c in 1:length(target_labels)) {
  start_e <- c + 1 
  if (start_e > length(target_labels)) next 
  
  for (e in start_e:length(target_labels)) {
    var_x <- target_labels[c]
    var_y <- target_labels[e]
    
    # 1. 【核心：计算相关性】
    # 提取两列数据并去掉缺失值
    x_data <- lab_wider1[[var_x]]
    y_data <- lab_wider1[[var_y]]
    
    # 至少需要一些非空数据点才能计算
    if (sum(!is.na(x_data) & !is.na(y_data)) < 5) {
      message(paste("跳过:", var_x, "vs", var_y, "(有效数据点太少)"))
      next
    }
    
    # 使用 tryCatch 防止因数据问题（如标准差为0）导致报错中断循环
    test_result <- tryCatch({
      cor.test(x_data, y_data, method = "pearson", use = "complete.obs")
    }, error = function(err) return(NULL))
    
    if (is.null(test_result)) next
    
    r_val <- test_result$estimate
    p_val <- test_result$p.value
    
    
    # 如果满足显著性条件，在控制台打印醒目的提示
    if (r_val > r_threshold) {
      # 使用 repeat 符号让提示更显眼
      cat("\n", rep("!", 20), "\n")
      message(paste0("发现显著相关性!!! [ ", var_x, " ] vs [ ", var_y, " ]"))
      message(paste0("R = ", round(r_val, 3), " | P-value = ", format.pval(p_val, digits = 3)))
      cat(rep("!", 20), "\n\n")
      
      # 可选：如果相关，给文件名加个前缀，方便在文件夹里一眼看到
      file_prefix <- "SIG_"
      # 3. 【绘图与保存】
      p3 <- ggplot(data = lab_wider1, mapping = aes(x = .data[[var_x]], y = .data[[var_y]])) +
        geom_point(color = "blue", alpha = 0.3, size = 1.5) +
        geom_smooth(color = "red", method = "lm", size = 0.9, formula = y ~ x) +
        stat_cor(method = "pearson", color = "black", digits = 3) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10, face = "bold", color = "black"),
              axis.ticks = element_line(color = "black", linewidth = 1),
              axis.ticks.length = unit(0.15, "cm"),
              legend.position = "none") + 
        xlab(var_x) +
        ylab(var_y)
      
      # 保存图片
      safe_name_x <- gsub("[^a-zA-Z0-9]", "_", var_x)
      safe_name_y <- gsub("[^a-zA-Z0-9]", "_", var_y)
      file_name <- paste0(file_prefix, safe_name_x, "_vs_", safe_name_y, "_corr.png")
      
      ggsave(filename = file_name, plot = p3, path = folder_name, 
             width = 6, height = 6, units = "in", dpi = 300)
      
    } else {
      file_prefix <- ""
    }
    
    
  }
}

# 住院天数大于2小于40的中心源性休克患者log #---------------------------------------------------
folder_name <- "D:\\Lab_project\\2026work\\Markfu\\PLOT\\sepsis\\correlation\\住院天数大于2小于40的患者脓毒症相关指标相关性\\心源性休克患者住院期间平均值log"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
  print("文件夹创建成功")
} else {
  print("文件夹已存在")
}


lab_wider1 <- df_lab_log %>%
  filter(shock == 1) %>% 
  group_by(hadm_id, label) %>%
  summarise(avg_value = mean(as.numeric(log_valuenum) , na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = avg_value)


# 设定阈值
p_threshold <- 0.05      # 显著性水平
r_threshold <- 0.3       # 相关系数绝对值阈值（可选，0.3以上通常认为有相关性）

# 循环开始
for (c in 1:length(target_labels)) {
  start_e <- c + 1 
  if (start_e > length(target_labels)) next 
  
  for (e in start_e:length(target_labels)) {
    var_x <- target_labels[c]
    var_y <- target_labels[e]
    
    # 1. 【核心：计算相关性】
    # 提取两列数据并去掉缺失值
    x_data <- lab_wider1[[var_x]]
    y_data <- lab_wider1[[var_y]]
    
    # 至少需要一些非空数据点才能计算
    if (sum(!is.na(x_data) & !is.na(y_data)) < 5) {
      message(paste("跳过:", var_x, "vs", var_y, "(有效数据点太少)"))
      next
    }
    
    # 使用 tryCatch 防止因数据问题（如标准差为0）导致报错中断循环
    test_result <- tryCatch({
      cor.test(x_data, y_data, method = "pearson", use = "complete.obs")
    }, error = function(err) return(NULL))
    
    if (is.null(test_result)) next
    
    r_val <- test_result$estimate
    p_val <- test_result$p.value
    
    
    # 如果满足显著性条件，在控制台打印醒目的提示
    if (r_val > r_threshold) {
      # 使用 repeat 符号让提示更显眼
      cat("\n", rep("!", 20), "\n")
      message(paste0("发现显著相关性!!! [ ", var_x, " ] vs [ ", var_y, " ]"))
      message(paste0("R = ", round(r_val, 3), " | P-value = ", format.pval(p_val, digits = 3)))
      cat(rep("!", 20), "\n\n")
      
      # 可选：如果相关，给文件名加个前缀，方便在文件夹里一眼看到
      file_prefix <- "SIG_"
      # 3. 【绘图与保存】
      p3 <- ggplot(data = lab_wider1, mapping = aes(x = .data[[var_x]], y = .data[[var_y]])) +
        geom_point(color = "blue", alpha = 0.3, size = 1.5) +
        geom_smooth(color = "red", method = "lm", size = 0.9, formula = y ~ x) +
        stat_cor(method = "pearson", color = "black", digits = 3) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
              axis.title = element_text(size = 12, face = "bold"),
              axis.text = element_text(size = 10, face = "bold", color = "black"),
              axis.ticks = element_line(color = "black", linewidth = 1),
              axis.ticks.length = unit(0.15, "cm"),
              legend.position = "none") + 
        xlab(var_x) +
        ylab(var_y)
      
      # 保存图片
      safe_name_x <- gsub("[^a-zA-Z0-9]", "_", var_x)
      safe_name_y <- gsub("[^a-zA-Z0-9]", "_", var_y)
      file_name <- paste0(file_prefix, safe_name_x, "_vs_", safe_name_y, "_corr.png")
      
      ggsave(filename = file_name, plot = p3, path = folder_name, 
             width = 6, height = 6, units = "in", dpi = 300)
      
    } else {
      file_prefix <- ""
    }
    
    
  }
}

# 7.实验室指标基线值_人口学信息相关性分析 #-------------------------------
{
  # 第一次icu病程数据 #---------------------------------------------------
  setwd(file.path(Output_derived_data,"data_filtered"))
  df_lab_log <- readRDS("icu_filtered_mimic4_Sepsislog.rds")
  folder_name <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\correlation\\548第一次住icu的病程记录\\基线值"
  if (!dir.exists(folder_name)) { dir.create(folder_name, recursive = TRUE) }
  
  unique(df_lab_log$label)
  colnames(df_lab_log)
  df_lab_log1 <- df_lab_log %>% mutate(race_group = case_when(
    grepl("WHITE|PORTUGUESE", race, ignore.case = TRUE) ~ 1,
    grepl("BLACK", race, ignore.case = TRUE) ~ 2,
    grepl("ASIAN", race, ignore.case = TRUE) ~ 3,
    TRUE ~ 4
  )) %>% 
  select(hadm_id,gender,age,race_group,charttime_rel_icu,valuenum,label)
  lab_wider_baseline <- df_lab_log1 %>%
    group_by(hadm_id, label) %>%
    arrange(charttime_rel_icu) %>%
    slice(1) %>%
    ungroup() %>% 
    select(-charttime_rel_icu) %>% 
    pivot_wider(names_from = label, values_from = valuenum)
  colnames(lab_wider_baseline)
 
  covname <- colnames(lab_wider_baseline[,c("age"                            ,                  
                                            "Anion Gap"                      ,"Base Excess"                      ,               
                                            "Calculated Total CO2"           ,"Chloride"                       ,
                                             "Creatine Kinase (CK)"          , "Creatine Kinase, MB Isoenzyme"  ,
                                             "Creatinine"                    , "Glucose"                        ,
                                             "Hematocrit"                    , "Hemoglobin"                     ,
                                             "INR(PT)"                       , "Lactate"                        ,
                                             "Lymphocytes"                   , "Magnesium"                      ,
                                             "Monocytes"                     , "Neutrophils"                    ,
                                             "PT"                            , "PTT"                            ,
                                             "Phosphate"                     , "Platelet Count"                 ,
                                             "Potassium"                     , "RDW"                            ,
                                             "Red Blood Cells"               , "Sodium"                         ,
                                             "Troponin T"                    , "Urea Nitrogen"                  ,
                                             "White Blood Cells"             , "pCO2"                           ,
                                             "pH"                            , "pO2"                            ,
                                             "Alanine Aminotransferase (ALT)", "Asparate Aminotransferase (AST)",
                                             "Bilirubin, Total"              , "Free Calcium"                   ,
                                             "Lactate Dehydrogenase (LD)"    , "Albumin"   )])
  
  catcovname <- colnames(lab_wider_baseline[,c("gender","race_group")]) 
  cov <- lab_wider_baseline %>%  mutate(across(all_of(covname), as.numeric), across(all_of(catcovname), as.factor))

  # 连续协变量之间是否存在相关性预筛选
  for (c in 1:length(covname)) {
    if (c < length(covname)){
      d <- c+1
    }else{
      d <- c
    }
    for (e in d:length(covname)) {
      p3 <- ggplot(data=cov, mapping = aes(x=get((covname[c])),y =get((covname[e]))))+
        geom_point(color="blue",alpha=0.3,size=1.5)+
        geom_smooth(color="red",method = "lm",size=0.9)+
         ggpubr::stat_cor(method = "pearson",color="black",digits = 3)+ #添加相关性系数 根据画图需要选择是否保留
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth =1.2),
              axis.title = element_text(size = 12,face = "bold"),
              axis.title.x = element_text(vjust = -0.8),
              axis.title.y = element_text(vjust = 1.5),
              axis.text = element_text(size=10,face = "bold",color = "black"),
              axis.ticks = element_line(color="black",linewidth =1),
              axis.ticks.length=unit(0.15,"cm"),
              legend.position = "none") + 
        xlab(paste(covname[c]))+
        ylab(paste(covname[e]))
      setwd(folder_name)
      png(filename = paste(covname[c],covname[e],"corr.png",sep = "_"),
          width = 800,height = 800,res = 240)
      print(p3)
      dev.off()
    }
  }
  # 分类协变量与连续协变量分析
  for (j in 1:length(catcovname)) {
    for (k in 1:length(covname)) {
      p4 <- ggplot(data=cov,mapping = aes(x=get((catcovname[j])),y =get((covname[k]))))+
        geom_boxplot()+
         ggpubr::stat_compare_means(method = "t.test",aes(label=paste0("p=",..p.format..)))+ #添加显著性 按需要删减
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, linewidth =1.2),
              axis.title = element_text(size = 12,face = "bold"),
              axis.title.x = element_text(vjust = -0.8),
              axis.title.y = element_text(vjust = 1.5),
              axis.text = element_text(size=10,face = "bold",color = "black"),
              axis.ticks = element_line(color="black",linewidth =1),
              axis.ticks.length=unit(0.15,"cm"),
              legend.position = "none") + 
        xlab(paste(catcovname[j]))+
        ylab(paste(covname[k]))
      setwd(folder_name)
      png(filename = paste(catcovname[j],covname[k],"corr.png",sep = "_"),
          width = 800,height = 800,res = 240)
      print(p4)
      dev.off()
    }
  }
  ### 绘制相关性总图
  ggparis_data <- cov[,c(2:42)] #根据数据集修改
  ggparis_data$GENDER <- as.factor(ggparis_data$gender) #修改分类协变量的数据格式
  ggparis_data$ETHNIC <- as.factor(ggparis_data$race_group) #修改分类协变量的数据格式 
  
  p5 <- ggpairs(ggparis_data,
                lower = list(continuous = wrap(lowerFn)),
                upper = list(continuous = wrap("cor",size = 3,digits = 3),
                             combo = wrap("box",outlier.size = 0.3,size = 0.3)),
                columnLabels = colnames(ggparis_data)) + 
    theme(axis.text = element_text(size = 5),
          axis.title = element_text(size = 5))
  setwd(folder_name)
  png(filename = paste("ggparis_corr.png"),
      width = 6000,height = 6000,res = 240)
  print(p5)
  dev.off()
}

# 对race_group四个分组的相关性分析
{
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(rstatix) 
  
  # --- 1. 定义主题与配色 ---
  mytheme <- theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1.2),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", face = "bold"),
    legend.position = "none"
  )

  # 确保输出文件夹存在
  output_path <- file.path(folder_name, "T_Test_Trend_Plots")
  if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
  
  # --- 2. 开始自动化绘图循环 ---
  for (cont_var in covname) {
    cat("正在生成 T 检验趋势图: ", cont_var, "\n")
    
    # A. 数据清洗：剔除 NA
    plot_df <- cov %>% 
      select(race_group, all_of(cont_var)) %>% 
      filter(!is.na(.data[[cont_var]]), !is.na(race_group))
    
    # 样本量太小时跳过
    if (nrow(plot_df) < 10) next
    
    tryCatch({
      # B. 计算两组间的显著性 (核心修改：t_test)
      # 自动执行 race_group 内各组之间的两两 T 检验
      stat_test <- plot_df %>%
        t_test(as.formula(paste0("`", cont_var, "` ~ race_group"))) %>%
        add_significance() %>%
        add_xy_position(x = "race_group") 
      # %>%
        # 可选：如果你只想看显著的对比，取消下面这行的注释
        # filter(p < 0.05) 
      
      # C. 绘图
      p <- ggplot(plot_df, aes(x = race_group, y = .data[[cont_var]])) +
        # 1. 箱线图上下误差线
        stat_boxplot(geom = "errorbar", width = 0.2, linewidth = 0.8) +
        
        # 2. 箱线图本体
        geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", color = "black") +
        
        # 3. 抖动散点 (pch=21 为带边框的圆点)
        geom_jitter(aes(fill = race_group), pch = 21, size = 3, width = 0.2, alpha = 0.7) +
        
        # 4. 跨组线性趋势线 (group=1 强行跨类)
        geom_smooth(method = "lm", formula = y ~ x, 
                    aes(group = 1), 
                    color = "black", linetype = "dashed", size = 1, se = TRUE) +
        
        # 5. 添加 T 检验显著性标记 (stat_pvalue_manual)
        stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
        
        # 6. 添加相关性统计 (R2 和 P)
        stat_cor(aes(group = 1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
                 method = "pearson", 
                 label.x.npc = "left", 
                 label.y.npc = 1, # 锁定在绘图区顶部
                 size = 4.5, fontface = "bold") +
        
        # 7. 添加回归方程
        stat_regline_equation(aes(group = 1), 
                              label.x.npc = "left", 
                              label.y.npc = 0.95, # 稍微偏下一点
                              size = 4.5, fontface = "bold") +
        
        # 8. 细节修饰
        scale_fill_manual(values = c("#db6968", "#4d97cd", "#f8984e", "#459943")) +
        # 关键：扩展 Y 轴顶部空间 (30%)，防止 T 检验连线被遮挡
        scale_y_continuous(expand = expansion(mult = c(0.1, 0.5))) + 
        labs(x = "Race Group", y = cont_var) +
        mytheme
      
      output_path <- folder_name
      # D. 保存
      safe_name <- gsub("[[:punct:]]| ", "_", cont_var)
      ggsave(filename = file.path(output_path, paste0("TTestTrend_", safe_name, ".png")),
             plot = p, width = 8, height = 8, dpi = 300)
      
    }, error = function(e) {
      cat("  ❌ ", cont_var, " 出错: ", e$message, "\n")
    })
  }
  
  cat("\n✅ 所有 T 检验趋势图已导出至:", output_path, "\n")
}
 
# 输出有显著性的表
{
  library(dplyr)
  library(tidyr)
  
  # --- 1. 设定阈值 ---
  p_threshold <- 0.05
  r_threshold <- 0.3  # 仅针对连续变量相关性
  
  cont_results <- list()
  cat_results <- list()
  
  # --- 2. 连续变量 vs 连续变量 (保持 cor.test，因为 t.test 不适用) ---
  cat("正在计算连续变量相关性 (Pearson)...\n")
  cont_idx <- 1
  for (i in 1:(length(covname) - 1)) {
    for (j in (i + 1):length(covname)) {
      var1 <- covname[i]
      var2 <- covname[j]
      
      tmp_data <- cov %>% select(all_of(c(var1, var2))) %>% na.omit()
      
      if (nrow(tmp_data) > 10) {
        try({
          test <- cor.test(tmp_data[[1]], tmp_data[[2]], method = "pearson")
          if (test$p.value < p_threshold && abs(test$estimate) >= r_threshold) {
            cont_results[[cont_idx]] <- data.frame(
              Variable_A = var1,
              Variable_B = var2,
              R_Value = round(test$estimate, 4),
              P_Value = test$p.value,
              Method = "Pearson Correlation"
            )
            cont_idx <- cont_idx + 1
          }
        }, silent = TRUE)
      }
    }
  }
  
  # --- 3. 分类变量 vs 连续变量 (使用 t.test) ---
  cat("正在使用 t.test 计算分类与连续变量相关性...\n")
  cat_idx <- 1
  for (i in 1:length(catcovname)) {
    for (j in 1:length(covname)) {
      cat_var <- catcovname[i]
      cont_var <- covname[j]
      
      # 准备数据
      tmp_data <- cov %>% select(all_of(c(cat_var, cont_var))) %>% na.omit()
      
      # 检查分类变量的水平数量
      levels_count <- length(unique(tmp_data[[cat_var]]))
      
      if (nrow(tmp_data) > 10) {
        if (levels_count == 2) {
          # 仅当分类变量恰好有两组时（如性别），适用 t.test
          try({
            form <- as.formula(paste0("`", cont_var, "` ~ `", cat_var, "`"))
            t_res <- t.test(form, data = tmp_data)
            
            if (t_res$p.value < p_threshold) {
              cat_results[[cat_idx]] <- data.frame(
                Categorical_Var = cat_var,
                Continuous_Var = cont_var,
                P_Value = t_res$p.value,
                Method = "t.test",
                Note = "Groups: 2"
              )
              cat_idx <- cat_idx + 1
            }
          }, silent = TRUE)
        } else if (levels_count > 2) {
          # 如果超过两组（如 Race），t.test 不适用
          # 此处仅记录，不进行 t.test 运算以防报错
          # message(paste("跳过", cat_var, ": 组数 > 2，不适用 t.test"))
        }
      }
    }
  }
  
  # --- 4. 汇总并导出 ---
  
  # 转化为数据框
  df_continuous_corr <- bind_rows(cont_results)
  df_categorical_ttest <- bind_rows(cat_results)
  
  # 排序
  if(nrow(df_continuous_corr) > 0) {
    df_continuous_corr <- df_continuous_corr %>% arrange(P_Value)
  }
  if(nrow(df_categorical_ttest) > 0) {
    df_categorical_ttest <- df_categorical_ttest %>% arrange(P_Value)
  }
  
  # 打印反馈
  print("--- 连续变量相关性 (P < 0.05) ---")
  print(head(df_continuous_corr))
  
  print("--- 分类变量 t.test 显著差异 (P < 0.05) ---")
  print(head(df_categorical_ttest))
  
  # 保存
  write.csv(df_continuous_corr, "Corr_Continuous_Pearson.csv", row.names = FALSE)
  write.csv(df_categorical_ttest, "Corr_Categorical_TTest.csv", row.names = FALSE)
  
  cat("\n✅ 分析完成！已生成：\n 1. Corr_Continuous_Pearson.csv\n 2. Corr_Categorical_TTest.csv\n")
}
# 输出有显著性的图
# 输出有显著性的图 (连续变量 Pearson / 性别 t.test)
{
  library(ggplot2)
  library(dplyr)
  library(ggpubr) # 用于方便地添加显著性标注
  
  # --- 设定阈值 ---
  p_threshold <- 0.05      # 显著性水平
  r_threshold <- 0.3       # 相关系数绝对值阈值
 
  folder_name <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\correlation\\548第一次住icu的病程记录\\基线值\\有显著相关性"
  # 确保输出目录存在
  if(!dir.exists(folder_name)) dir.create(folder_name, recursive = TRUE)
  setwd(folder_name)
  
  # 1. 连续协变量之间：基于 P 值和 R 值的筛选 (Pearson)
  cat("开始分析连续变量相关性...\n")
  for (c in 1:(length(covname) - 1)) {
    for (e in (c + 1):length(covname)) {
      
      # 提取数据并剔除 NA
      test_data <- cov %>% select(all_of(c(covname[c], covname[e]))) %>% na.omit()
      
      if(nrow(test_data) < 5) next
      
      tryCatch({
        ctest <- cor.test(test_data[[1]], test_data[[2]], method = "pearson")
        r_val <- ctest$estimate
        p_val <- ctest$p.value
        
        # 阈值判定
        if (p_val < p_threshold && abs(r_val) >= r_threshold) {
          
          p3 <- ggplot(data = cov, mapping = aes(x = .data[[covname[c]]], y = .data[[covname[e]]])) +
            geom_point(color = "#80B1D3", alpha = 0.8, size = 1.5) +
            ggpubr::stat_cor(method = "pearson",color="black",digits = 3)+ #添加相关性系数 根据画图需要选择是否保留
            geom_smooth(color = "#FB8072", method = "lm", size = 0.9) +
            labs(
                 # title = paste0("Pearson R=", round(r_val, 3), ", P=", format.pval(p_val, digits = 3)),
                 x = covname[c], y = covname[e]) +
            theme_bw() + 
            theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.2),
                  axis.title = element_text(size = 12, face = "bold"),
                  axis.text = element_text(size = 10, face = "bold", color = "black"))
          
          safe_name_c <- gsub("[[:punct:]]| ", "_", covname[c])
          safe_name_e <- gsub("[[:punct:]]| ", "_", covname[e])
          
          setwd(folder_name)
          png(filename = paste("SIG_CONT", safe_name_c, safe_name_e, "corr.png", sep = "_"),
              width = 1000, height = 1000, res = 200)
          print(p3)
          dev.off()
          cat("  [已保存] 连续相关:", covname[c], "vs", covname[e], "\n")
        }
      }, error = function(e) { message("跳过: ", covname[c], " ", covname[e]) })
    }
  }
  
  # 2. 性别 (Gender) vs 连续协变量：基于 t.test P 值的筛选
  cat("\n开始分析性别差异 (t.test)...\n")
  # 目标变量固定为 gender
  target_cat <- "gender" 
  
  for (k in 1:length(covname)) {
    cont_var <- covname[k]
    
    # 提取数据并剔除 NA
    test_data <- cov %>% select(all_of(c(target_cat, cont_var))) %>% na.omit()
    
    # 确保两组都有数据且样本量足够
    if(length(unique(test_data[[target_cat]])) < 2 | nrow(test_data) < 10) next
    
    tryCatch({
      # 执行 t 检验
      t_res <- t.test(as.formula(paste0("`", cont_var, "` ~ ", target_cat)), data = test_data)
      p_val <- t_res$p.value
      
      # 阈值判定
      if (!is.na(p_val) && p_val < p_threshold) {
        
        p4 <- ggplot(data = test_data, aes(x = .data[[target_cat]], y = .data[[cont_var]])) +
          geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", color = "black") +
          # 添加抖动点增加可视化信息
          geom_jitter(aes(fill = gender), pch = 21, size = 2, width = 0.2, alpha = 0.7) +
          ggpubr::stat_compare_means(method = "t.test",aes(label=paste0("p=",..p.format..)))+ #添加显著性 按需要删减
          # 添加显著性标记
          # stat_compare_means(method = "t.test", label = "p.signif", label.x = 1.5) +
          labs(
            # title = paste0("t-test P=", format.pval(p_val, digits = 3)),
            #    subtitle = paste("Comparison of", cont_var, "by Gender"),
               x = "Gender", y = cont_var) +
          theme_bw() +
          scale_fill_manual(values = c("F" = "#FB8072", "M" = "#80B1D3")) + # 建议配色
          theme(panel.border = element_rect(colour = "black", fill=NA, linewidth = 1.2),
                axis.title = element_text(size = 12, face = "bold"),
                axis.text = element_text(size = 10, face = "bold", color = "black"),
                legend.position = "none")
        
        safe_name_k <- gsub("[[:punct:]]| ", "_", cont_var)
        
        setwd(folder_name)
        png(filename = paste("SIG_GENDER", safe_name_k, "ttest.png", sep = "_"),
            width = 1000, height = 1000, res = 200)
        print(p4)
        dev.off()
        cat("  [已保存] 性别显著差异:", cont_var, "\n")
      }
    }, error = function(e) { message("跳过性别分析: ", cont_var) })
  }
}
unique(test_data$gender)
# 7.实验室指标箱型图 # ---------------------------------------------------------
# 住院天数大于2小于40的患者 #---------------------------------------------------
{
  folder_name <- "D:\\Lab_project\\2026work\\Markfu\\PLOT\\sepsis\\correlation\\住院天数大于2小于40的患者脓毒症相关指标相关性\\lab_shock箱型图"
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
    print("文件夹创建成功")
  } else {
    print("文件夹已存在")
  }
  
  for (i in 1:length(target_labels)) {
    plot_data <- df_lab %>% mutate(shock = as.factor(shock),
                                   valuenum = as.numeric(valuenum)) %>% filter(!is.na(valuenum)) %>% 
      filter(label == target_labels[i] )
    
    p4 <- ggplot(data=plot_data,mapping = aes(x=shock, y = valuenum, color = shock))+
      geom_boxplot()+
      ggpubr::stat_compare_means(method = "t.test",aes(label=paste0("p=",..p.format..)))+ #添加显著性 按需要删减
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth =1.2),
            axis.title = element_text(size = 12,face = "bold"),
            axis.title.x = element_text(vjust = -0.8),
            axis.title.y = element_text(vjust = 1.5),
            axis.text = element_text(size=10,face = "bold",color = "black"),
            axis.ticks = element_line(color="black",linewidth =1),
            axis.ticks.length=unit(0.15,"cm"),
            legend.position = "none") + 
      xlab("Cardiogenic shock")+
      ylab(paste(target_labels[i]))
    setwd(folder_name)
    png(filename = paste("shock",target_labels[i],"corr.png",sep = "_"),
        width = 800,height = 800,res = 240)
    print(p4)
    dev.off()
  }
  
  
}

# 住院天数大于2小于40的患者log #---------------------------------------------------
log_needed_vars <- c(
  "NTproBNP", "proBNP, Pleural",
  "Asparate Aminotransferase (AST)", 
  "Alanine Aminotransferase (ALT)",
  "Creatine Kinase (CK)", 
  "Bilirubin, Total", "Bilirubin, Direct",
  "Troponin T","Lactate Dehydrogenase (LD)","Lactate" ,
  "Creatine Kinase, MB Isoenzyme",
  "White Blood Cells"
)

{
  folder_name <- "D:\\Lab_project\\2026work\\Markfu\\PLOT\\sepsis\\correlation\\住院天数大于2小于40的患者脓毒症相关指标相关性\\lab_shock箱型图log"
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
    print("文件夹创建成功")
  } else {
    print("文件夹已存在")
  }
  
  for (i in 1:length(log_needed_vars)) {
    plot_data <- df_lab_log %>% mutate(shock = as.factor(shock),
                                       log_valuenum = as.numeric(log_valuenum)) %>% filter(!is.na(log_valuenum)) %>% 
      filter(label == log_needed_vars[i] )
    
    p4 <- ggplot(data=plot_data,mapping = aes(x=shock, y = log_valuenum, color = shock))+
      geom_boxplot()+
      ggpubr::stat_compare_means(method = "t.test",aes(label=paste0("p=",..p.format..)))+ #添加显著性 按需要删减
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth =1.2),
            axis.title = element_text(size = 12,face = "bold"),
            axis.title.x = element_text(vjust = -0.8),
            axis.title.y = element_text(vjust = 1.5),
            axis.text = element_text(size=10,face = "bold",color = "black"),
            axis.ticks = element_line(color="black",linewidth =1),
            axis.ticks.length=unit(0.15,"cm"),
            legend.position = "none") + 
      xlab("Cardiogenic shock")+
      ylab(paste(target_labels[i]))
    setwd(folder_name)
    png(filename = paste("shock",log_needed_vars[i],"corr.png",sep = "_"),
        width = 800,height = 800,res = 240)
    print(p4)
    dev.off()
  }
  
  
}


