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
  library(tidyverse)
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
setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
df_diagnose <- read.xlsx("mimic4_sepsis_seqnum1_msm_诊断.xlsx")
colnames(df_diagnose)

# 查看每个hadm_id的seq_num==1的诊断详情
df_diagnose1 <- df_diagnose %>%
  filter(seq_num == 1) %>%
  select(hadm_id, icd_code, icd_version, long_title)  

# 完整统计：每种诊断的患者数、占比
diagnose_stats <- df_diagnose1 %>%
  group_by(icd_code, long_title) %>%
  summarise(n_patients = n()) %>%
  ungroup() %>%
  mutate(
    percentage = n_patients / sum(n_patients) * 100,
    cumulative_percent = cumsum(percentage)
  ) %>%
  arrange(desc(n_patients))

# 查看结果
print("=== 诊断统计前20名 ===")
print(head(diagnose_stats))

print(paste("=== 总诊断种类数: ", nrow(diagnose_stats), " ==="))
print(paste("=== 总患者（住院）人次: ", sum(diagnose_stats$n_patients), " ==="))

setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
write.xlsx(diagnose_stats,"mimic4_sepsis_seqnum1_msm_脓毒症诊断种类.xlsx")


# 加载必要的包
library(openxlsx)
library(dplyr)
library(stringr)

# 1. 读取数据
df_diagnose <- read.xlsx("mimic4_sepsis_seqnum1_msm_诊断.xlsx")

# 2. 定义心血管疾病相关的关键词 (新增高血压相关)
# 涵盖：心衰、心梗、房颤、心脏/心血管、冠状动脉、缺血、心绞痛 + 高血压/血压升高
cvd_keywords <- "heart failure|myocardial infarction|atrial fibrillation|cardiac|cardiovascular|coronary|ischemic heart|angina|hypertension|hypertensive|high blood pressure"

# 3. 标记每一条诊断记录是否为心血管疾病
df_diagnose <- df_diagnose %>%
  mutate(
    # 将长标题转换为小写后进行匹配，忽略大小写差异
    is_cvd = str_detect(tolower(long_title), cvd_keywords)
  )

# 4. 按 hadm_id 汇总，判断每个患者是否合并心血管疾病
patient_cvd_summary <- df_diagnose %>%
  group_by(hadm_id) %>%
  summarise(
    # 只要该患者的任意一条诊断匹配上了心血管疾病，就记为 TRUE
    has_cvd = any(is_cvd, na.rm = TRUE)
  )

# 5. 打印统计结果
total_patients <- nrow(patient_cvd_summary)
cvd_patients <- sum(patient_cvd_summary$has_cvd)
cvd_rate <- round(cvd_patients / total_patients * 100, 2)

cat("================ 重新筛选统计结果 ================\n")
cat("总入室人次 (hadm_id 数量):", total_patients, "\n")
cat("合并心血管疾病(含高血压)的患者人数:", cvd_patients, "\n")
cat("合并心血管疾病(含高血压)的人群比例:", cvd_rate, "%\n")
cat("==================================================\n\n")

# 6. 查看该队列中最常见的心血管疾病具体是哪些 (Top 15)
# 这次可以多看几个，因为高血压的数量通常非常大
top_cvd_diseases <- df_diagnose %>%
  filter(is_cvd == TRUE) %>%
  count(long_title, sort = TRUE) %>%
  head(15)

cat("队列中最常见的前15种心血管疾病诊断（检查是否成功抓取高血压）：\n")
print(top_cvd_diseases)

# 7. (可选) 将打好标签的患者数据导出来，以备后续 msm 建模使用
# 建议导出时把列名改成比较通用的名字，比如 cvd_comorbidity
patient_cvd_summary <- patient_cvd_summary %>%
  mutate(cvd_comorbidity = ifelse(has_cvd, 1, 0)) %>% # 转化为 1/0 方便建模
  select(hadm_id, cvd_comorbidity)

# write.xlsx(patient_cvd_summary, "patient_cvd_comorbidity.xlsx")


{
  # 第一步：筛选出所有心脏相关的诊断（不限制seq_num）
  heart_diseases_all <- df_diagnose %>%
    mutate(
      # 给每个诊断打上简化的心脏病类别标签
      heart_type = case_when(
        # 心力衰竭
        (icd_version == 9 & icd_code %in% c("428", "4280", "4281", "4282", "4283", "4284")) |
          (icd_version == 10 & substr(icd_code, 1, 3) == "I50") ~ "心衰",
        
        # 心肌梗死
        (icd_version == 9 & substr(icd_code, 1, 3) == "410") |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I21", "I22")) ~ "心梗",
        
        # 冠状动脉疾病（不包括急性心梗）
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("411", "412", "413", "414")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I20", "I23", "I24", "I25")) ~ "冠心病",
        
        # 房颤
        (icd_version == 9 & icd_code == "427.31") |
          (icd_version == 10 & substr(icd_code, 1, 3) == "I48") ~ "房颤",
        
        # 其他心律失常
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("426", "427")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I44", "I45", "I46", "I47", "I49")) ~ "其他心律失常",
        
        # 心脏瓣膜病
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("394", "395", "396", "397", "424")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37", "I38")) ~ "瓣膜病",
        
        # 心肌病
        (icd_version == 9 & substr(icd_code, 1, 3) == "425") |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I42", "I43")) ~ "心肌病",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(heart_type))  # 只保留心脏相关诊断
  
  # 第二步：为每个hadm_id合并心脏病种类
  heart_combination_all <- heart_diseases_all %>%
    group_by(hadm_id) %>%
    summarise(
      heart_types = paste(sort(unique(heart_type)), collapse = "+"),
      n_heart_diseases = n_distinct(heart_type),
      all_icd_codes = paste(unique(icd_code), collapse = ", "),
      diagnoses_detail = paste(unique(long_title), collapse = " | ")
    ) %>%
    arrange(hadm_id)
  
  # 查看结果
  cat("=== 有心脏病诊断的住院数（不限首诊断）===\n")
  cat("住院数:", nrow(heart_combination_all), "\n")
  cat("\n=== 前20条记录示例 ===\n")
  print(head(heart_combination_all, 20))  
  }
df_diagnose_other <- df_diagnose %>% filter(!hadm_id %in% heart_combination_all$hadm_id)
setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
write.xlsx(df_diagnose_other,"mimic4_sepsis_seqnum1_msm_脓毒症诊断_没有筛选出心脏病人群.xlsx")


{
  library(dplyr)
  library(stringr)
  
  # 升级后的核心逻辑：增加高血压、心源性休克、心包疾病等路径
  heart_diseases_all_v2 <- df_diagnose %>%
    mutate(
      heart_type = case_when(
        # 1. 心力衰竭 (HF)
        (icd_version == 9 & substr(icd_code, 1, 3) == "428") |
          (icd_version == 10 & substr(icd_code, 1, 3) == "I50") ~ "心衰",
        
        # 2. 高血压 (Hypertension) - 新增
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("401", "402", "403", "404", "405")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I10", "I11", "I12", "I13", "I15")) ~ "高血压",
        
        # 3. 心肌梗死 (MI)
        (icd_version == 9 & substr(icd_code, 1, 3) == "410") |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I21", "I22")) ~ "心梗",
        
        # 4. 心源性休克 (Cardiogenic Shock) - 新增 [来自 PDF 证据]
        (icd_version == 9 & icd_code == "78551") |
          (icd_version == 10 & icd_code == "R570") ~ "心源性休克",
        
        # 5. 房颤与房扑 (AF/AFL)
        (icd_version == 9 & icd_code == "427.31") |
          (icd_version == 10 & substr(icd_code, 1, 3) == "I48") ~ "房颤",
        
        # 6. 其他心律失常
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("426", "427")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I44", "I45", "I46", "I47", "I49")) ~ "其他心律失常",
        
        # 7. 冠心病 (CAD)
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("411", "412", "413", "414")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I20", "I23", "I24", "I25")) ~ "冠心病",
        
        # 8. 心包疾病 (Pericardial Disease) - 新增 [来自 PDF 证据]
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("420", "423")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I30", "I31", "I32")) ~ "心包疾病",
        
        # 9. 瓣膜病 (Valvular Heart Disease)
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("394", "395", "396", "397", "424")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I05", "I06", "I07", "I08", "I34", "I35", "I36", "I37", "I38")) ~ "瓣膜病",
        
        # 10. 心肌病 (Cardiomyopathy) - 包含应激性心肌病 [来自 PDF 证据]
        (icd_version == 9 & substr(icd_code, 1, 3) == "425") |
          (icd_version == 10 & (substr(icd_code, 1, 3) %in% c("I42", "I43") | icd_code == "I5181")) ~ "心肌病",
        
        # 11. 外周血管与肺心病 - 新增
        (icd_version == 9 & substr(icd_code, 1, 3) %in% c("415", "416", "440", "453")) |
          (icd_version == 10 & substr(icd_code, 1, 3) %in% c("I26", "I27", "I28", "I70", "I82")) ~ "肺心病/血管病",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(heart_type))
  
  # 汇总每个 hadm_id 的心血管病谱
  cvd_summary_final <- heart_diseases_all_v2 %>%
    group_by(hadm_id) %>%
    summarise(
      cvd_types = paste(sort(unique(heart_type)), collapse = "+"),
      # 是否存在至少一种严重疾病 (不含单纯高血压)
      has_major_heart_disease = any(heart_type != "高血压"),
      # 是否存在广义心血管病
      has_any_cvd = TRUE
    )
  
  cat("=== 重新分类后的统计 ===\n")
  cat("包含广义心血管病（含高血压）的住院数:", nrow(cvd_summary_final), "\n")
}
df_diagnose_other1 <- df_diagnose %>% filter(!hadm_id %in% cvd_summary_final$hadm_id)
setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
write.xlsx(cvd_summary_final,"mimic4_sepsis_seqnum1_msm_脓毒症心脏病诊断.xlsx")
write.xlsx(df_diagnose_other1,"mimic4_sepsis_seqnum1_msm_脓毒症诊断_没有筛选出心脏病人群1.xlsx")


{
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  cat("========== 1. 宏观层面：单纯高血压 vs 结构性心脏病 ==========\n")
  macro_summary <- cvd_summary_final %>%
    count(has_major_heart_disease) %>%
    mutate(
      分类 = ifelse(has_major_heart_disease, "合并重大心血管疾病 (Major CVD)", "单纯高血压 (Hypertension Only)"),
      人数 = n,
      比例_百分比 = round(n / sum(n) * 100, 2)
    ) %>%
    select(分类, 人数, 比例_百分比) %>%
    arrange(desc(人数))
  
  print(as.data.frame(macro_summary))
  
  
  cat("\n========== 2. 疾病负担：单个患者合并的疾病数量 (Comorbidity Burden) ==========\n")
  burden_summary <- cvd_summary_final %>%
    mutate(
      # 通过计算 "+" 号的数量再加 1，得出该患者有几种心血管疾病
      疾病数量 = str_count(cvd_types, "\\+") + 1
    ) %>%
    count(疾病数量, name = "患者人数") %>%
    mutate(
      比例_百分比 = round(患者人数 / sum(患者人数) * 100, 2)
    )
  
  print(as.data.frame(burden_summary))
  
  
  cat("\n========== 3. 单病种频率：各种具体心血管疾病的发病率 ==========\n")
  # 这里使用了 separate_rows，它是非常强大的长宽表转换函数，能把 "高血压+心衰" 拆成两行分别统计
  individual_disease_summary <- cvd_summary_final %>%
    separate_rows(cvd_types, sep = "\\+") %>%
    count(cvd_types, name = "患病人数", sort = TRUE) %>%
    mutate(
      # 注意：这里的基数是总人数 (nrow)，这意味着一个人可以贡献多个疾病
      占所有CVD患者的比例_百分比 = round(患病人数 / nrow(cvd_summary_final) * 100, 2)
    )
  
  print(as.data.frame(individual_disease_summary))
  
  
  cat("\n========== 4. 疾病组合模式：最常见的疾病组合 Top 10 ==========\n")
  combo_summary <- cvd_summary_final %>%
    count(cvd_types, name = "患者人数", sort = TRUE) %>%
    mutate(
      比例_百分比 = round(患者人数 / sum(患者人数) * 100, 2)
    ) %>%
    head(10) # 只看前 10 种最常见的组合
  
  print(as.data.frame(combo_summary))
}
# ========== 1. 宏观层面：单纯高血压 vs 结构性心脏病 ==========
#   分类 人数 比例_百分比
# 1 合并重大心血管疾病 (Major CVD)  525        97.4
# 2 单纯高血压 (Hypertension Only)   14         2.6
# 
# ========== 2. 疾病负担：单个患者合并的疾病数量 (Comorbidity Burden) ==========
#   疾病数量 患者人数 比例_百分比
# 1        1       47        8.72
# 2        2      108       20.04
# 3        3      111       20.59
# 4        4      131       24.30
# 5        5       76       14.10
# 6        6       45        8.35
# 7        7       16        2.97
# 8        8        3        0.56
# 9        9        2        0.37
# 
# ========== 3. 单病种频率：各种具体心血管疾病的发病率 ==========
#   cvd_types 患病人数 占所有CVD患者的比例_百分比
# 1         高血压      397                      73.65
# 2           心衰      315                      58.44
# 3         冠心病      283                      52.50
# 4   其他心律失常      272                      50.46
# 5           心梗      212                      39.33
# 6  肺心病/血管病      119                      22.08
# 7         瓣膜病      103                      19.11
# 8           房颤       94                      17.44
# 9     心源性休克       59                      10.95
# 10        心肌病       54                      10.02
# 11      心包疾病       16                       2.97
# 
# ========== 4. 疾病组合模式：最常见的疾病组合 Top 10 ==========
#   cvd_types 患者人数 比例_百分比
# 1                             高血压+冠心病+心衰       20        3.71
# 2                                  高血压+冠心病       18        3.34
# 3                        高血压+冠心病+心梗+心衰       18        3.34
# 4                            高血压+其他心律失常       17        3.15
# 5                                   其他心律失常       15        2.78
# 6                                         高血压       14        2.60
# 7                     高血压+冠心病+其他心律失常       14        2.60
# 8           高血压+冠心病+其他心律失常+心梗+心衰       14        2.60
# 9                高血压+冠心病+其他心律失常+心衰       13        2.41
# 10 肺心病/血管病+高血压+冠心病+其他心律失常+心衰       11        2.04