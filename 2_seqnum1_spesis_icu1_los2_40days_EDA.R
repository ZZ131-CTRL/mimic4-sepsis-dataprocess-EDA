###################################################
# 脓毒症患者第一次进入icu病程CT-HMM数据EDA分析
# DATE:20260325
# 注意这里数据是筛选了 seqnum1为脓毒症且病程2-40天患者数据处理
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
# session1 基础设置和function定义 #---------------------------------------------------------------- 
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


# session2 数据检视 #----------------------------------------------------------------
# 1. 提取首诊 ICU 期间的 labevents 记录
{
  setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
  df_final_paths <- readRDS("lab_sepsis_seqnum1_finalpath.rds")
  all_data_rel <- readRDS("lab_sepsis_seqnum1_filtered.rds") 

  df_final_paths <- df_final_paths %>% rename(first_intime = intime, first_outtime = outtime)
  colnames(df_final_paths)
  colnames(all_data_rel)
  df_first_icu_course <- all_data_rel %>%
    # 先确保实验室数据的 ID 和时间格式正确
    mutate(across(c(subject_id, hadm_id), as.numeric) ,charttime = as.POSIXct(charttime)
    ) %>%
    
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
  
  setwd(Output_derived_data)
  df_diagnose <- read.xlsx("mimic4_sepsis_heart_诊断.xlsx") %>% filter(hadm_id %in% df_first_icu_course$hadm_id)
  length(unique(df_diagnose$hadm_id))
  setwd("D:\\Lab_project\\2026work\\sepsis\\DATA\\sepsis\\derived_data\\data_filtered")
  write.xlsx(df_diagnose,"mimic4_sepsis_seqnum1_msm_诊断.xlsx")
  
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
unique(cov1$race)
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

# 4. 床边监护指标数据统计 #-----------------------------------------------------------------
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

# 5. 床边监护指标统计表 #----------------------------------------------------
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

# 6. 实验室指标ct图 #-------------------------------------------------------------
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
base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_PreShock"
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
      # filter(charttime_rel_icu <= 13) %>%
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

# 休克发生前两组患者的实验室指标对比13days #------------------------------------------
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
# 非休克患者死亡存活组的实验室指标对比 #------------------------------------------
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
  base_path <- "D:\\Lab_project\\2026work\\sepsis\\PLOT\\sepsis\\CT_curve\\ICU_Biomarker_Plots_noShock_survival_death"
  
  for (group_name in names(label_groups)) {
    output_dir <- file.path(base_path, group_name)
    if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
    
    target_list <- label_groups[[group_name]]
    
    for (target_biomarker in target_list) {
      message("正在处理: ", target_biomarker)
      
      # --- 2. 数据准备 ---
      plot_df <- df_lab_log %>%
        filter(label == target_biomarker & !is.na(valuenum)) %>%
        filter(shock == 0) %>%                      # 仅保留no休克患者  
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
        model_linear <- lmer(valuenum ~ Group * charttime_rel_icu + (1 | subject_id), data = plot_df)
        coef_summary <- summary(model_linear)$coefficients
        # 提取 GroupDeceased 的 P 值
        p_linear_val <- format_p(coef_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # 对数混合模型
      try({
        model_log <- lmer(log10(valuenum + 0.001) ~ Group * charttime_rel_icu + (1 | subject_id), data = plot_df)
        coef_log_summary <- summary(model_log)$coefficients
        p_log_val <- format_p(coef_log_summary["GroupDeceased", "Pr(>|t|)"])
      }, silent = TRUE)
      
      # --- 4. 绘图 ---
      # 存活组 Survivor(3) 用绿色，死亡组 Deceased(4) 用红色
      color_values <- c("Deceased" = "#E41A1C", "Survivor" = "#377EB8")
      
      p_linear_plot <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Time (days)", y = target_biomarker, 
             title = paste0("Linear (Group P: ", p_linear_val, ")")) +
        theme_cor
      
      p_log_plot <- ggplot(plot_df, aes(x = charttime_rel_icu, y = valuenum, color = Group, fill = Group)) +
        geom_point(alpha = 0.2, size = 0.8, stroke = 0) + 
        geom_smooth(method = "loess", size = 1, alpha = 0.2) + 
        scale_y_log10() +
        scale_color_manual(values = color_values) +
        scale_fill_manual(values = color_values) +
        labs(x = "Time (days)", y = paste0("Log10 ", target_biomarker),
             title = paste0("Log-scale (Group P: ", p_log_val, ")")) +
        theme_cor
      
      # 合并保存
      combined_plot <- (p_linear_plot | p_log_plot) + 
        plot_layout(guides = 'collect') & 
        theme(legend.position = "bottom")
      
      safe_name <- gsub("[^[:alnum:]]", "_", target_biomarker)
      ggsave(file.path(output_dir, paste0(safe_name, "_noShock_Survivor_vs_Death.png")), 
             combined_plot, width = 12, height = 6, dpi = 300)
    }
  }
}

# 7. 床边监护ct图 #-------------------------------------------------------------
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
        geom_smooth(color = "red", method = "loess", size = 0.9, formula = y ~ x) +
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


# 6.2 实验室指标基线值_人口学信息_相关性分析 #-------------------------------
# 数据处理 #---------------------------------------------------
{
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
        geom_smooth(color="red",method = "loess",size=0.9)+
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

# 对race_group四个分组的相关性分析 ------------------------------------------------
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

# 输出有显著性的表----------------------------------------------------------------
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

# 输出有显著性的图 (连续变量 Pearson / 性别 t.test)------------------------------
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
            geom_smooth(color = "#FB8072", method = "loess", size = 0.9) +
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

