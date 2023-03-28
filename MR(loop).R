# Gwas在线分析数据
# 设置工作目录
setwd("/Users/charles/Desktop/MR")

library(TwoSampleMR)
# --ieu在线 【一键精灵】
# 【本地自己导入，环境保持一致】
# 暴露
exp_dat_id <- c("ieu-a-2") #体重指数GwasID

# 结局 {这边可以填入多个结局}
outcomeids <- c("ieu-a-7") #冠心病GwasID
out_comes <- c("Coronary heart disease") #冠心病

#
dir.create(path = "Charles_result") #创建文件夹
startALL=Sys.time() #记录开始时间

# 循环，可遍历多种疾病
for (i in 1:length(outcomeids)) {   #length(outcomeids),只有一个疾病，所有长度为1
  outcomeid <- outcomeids[i]  # 冠心病id
  out_come <- out_comes[i]  # 冠心病
  # ----{1.1}--exposure_data-------读取暴露数据
  exp_data <- extract_instruments(
    outcomes = exp_dat_id,
    clump =TRUE  #去除连锁不平衡【独立性设置】，可单独设置clump_data
    )
  # ----{1.2}--weak instruments------
  # xxxxxx
  # ----{2}--outcome_data-------读取结局数据
  out_data <- extract_outcome_data(
    snps = exp_data$SNP,  #取交集设置
    outcomes = outcomeid
    ) 
  # ----{3}--harmonise_data-----协同（保存，MR操作基础）1、等位基因方向协同2、剔除不能判断方向的snp 
  dat <- TwoSampleMR::harmonise_data(
      exposure_dat = exp_data,
      outcome_dat = out_data 
      )
  # ----{4}--MR_res-----MR分析
  res=TwoSampleMR::mr(dat)
  # 
  # check the least SNP counts
  print(paste0(out_come," SNP_num:",res$nsnp[1]))
  # --main_result--
  results <- TwoSampleMR::generate_odds_ratios(res) #results比res更详细（产生OR值）
  # OR
  results$estimate <- paste0(
    format(round(results$or, 2), nsmall = 2), " (", 
    format(round(results$or_lci95, 2), nsmall = 2), "-",
    format(round(results$or_uci95, 2), nsmall = 2), ")"
    )
  resdata <- dat
  # Assumption 1 and 3
  names(resdata)
  Assumption13 <- subset(
    resdata,
    mr_keep==TRUE,
    select = c("SNP","pval.exposure",
          "pval.outcome", # "F_statistic", 统计强度计算
          "mr_keep")
          )
  # -----{5}--Sensitive_analysis------敏感性分析
  res_hete <- TwoSampleMR::mr_heterogeneity(dat) #异质性检测
  res_plei <- TwoSampleMR::mr_pleiotropy_test(dat) #多效性检测
  res_leaveone <- mr_leaveoneout(dat)  # 
  # 
  res_presso <- TwoSampleMR::run_mr_presso(dat,
                                   NbDistribution = 100) #迭代次数，设置成10000，发现离群值
  # [["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
  res_presso[[1]][[2]][[1]][["Pvalue"]]
  sink(paste0("Charles_result/",out_come,"_PRESSO.txt"),
       append=FALSE,split = FALSE) 
  print(res_presso)
  sink()
  print(res_presso)
  # ------Merge_results--------
  # Export
  openxlsx::write.xlsx(x = list(
    "main"=results,
    "Assumption13"=Assumption13,
    "pleiotropy"=res_plei,
    "heterogeneity"=res_hete,
    "leaveone"=res_leaveone
    ),
    overwrite = TRUE,
    paste0("Charles_result/",out_come,".xlsx")
    )
  
  # -----{5}--Sensitive_analysis------四张图
    p1 <- mr_scatter_plot(res, dat)  #散点图可视化
    p1[[1]]
    pdf(paste0("Charles_result/",out_come,"_scatter.pdf"))
    print(p1[[1]])
    dev.off()
    #
    res_single <- mr_singlesnp(dat)
    p2 <- mr_forest_plot(res_single) #森林图
    pdf(paste0("Charles_result/",out_come,"_forest.pdf"))
    print(p2[[1]])
    dev.off()
    # 
    p3 <- mr_funnel_plot(res_single) #异质性，漏斗图
    pdf(paste0("Charles_result/",out_come,"_funnel.pdf"))
    print(p3[[1]])
    dev.off()
    # 
    res_loo <- mr_leaveoneout(dat) #留一交叉法
    pdf(paste0("Charles_result/",out_come,"_leave_one_out.pdf"))
    print(mr_leaveoneout_plot(res_loo))
    dev.off()
}
endALL=Sys.time();print(endALL-startALL) #记录时间




