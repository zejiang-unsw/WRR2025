# Experiment on phases
rm(list=ls()) # remove all variables
graphics.off() # remove all figures

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

# load library ----
library(data.table)
library(ggplot2)
library(egg)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)

#display.brewer.all(colorblindFriendly = T, type=c("div","qual","seq","all")[2])
col_len <- 9
col_qual <- brewer.pal(n=col_len, name="Paired")
col_div <- brewer.pal(n=col_len, name="RdYlBu")
col_seq <- brewer.pal(n=col_len, name="Blues")

# global variables
mode <- switch(1, "DWT","MODWT")
wf <- switch(3, "haar","d4","d16")

flag.v <- switch(2, "","_v1")
# synthetic data ----
set.seed(2025-07-08)
ensemble <- 100

fd <- c("x","y")
n <- 5000
noise <- 0.1

ts.r <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = c(-2, -10, 0.2), n, s=noise)
ts.r <- WASP2.0::data.gen.Lorenz(sigma = 10, beta = 8/3, rho = 28, start = c(-13, -14, 47),
                                 time = seq(0, by=0.05, length.out = n), s=0)
data <- list(x=ts.r$z, dp=cbind(ts.r$x, ts.r$y))
plot.ts(cbind(data$x,data$dp), main="")

fig <- rgl::plot3d(ts.r$x, ts.r$y, ts.r$z,
            type = "l",
            size = 3,
            col = scales::alpha(col_qual[2], 0.6),
            box = TRUE,
            axes = FALSE,
            xlab = "",
            ylab = "",
            zlab = "")

fig %>% print()
rgl.postscript("plot.pdf", fmt="pdf")

data.list <- NULL
for(r in 1:ensemble){
  if(flag.v=="_v1"){
    z0 <- c(runif(1,-3,-1), runif(1,-11,-9), runif(1,0,0.5))
  } else if(flag.v=="") {
    z0 <- c(-2, -10, 0.2)
  }
  ts.r <- WASP2.0::data.gen.Rossler(a = 0.2, b = 0.2, w = 5.7, start = z0,
                                    n, s=noise)
  
  data.list[[r]] <- list(x=ts.r$z, dp=cbind(ts.r$x, ts.r$y))
  
}
sapply(data.list, function(ls) is.na(ls$x)) %>% sum()
data1 <- data.list[[sample(1:100,1)]]
plot.ts(cbind(data1$x,data1$dp), main="")

# evaluation ----
k.folds <- 2
folds <- cut(seq(1,n),breaks=k.folds,labels=FALSE)

file_dat <- paste0("../result/Rossler_",wf,"_n",n,"_sd",noise,
                   "_r",ensemble,flag.v,".Rdata")

if(!file.exists(file_dat)) {
  # Register parallel backend
  n_cores <- parallel::detectCores() - 1
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  clusterExport(cl, varlist=c("folds","data.list","wf"), envir=environment())
  clusterEvalQ(cl, {
    library(WASP2.0)
    library(FNN)
    library(dplyr)
  })
  
  # Parallel loop
  pred_df_list <- foreach(j = iter(1:ensemble,by=11), .combine = rbind) %:%
    foreach(k_fold = 1:k.folds, .combine = rbind) %dopar% {
      
  # k_fold = 1
  ind_cal <- which(folds!=k_fold)
  ind_val <- which(folds==k_fold)

  # WASP
  dwt <- dwt.wasp(data.list[[j]], wf, index_c = ind_cal, phase=FALSE,verbose = F)
  
  if(flag.v==""){
    dwt1 <- dwt.wasp(data.list[[j]], wf, index_c = ind_cal, phase=TRUE, 
                   k=5, verbose = F)
  } else {
    dwt1 <- WASP2.0::dwt.wasp(data.list[[j]], wf, index_c = ind_cal, phase=TRUE,
                              k_auto = TRUE, power_thresh = 0.8, k_max = n/2-1, #max_iter=3, 
                              #momentum = 0.3,
                              verbose = F)
  }
  
  pred_df <- NULL
  if(TRUE){
    i_dp <- 1:2
    # linear model 
    fit1 <- stats::lm(dwt$x[ind_cal] ~., data= dwt$dp[ind_cal,i_dp]%>% data.frame())
    fit2 <- stats::lm(dwt$x[ind_cal] ~., data= dwt$dp.n[ind_cal,i_dp]%>% data.frame())
    fit3 <- stats::lm(dwt1$x[ind_cal] ~., data= dwt1$dp.n[ind_cal,i_dp]%>% data.frame())
    
    pred_cal_lm <- data.frame(method="lm", group="cal", obs=dwt$x[ind_cal], 
                              mod1=fit1$fitted.values, mod2=fit2$fitted.values, 
                              mod3=fit3$fitted.values)
    
    mod1 <- stats::predict(fit1, newdata = dwt$dp[ind_val,i_dp] %>% data.frame())
    mod2 <- stats::predict(fit2, newdata = dwt$dp.n[ind_val,i_dp]%>% data.frame())
    mod3 <- stats::predict(fit3, newdata = dwt1$dp.n[ind_val,i_dp]%>% data.frame())
    
    #plot(dwt$x[ind_val], mod2)
    
    pred_val_lm <- data.frame(method="lm", group="val", obs=dwt$x[ind_val], 
                              mod1=mod1, mod2=mod2, mod3=mod3)
    # knn model
    k = 5
    #fit1 <- knn(dwt$x[ind_cal], z=dwt$dp[ind_cal,i_dp], zout=dwt$dp[ind_cal,i_dp], k=k)
    #fit2 <- knn(dwt$x[ind_cal], z=dwt$dp.n[ind_cal,i_dp], zout=dwt$dp.n[ind_cal,i_dp], k=k)
    
    fit1 <- FNN::knn.reg(train=dwt$dp[ind_cal,i_dp], y=dwt$x[ind_cal], k=k)
    fit2 <- FNN::knn.reg(train=dwt$dp.n[ind_cal,i_dp], y=dwt$x[ind_cal], k=k)
    fit3 <- FNN::knn.reg(train=dwt1$dp.n[ind_cal,i_dp], y=dwt1$x[ind_cal], k=k)
    
    pred_cal_knn <- data.frame(method="knn", group="cal", obs=dwt$x[ind_cal], 
                               mod1=fit1$pred, mod2=fit2$pred, mod3=fit3$pred)
    
    # fit1 <- knn(dwt$x[ind_cal], z=dwt$dp[ind_cal,i_dp], zout=dwt$dp[ind_val,i_dp], k=k)
    # fit2 <- knn(dwt$x[ind_cal], z=dwt$dp.n[ind_cal,i_dp], zout=dwt$dp.n[ind_val,i_dp], k=k)
    
    fit1 <- FNN::knn.reg(train=dwt$dp[ind_cal,i_dp] %>% as.matrix(), test=dwt$dp[ind_val,i_dp]%>% as.matrix(), 
                         y=dwt$x[ind_cal], k=k)
    fit2 <- FNN::knn.reg(train=dwt$dp.n[ind_cal,i_dp]%>% as.matrix(), test=dwt$dp.n[ind_val,i_dp]%>% as.matrix(), 
                         y=dwt$x[ind_cal], k=k)
    fit3 <- FNN::knn.reg(train=dwt1$dp.n[ind_cal,i_dp]%>% as.matrix(), test=dwt1$dp.n[ind_val,i_dp]%>% as.matrix(), 
                         y=dwt1$x[ind_cal], k=k)
    
    pred_val_knn <- data.frame(method="knn", group="val", obs=dwt$x[ind_val], 
                               mod1=fit1$pred, mod2=fit2$pred, mod3=fit3$pred)
    
    pred_df <- rbind(pred_df, data.frame(r=j, fold=k_fold, rbind(pred_cal_lm, pred_val_lm, 
                                                                                pred_cal_knn, pred_val_knn)))
    
  }
  
  return(pred_df)
}

  stopCluster(cl)
  
  ## save ----
  save(pred_df_list, file = file_dat)
} else {
  load(file_dat)
}

summary(pred_df_list)
# metrics ----
pred_dt <- pred_df_list %>% gather(model, pred, which(colnames(pred_df_list) %in% paste0("mod",1:3))) %>% 
  mutate_if(is.character, as.factor) %>% 
  data.table()
summary(pred_dt)

metric_df <- pred_dt[,.(RMSE=Metrics::rmse(obs,pred),
                        R = cor(obs,pred)), by=c("method", "group", "model","r")]
metric_df$group <- factor(metric_df$group, labels = c("Calibration","Validation"))
metric_df$model <- factor(metric_df$model, labels=c("x","x'","x'(φ)"))
summary(metric_df)

# plot ----
## RMSE & R ----
p_RMSE <- ggplot(metric_df, 
                 aes(x=model, y=RMSE,fill=model)) + 
  #geom_bar(position="dodge", stat="identity") + 
  geom_boxplot(outlier.colour="red") + 
  facet_grid(method~group) + 
  
  #scale_y_continuous(limits = c(0,5)) + 
  scale_fill_manual(values=c("black","red","blue")) + 
  
  labs(title="RMSE") + 
  #theme_article() +
  theme(
    legend.position = "none"
  )

p_RMSE #%>% print()

p_R <- ggplot(metric_df, 
              aes(x=model, y=R,fill=model)) + 
  #geom_bar(position="dodge", stat="identity") + 
  geom_boxplot(outlier.colour="red") + 
  facet_grid(method~group) + 
  
  #scale_y_continuous(limits = c(-1,1)) + 
  scale_fill_manual(values=c("black","red","blue")) + 
  
  labs(title="R") + 
  #theme_article() +
  theme(
    legend.position = "none"
  )

p_R #%>% print()

fig <- cowplot::plot_grid(p_RMSE, p_R, ncol=2)  
print(fig)

# filen <- paste0("Figure_Rossler_",wf,"_n",n,"_sd",noise,"_r",ensemble,
#                 flag.v,".png")
# png(filen, height=10, width=12, units="cm", bg = "transparent", res=300)
# fig %>% print()
# dev.off()

