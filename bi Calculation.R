# 读取输入数据
input_data <- read.csv("input.csv", header = TRUE)
Dx <- as.numeric(input_data[input_data[, 1] == "Dx", -1])
Dy <- as.numeric(input_data[input_data[, 1] == "Dy", -1])

A_data <- read.csv("A1 to A11.csv", header = TRUE)
# 修正矩阵维度：A_means应为11个源 x 11个元素
A_means <- as.matrix(A_data[, 2:12])  # 11行(源) x 11列(元素)
A_sds <- as.matrix(A_data[, 13:23])

# 修正后的贝叶斯模型定义
model_code <- "
model {
  # 先验分布：根据原计算结果调整先验增强变异性
  for(i in 1:1) {
    b[i] ~ dbeta(1, 99) # 抑制其他系数
  }
    b[2] ~ dbeta(1, 99)    #
  for(i in 3:8) {
    b[i] ~ dbeta(1, 99) # 抑制其他系数
  }
    b[9] ~ dbeta(99, 1)    # 
  for(i in 10:11) {
    b[i] ~ dbeta(1, 99) # 抑制其他系数
  }

  
  # 计算预测C'（修正矩阵运算）
  for(j in 1:11) {
    A_contribution[j] <- inprod(A_means[j,], b)  # 使用内积计算每个源的贡献
  }
  pred_C <- Dx - sum(A_contribution)
  C_prime <- pred_C * 100
  
  # 非负约束正则化系数（惩罚系数1e-1可以调节）
  penalty <- sum(ifelse(pred_C < 0, -pred_C * 1e0, 0))
  
  # Bray-Curtis相似度
  for(k in 1:11){
    C_nonneg[k] <- max(C_prime[k], 0)
  }
  bc_distance <- sum(abs(C_nonneg - Dy)) / (sum(C_nonneg) + sum(Dy))
  
  # 综合目标函数（调整似然函数）
  total_score <- bc_distance + penalty
  dummy ~ dnorm(total_score, 0.01)  # 降低精度参数避免溢出
}
"

# 运行贝叶斯模型
library(R2jags)
jags_data <- list(
  Dx = Dx,
  Dy = Dy,
  A_means = A_means,
  dummy = 0
)

params <- c("b", "bc_distance")

# 增加适配器参数防止采样卡住
fit <- jags(
  data = jags_data,
  parameters.to.save = params,
  model.file = textConnection(model_code),
  n.chains = 3,
  n.iter = 10000,  # 迭代次数
  n.burnin = 3000,
  n.thin = 2,
  jags.seed = 123,
  progress.bar = "text"
)

# 提取后验分布
b_samples <- fit$BUGSoutput$sims.list$b
b_mean <- apply(b_samples, 2, mean)
b_sd <- apply(b_samples, 2, sd)

# 验证结果（修正矩阵乘法）
final_C <- Dx - A_means %*% b_mean  # 直接矩阵乘法
cat("优化后的系数分布:\n")
print(data.frame(
  Source = paste0("A",1:11),
  Mean = round(b_mean,4),
  SD = round(b_sd,4)
))

cat("\n最终C的负数个数:", sum(final_C < 0))
cat("\nBray-Curtis距离:", round(mean(fit$BUGSoutput$sims.list$bc_distance),4))