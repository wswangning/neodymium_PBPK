# ============================================================
# sensitivity_final.R - 小鼠 Nd PBPK 全局敏感性分析
# ============================================================

library(sensitivity)

# ---- 1. 解析模型函数 ----
model_auc <- function(Ka, F, CL_renal, Kp_brain, Kp_liver) {
  BW   <- 0.022
  Dose <- 344 * BW
  CL_animal <- CL_renal * BW
  auc_plasma <- (F * Dose) / CL_animal
  auc_brain  <- Kp_brain * auc_plasma
  if (!is.finite(auc_brain) || auc_brain <= 0) return(NA)
  return(auc_brain)
}

# ---- 2. 包装函数（参数为命名向量） ----
wrapper <- function(X) {
  apply(X, 1, function(r) {
    # r 为行向量，按列顺序：Ka, F, CL_renal, Kp_brain, Kp_liver
    model_auc(Ka=r[1], F=r[2], CL_renal=r[3], Kp_brain=r[4], Kp_liver=r[5])
  })
}

# ---- 3. 生成样本 ----
set.seed(123)
n <- 1000

# 参数范围（基于您的基准值 ±50%）
X1 <- data.frame(
  Ka       = runif(n, 0.4, 1.2),
  F        = runif(n, 0.1, 0.2),
  CL_renal = runif(n, 0.005, 0.02),
  Kp_brain = runif(n, 0.3, 0.7),
  Kp_liver = runif(n, 1.5, 3.5)
)
X2 <- data.frame(
  Ka       = runif(n, 0.4, 1.2),
  F        = runif(n, 0.1, 0.2),
  CL_renal = runif(n, 0.005, 0.02),
  Kp_brain = runif(n, 0.3, 0.7),
  Kp_liver = runif(n, 1.5, 3.5)
)

# ---- 4. 计算输出 ----
cat("计算输出值...\n")
y1 <- wrapper(X1)
y2 <- wrapper(X2)

# ---- 5. 数据清洗 ----
valid <- is.finite(y1) & is.finite(y2)
cat(sprintf("有效样本: %d / %d\n", sum(valid), n))

if (sum(valid) < 100) stop("有效样本不足，请检查参数范围")

X1 <- X1[valid, ]
X2 <- X2[valid, ]
y1 <- y1[valid]
y2 <- y2[valid]

# ---- 6. Sobol 分析 ----
cat("正在执行 Sobol-Jansen 分析...\n")
sob <- soboljansen(model = wrapper, X1 = X1, X2 = X2, nboot = 100)
print(sob)

# ---- 7. 保存结果 ----
res <- sob$S
res <- cbind(Parameter = rownames(res), res)
write.csv(res, "Table_S6_Sensitivity_Analysis.csv", row.names = FALSE)

cat("✅ 完成！结果已保存。\n")