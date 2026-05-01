# Nd_PBPK_mrgsolve.R  修正版
# 小鼠钕（Nd）PBPK模型 - deSolve 实现

library(deSolve)

# ======================== 参数定义 ========================
BW      <- 0.022         # kg
Qc      <- 5.4           # L/h/kg
Hct     <- 0.45
Vplasma <- 0.07 * BW * (1 - Hct)

# 组织容积 (L)
V_lung   <- 0.006 * BW
V_brain  <- 0.019 * BW
V_heart  <- 0.005 * BW
V_spleen <- 0.005 * BW
V_liver  <- 0.055 * BW
V_kidney <- 0.015 * BW
V_rich   <- 0.10  * BW
V_slow   <- 0.40  * BW
V_fat    <- 0.10  * BW
V_bone   <- 0.07  * BW

# 组织血流 (L/h)
Q_total  <- Qc * BW
Q_lung   <- 0.02 * Q_total
Q_brain  <- 0.02 * Q_total
Q_heart  <- 0.04 * Q_total
Q_spleen <- 0.02 * Q_total
Q_liver  <- 0.25 * Q_total
Q_kidney <- 0.19 * Q_total
Q_rich   <- 0.18 * Q_total
Q_slow   <- 0.18 * Q_total
Q_fat    <- 0.07 * Q_total
Q_bone   <- 0.05 * Q_total

# 化合物参数
Ka       <- 0.8
F_bio    <- 0.15
CL_renal <- 0.18 * 60 / 1000   # L/h/kg
CL       <- CL_renal * BW       # L/h

# 分配系数
Kp_lung   <- 1.5
Kp_brain  <- 0.5
Kp_heart  <- 1.2
Kp_spleen <- 1.8
Kp_liver  <- 2.5
Kp_kidney <- 3.0
Kp_rich   <- 1.0
Kp_slow   <- 0.8
Kp_fat    <- 0.3
Kp_bone   <- 0.9

# ======================== 微分方程 ========================
pbpk_ode <- function(time, state, parameters) {
  A_abs    <- state[1]
  A_plasma <- state[2]
  A_lung   <- state[3]
  A_brain  <- state[4]
  A_heart  <- state[5]
  A_spleen <- state[6]
  A_liver  <- state[7]
  A_kidney <- state[8]
  A_rich   <- state[9]
  A_slow   <- state[10]
  A_fat    <- state[11]
  A_bone   <- state[12]
  
  C_plasma <- A_plasma / Vplasma
  
  # 吸收室
  dA_abs <- -Ka * A_abs
  
  # 血浆室 (所有回流血流量 + 吸收 - 流出 - 肾清除)
  dA_plasma <- Ka * A_abs +
    Q_lung   * A_lung   / (Kp_lung   * V_lung)   +
    Q_brain  * A_brain  / (Kp_brain  * V_brain)  +
    Q_heart  * A_heart  / (Kp_heart  * V_heart)  +
    Q_spleen * A_spleen / (Kp_spleen * V_spleen)  +
    Q_liver  * A_liver  / (Kp_liver  * V_liver)  +
    Q_kidney * A_kidney / (Kp_kidney * V_kidney)  +
    Q_rich   * A_rich   / (Kp_rich   * V_rich)   +
    Q_slow   * A_slow   / (Kp_slow   * V_slow)   +
    Q_fat    * A_fat    / (Kp_fat    * V_fat)    +
    Q_bone   * A_bone   / (Kp_bone   * V_bone)   -
    Q_total * C_plasma -
    CL * C_plasma
  
  # 各组织 (流出 = 血流 * (动脉浓度 - 静脉浓度/Kp))
  dA_lung   <- Q_total * (C_plasma - A_lung   / (Kp_lung   * V_lung))
  dA_brain  <- Q_brain * (C_plasma - A_brain  / (Kp_brain  * V_brain))
  dA_heart  <- Q_heart * (C_plasma - A_heart  / (Kp_heart  * V_heart))
  dA_spleen <- Q_spleen* (C_plasma - A_spleen / (Kp_spleen * V_spleen))
  dA_liver  <- Q_liver * (C_plasma - A_liver  / (Kp_liver  * V_liver))
  dA_kidney <- Q_kidney* (C_plasma - A_kidney / (Kp_kidney * V_kidney))
  dA_rich   <- Q_rich  * (C_plasma - A_rich   / (Kp_rich   * V_rich))
  dA_slow   <- Q_slow  * (C_plasma - A_slow   / (Kp_slow   * V_slow))
  dA_fat    <- Q_fat   * (C_plasma - A_fat    / (Kp_fat    * V_fat))
  dA_bone   <- Q_bone  * (C_plasma - A_bone   / (Kp_bone   * V_bone))
  
  return(list(c(dA_abs, dA_plasma, dA_lung, dA_brain, dA_heart,
                dA_spleen, dA_liver, dA_kidney, dA_rich, dA_slow,
                dA_fat, dA_bone)))
}

# ======================== 模拟函数 ========================
simulate_dose <- function(dose_mg_kg = 344, times = seq(0, 24, by = 0.1)) {
  A_abs0 <- dose_mg_kg * BW * F_bio
  state <- c(A_abs = A_abs0, A_plasma = 0,
             A_lung = 0, A_brain = 0, A_heart = 0,
             A_spleen = 0, A_liver = 0, A_kidney = 0,
             A_rich = 0, A_slow = 0, A_fat = 0, A_bone = 0)
  
  out <- ode(y = state, times = times, func = pbpk_ode, parms = NULL)
  out <- as.data.frame(out)
  
  # 添加组织浓度
  out$C_plasma <- out$A_plasma / Vplasma
  out$C_lung   <- out$A_lung   / V_lung
  out$C_brain  <- out$A_brain  / V_brain
  out$C_heart  <- out$A_heart  / V_heart
  out$C_spleen <- out$A_spleen / V_spleen
  out$C_liver  <- out$A_liver  / V_liver
  out$C_kidney <- out$A_kidney / V_kidney
  out$C_rich   <- out$A_rich   / V_rich
  out$C_slow   <- out$A_slow   / V_slow
  out$C_fat    <- out$A_fat    / V_fat
  out$C_bone   <- out$A_bone   / V_bone
  
  return(out)
}

# ======================== 运行三个剂量 ========================
doses <- c(344, 688, 1375)
sim   <- lapply(setNames(doses, doses), simulate_dose)

# ======================== 绘图 ========================
par(mfrow = c(2,2))

# 血浆
plot(sim[["344"]]$time, sim[["344"]]$C_plasma, type="l",
     xlab="Time (h)", ylab="mg/L", main="Plasma")
lines(sim[["688"]]$time, sim[["688"]]$C_plasma, col=2)
lines(sim[["1375"]]$time, sim[["1375"]]$C_plasma, col=3)
legend("topright", legend=c("344","688","1375"), col=1:3, lty=1)

# 脑
plot(sim[["344"]]$time, sim[["344"]]$C_brain, type="l",
     xlab="Time (h)", ylab="mg/L", main="Brain")
lines(sim[["688"]]$time, sim[["688"]]$C_brain, col=2)
lines(sim[["1375"]]$time, sim[["1375"]]$C_brain, col=3)

# 肝
plot(sim[["344"]]$time, sim[["344"]]$C_liver, type="l",
     xlab="Time (h)", ylab="mg/L", main="Liver")
lines(sim[["688"]]$time, sim[["688"]]$C_liver, col=2)
lines(sim[["1375"]]$time, sim[["1375"]]$C_liver, col=3)

# 肾
plot(sim[["344"]]$time, sim[["344"]]$C_kidney, type="l",
     xlab="Time (h)", ylab="mg/L", main="Kidney")
lines(sim[["688"]]$time, sim[["688"]]$C_kidney, col=2)
lines(sim[["1375"]]$time, sim[["1375"]]$C_kidney, col=3)

# 导出CSV（可选）
# pred344 <- sim[["344"]]; write.csv(pred344, "prediction_344.csv", row.names=FALSE)