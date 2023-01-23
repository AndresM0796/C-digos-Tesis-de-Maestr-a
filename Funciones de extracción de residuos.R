### Estas funciones deben ser ejecutadas para que las aplicaciones funcionen 
# correctamente.


#### Funcion para extraer residuos de un MLG ####
residuos_glm <- function(mod){
  require(statmod)
  r_p <- resid(mod, type = "pearson") # Pearson
  r_d <- resid(mod) # desvio
  r_q <- qresid(mod) # cuantil
  
  h_i <- hatvalues(mod) # leverages
  phi.P <- summary(mod)$dispersion # phi estimado con Pearson
  
  r_sp <- rstandard(mod, type = "pearson") # Pearson est.
  r_sd <- rstandard(mod) # desvio est.
  r_sq <- r_q/sqrt(1-h_i) # cuantil est.
  
  r_G <- rstudent(mod) # Williams
  
  mu_hat <- fitted(mod) # mu estimado
  eta_hat <- predict(mod, type = "link") # eta estimado
  z <- resid(mod,type = "working") + eta_hat # respuesta linealizada
  c_i <- cooks.distance(mod)
  
  res.data <- data.frame(r_p,r_d,r_q,h_i,r_sp,r_sd,r_sq,r_G,
                         mu_hat, eta_hat, z, c_i)
  return(res.data)
}

#### Funcion para extraer residuos en el nivel 1 de un MLM ####
residuos_LMM_1 <- function(mod.lme4){
  # mod es un modelo de lme4
  require(HLMdiag)
  require(lme4)
  require(nlme)
  
  res_fm1 <- hlm_resid(mod.lme4, level = 1, type = "LS",standardize = F)
  res_fm1s <- hlm_resid(mod.lme4, level = 1, type = "LS",standardize = T)
  res_fmsemi <- hlm_resid(mod.lme4, level = 1, type = "LS", standardize = "semi")
  
  r1 <- residuals(mod.lme4)
  r1.est <- res_fm1s$.std.resid
  r1.sem <- res_fmsemi$.semi.ls.resid
  r_m <- res_fm1$.mar.resid
  r_ms <- res_fm1s$.chol.mar.resid
  y_hat <- predict(mod.lme4)
  ci_1 <- cooks.distance(mod.lme4, level = 1)
  mdf_1 <- mdffits(mod.lme4, level = 1)
  covr_1 <- covratio(mod.lme4, level = 1)
  covt_1 <- covtrace(mod.lme4, level = 1)
  rvc_1 <- rvc(mod.lme4, level = 1)
  lev_1 <- leverage(mod.lme4, level = 1)
  
  res.data.lvl1 <- data.frame(r1, r1.est, r1.sem, r_m, r_ms, y_hat, ci_1,
                              mdf_1, covr_1, covt_1, rvc_1, lev_1)
  return(res.data.lvl1)
}


#### Funcion para obtener medidas de influencia en el nivel 2 de un MLM ####
residuos_LMM_2 <- function(mod.lme4){
  require(HLMdiag)
  require(lme4)
  require(nlme)
  grupo <- names(summary(mod.lme4)$ngrps)
  ci_2 <- cooks.distance(mod.lme4, level = grupo)
  mdf_2 <- mdffits(mod.lme4, level = grupo)
  covr_2 <- covratio(mod.lme4, level = grupo)
  covt_2 <- covtrace(mod.lme4, level = grupo)
  rvc_2 <- rvc(mod.lme4, level = grupo)
  lev_2 <- leverage(mod.lme4, level = grupo)
  
  res.data.lvl2 <- data.frame(ci_2, mdf_2, covr_2, covt_2, rvc_2, lev_2)
  return(res.data.lvl2)
}
