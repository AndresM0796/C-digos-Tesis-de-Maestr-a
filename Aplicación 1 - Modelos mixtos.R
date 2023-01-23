#### Paquetes a cargar ####
library(MEMSS)
library(nlme)
library(MVA)
library(lattice)
library(mlmRev)
library(dplyr)
library(lme4)
library(ggplot2)
library(RLRsim) 
library(HLMdiag)
library(plotly)
library(ggpubr)
library(latex2exp)
library(stargazer)
library(xtable)

#### Base de datos ####
data(Dialyzer)
names(Dialyzer)

# Crear una copia de la base de datos para las figuras
Dialyzer2 <- Dialyzer %>% 
  dplyr::select(-c("index")) %>% 
  rename("Dializador" = "Subject", "Presion" = "pressure", "RUF" = "rate") 
Dialyzer2$QB = ifelse(Dialyzer2$QB == "200", "200 dL/min", "300 dL/min")


## Figura 4.1
Dialyzer2 %>%
  group_by(Dializador) %>% 
  ggplot(aes(Presion, RUF, group = Dializador, color = Dializador)) +
  geom_line() + geom_point() + facet_wrap(~QB) +
  xlab("Presión (dmHg)") + ylab("Razón de ultrafiltración (ml/hr)")

#### modelo en (4.1) ####
# el modelo se ajusta con maxima verosimilitud para realizar comparaciones 
# con el LRT
mod1 <- lme(rate ~ pressure + QB, random = ~ 1|Subject,
            data = Dialyzer, method = "ML") 

## estimacion del modelo con gls (ecuacion (4.2))
mod_gls <- gls(rate ~ pressure + QB, data = Dialyzer, 
               method = "ML", na.action = "na.omit")

## Test de componentes de varianza con RLRsim
exactRLRT(mod1)

## Comparacion con el LRT (mod_gls vs mod1)
anova(mod_gls, mod1)


#### modelos con terminos polinomiales #####

## modelo mixto ecuacion (4.3)
mod2 <- lme(rate ~ pressure + I(pressure^2)+ I(pressure^3) + QB,
            random = ~ 1|Subject,
            data = Dialyzer, method = "ML")

## modelo gls ecuacion (4.4)
mod_gls2 <- gls(rate ~ pressure + I(pressure^2)+ I(pressure^3) + QB,
                data = Dialyzer, method = "ML", na.action = "na.omit")


## test de razon de verosimilitud
anova(mod_gls, mod_gls2, mod2)



## termino cubico y pendiente aleatoria en pressure, ecuacion (4.5)
mod3 <- lme(rate ~ pressure + I(pressure^2)+ I(pressure^3) + QB,
            random = ~ pressure|Subject,
            data = Dialyzer, method = "ML")

anova(mod2, mod3)


## termino cuartico
mod4 <- lme(rate ~ pressure + I(pressure^2)+ I(pressure^3) + I(pressure^4) + QB,
            random = ~ pressure|Subject,
            data = Dialyzer, method = "ML")

anova(mod3, mod4) # no mejora

# Comparacion de los 3 modelos
anova(mod2, mod3, mod4)

# Comparacion del modelo mixto y el modelo gls con terminos polinomiales
anova(mod_gls2, mod3)


#### graficos de ajuste ####
# Funcion de panel para el gráfico de ajuste
pfun <- function(x,y){
  panel.xyplot(x,y[1:length(x)])
  panel.lines(x,y[1:length(x)+length(x)], lty = 1)
}


## Figura 4.2
Dialyzer$pred <- predict(mod3)
plot(xyplot(cbind(rate,pred) ~ pressure | Subject,
            data = Dialyzer, panel = pfun, 
            ylab = "Razón de ultrafiltración (ml/hr)"))


#### modelo definitivo ####
## en nlme
mod.nlme <- lme(rate ~ pressure + I(pressure^2)+ I(pressure^3) + QB,
                random = ~ pressure|Subject,
                data = Dialyzer, method = "REML")

## en lme4
rate.mod <- lmer(rate ~ pressure + I(pressure^2) + I(pressure^3)+ QB
                 +(pressure|Subject),data = Dialyzer)

# Las estimaciones de ambos modelos son iguales


#### Calculo de los residuos ####
res.data.lvl1 <- residuos_LMM_1(rate.mod)
res.data.lvl2 <- residuos_LMM_2(rate.mod)
r_2 <- data.frame(ranef(rate.mod)$Subject) # efectos aleatorios

res.data.lvl1$pressure <- Dialyzer$pressure
res.data.lvl1$pressure2 <- Dialyzer$pressure^2
res.data.lvl1$pressure3 <- Dialyzer$pressure^3

#### graficos ####
## Figura 4.3
res.data.lvl1 %>% 
  ggplot(aes(y_hat, r1.est)) + geom_point(cex=2) + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + theme_gray()+
  xlab(expression(hat(y))) + ylab("Residuos nivel 1 est.") 
res.data.lvl1 %>% 
  ggplot(aes(y_hat, r1.sem)) + geom_point(cex=2) + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + theme_gray()+
  xlab(expression(hat(y))) + ylab("Residuos semiestandarizados")

## Figura 4.4
res.data.lvl1 %>% 
  ggplot(aes(pressure, r1.sem)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + 
  ylab("Residuos semiestandarizados") + xlab("Presión (dmHg)")
res.data.lvl1 %>% 
  ggplot(aes(pressure2, r1.sem)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) +
  ylab("Residuos semiestandarizados") + xlab("Presión^2 (dmHg^2)")
res.data.lvl1 %>% 
  ggplot(aes(pressure, r1.sem)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + 
  ylab("Residuos semiestandarizados") + xlab("Presión^3 (dmHg^3)")

## Figura 4.5
qqnorm(res.data.lvl1$r1.sem, las=1, pch = 16, xlab = "Cuantiles teóricos",
       ylab = "Cuantiles de la muestra", main = "Q-Q plot Normal")
qqline(res.data.lvl1$r1.sem, col = "red", lty = 2)


## Figura 4.6
residdiag3.nlme(mod.nlme, limit=2, plotid=1)

## Figura 4.7
residdiag3.nlme(mod.nlme, limit=2, plotid=3)

## Figura 4.8
res.data.lvl2 %>% 
  dotplot_diag(x = sigma2, cutoff = "internal", name = "rvc") +
  ylab(TeX(r'(RVC $\sigma^2$)'))
res.data.lvl2 %>% 
  dotplot_diag(x = D11, cutoff = "internal", name = "rvc") +
  ylab(TeX(r'(RVC $\sigma_{0 0}$)'))
res.data.lvl2 %>% 
  dotplot_diag(x = D21, cutoff = "internal", name = "rvc") +
  ylab(TeX(r'(RVC $\sigma_{0 1}$)'))
res.data.lvl2 %>% 
  dotplot_diag(x = D22, cutoff = "internal", name = "rvc") +
  ylab(TeX(r'(RVC $\sigma_{1 1}$)'))

## Figura 4.9
res.data.lvl1 %>% 
  ggplot(aes(1:dim(res.data.lvl1)[1], r1.est)) + geom_point() + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos nivel 1 est.")+
  xlab("índice") +
  geom_text(aes(label=ifelse(abs(r1.est)>2,as.character(1:140),'')),hjust=0,vjust=0)

## Figura 4.10
residdiag3.nlme(mod.nlme, limit=2, plotid=4)

## Figura 4.11
res.data.lvl2 %>% 
  dotplot_diag(x = fixef, cutoff = "internal", name = "leverage") +
  ylab("Leverage efectos fijos")
res.data.lvl2 %>% 
  dotplot_diag(x = ranef.uc, cutoff = "internal", name = "leverage") +
  ylab("Leverage efectos aleatorios")

## Figura 4.12
res.data.lvl2 %>% 
  dotplot_diag(x = ci_2, cutoff = "internal", name = "cooks.distance") +
  ylab("Cook's distance") + xlab("school")
res.data.lvl1 %>% 
  dotplot_diag(x = ci_1, cutoff = "internal", name = "cooks.distance",
               xlim = c(0,0.15))+ylab("Cook's distance") + xlab("school")

## Figura 4.13
residdiag3.nlme(mod.nlme, limit=2, plotid=7)

## Figura 4.14
residdiag3.nlme(mod.nlme, limit=2, plotid=8)

## Figura 4.15
residdiag3.nlme(mod.nlme, limit=2, plotid=9)

### comparaciones en las estimaciones
## modelo original
summary(rate.mod)

## modelo sin los puntos 84,83 y 7
library(influence.ME)
mod.new <- exclude.influence(rate.mod, obs = c(7,83,84))
summary(mod.new)

## Figura 4.16
res.data.lvl2 %>% 
  dotplot_diag(x = covt_2, cutoff = "internal", name = "covtrace") +
  ylab("COVTRACE")
res.data.lvl1 %>% 
  dotplot_diag(x = covt_1, cutoff = "internal", name = "covtrace") +
  ylab("COVTRACE")