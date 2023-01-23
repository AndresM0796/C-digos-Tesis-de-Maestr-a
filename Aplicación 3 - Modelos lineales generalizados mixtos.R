#### Paquetes a cargar ####
library(mlmRev)
library(faraway)
library(gamair)
library(MASS)
library(lme4)
library(nlme)
library(dplyr)
library(ggplot2)

#### Base de datos ####
data("Mmmec")
datos <- Mmmec

#### Modelos ajustados ####
## modelo de Poisson
mod_glm <- glm(deaths ~ uvb + offset(log(expected)), poisson, data = datos)
summary(mod_glm)

## calculo de los residuos 
res.data <- residuos_glm(mod_glm)

## Figura 4.24
res.data %>% 
  ggplot(aes(eta_hat, r_sq)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos de cuantil est")+
  xlab(expression(hat(eta)))

## modelo binomial negativa
mod_nb <- glm.nb(deaths ~ uvb + offset(log(expected)),datos)
summary(mod_nb)

## calculo de los residuos
res.data <- residuos_glm(mod_nb)

## Figura 4.25
res.data %>% 
  ggplot(aes(eta_hat, r_sq)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos de cuantil est")+
  xlab(expression(hat(eta)))

## Modelo lineal generalizado mixto, con cuadratura de Gauss
modgh <- glmer(deaths ~ uvb + offset(log(expected)) + (1|region),
               family = poisson, nAGQ = 25,data = datos)
summary(modgh)

### test para determinar si son necesarios efectos aleatorios
## residuos de desvio para el modelo binomial negativa
rfn <- residuals(mod_nb,type="d")
datos$rfn <- rfn

mod0 <- gls(rfn~1, datos)
mod <- lme(rfn~1, random= ~1|region, datos)

anova(mod0, mod)

## Residuos del MLGM
datos$resid <- residuals(modgh, type = "pearson")
datos$fitted <- fitted(modgh)

## Figura 4.26
datos %>% 
  ggplot(aes(sqrt(fitted), resid)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos de Pearson")+
  xlab(expression(sqrt(hat(mu))))

## Efectos aleatorios estimados
r_2 <- data.frame(ranef(modgh)$region)

## Figura 4.27
qqnorm(r_2$X.Intercept., las=1, pch = 16, xlab = "Cuantiles teóricos",
       ylab = "Cuantiles de la muestra", main = "Q-Q plot normal")
qqline(r_2$X.Intercept., col = "red", lty = 2)
