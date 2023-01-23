#### Paquetes a cargar ####
library(glmtoolbox)
library(statmod)
library(faraway)
library(plotly)
library(ggplot2)
library(dplyr)

#### Base de datos ####
data("skincancer")

#### Grafico relacion incidencia de casos ####
skincancer2 <- skincancer %>% 
  rename(ciudad = city) %>% 
  mutate(rate = (cases/population)*10000)

## Figura 4.17
ggplot(data = skincancer2, aes(age,rate, group = ciudad, color = ciudad)) + 
  geom_line(lwd = 1.3) + theme_bw() + xlab("Rango de edad") + ylab("Casos/10000")


#### modelos ajustados ####
## ciudad como predictor 
mod0 <- glm(cases ~ city, offset=log(population), family=poisson("log"),
            data=skincancer)

## ciudad y edad como predictores
mod1 <- glm(cases ~ city+age, offset=log(population), family=poisson("log"),
            data=skincancer)

## interaccion entre ciudad y edad
mod2 <- glm(cases ~ city*age, offset=log(population), family=poisson("log"),
            data=skincancer)

## comparacion de los modelos con el LRT
anova(mod0,mod1,mod2, test = "Chi")

## resumen modelo seleccionado
summary(mod1)


#### Calculo de los residuos ####
mod1 <- glm(cases ~ city+age, offset=log(population), family=poisson("log"),
            data=skincancer)
res.data <- residuos_glm(mod1)

#### graficos ####
## Figura 4.18
res.data %>% 
  ggplot(aes(eta_hat, r_sq)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos de cuantil est")+
  xlab(expression(hat(eta)))

## Figura 4.19
res.data %>% 
  ggplot(aes(sqrt(mu_hat), r_sq)) + geom_point(cex=2) + theme_gray() + 
  geom_hline(yintercept = 0, col = "red", lty = 2) + ylab("Residuos de cuantil est")+
  xlab(expression(sqrt(hat(mu))))

## Figura 4.20
res.data %>% 
  ggplot(aes(sqrt(mu_hat), abs(r_sd))) + geom_point(cex=2) + theme_gray() + 
  ylab("|Residuos de desvío est|") + xlab(expression(log(hat(mu))))

## Figura 4.21
res.data %>% 
  ggplot(aes(eta_hat, z)) + geom_point(cex=2) + theme_gray() + 
  geom_smooth(method = "lm", se = F, col = "red", lty = 2, lwd = 0.7) + 
  ylab("Respuesta linealizada") + xlab(expression(hat(eta)))

## Figura 4.22
qqnorm(res.data$r_q, las=1, pch = 16,cex=1, xlab = "Cuantiles teóricos",
       ylab = "Cuantiles de la muestra", main = "Q-Q plot normal")
qqline(res.data$r_q, col = "red", lty = 2)


## Figura 4.23
par(mfrow = c(1,2))
halfnorm(res.data$c_i) # halfnorm
plot(mod1,4); abline(h=4/dim(res.data)[1], col = "red", lty=2)

## Modelo sin las observaciones 11 y 15
infl <- c(11,15)
mod.infl <- update(mod1, subset=(-infl))

## Comparacion de coeficientes
round(coef(mod1),3)
round(coef(mod.infl),3)
