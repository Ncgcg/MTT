library(tidyverse)
library(nlme)
library(lsmeans)
library(rstatix)
library(patchwork)

data <- read_csv('MTT.csv')
mod <- aov(data = data, log10(od) ~ h + line)
TukeyHSD(mod)
summary(mod)

data |> filter(h == 210) |> tukey_hsd(od ~ line)
data |> filter(h == 210) |> WRS2::lincon(formula = od ~ line)
data |> filter(h == 210, line == 'TKS4b+') |> shapiro_test(od)
data |> filter(h == 210) |> levene_test(od ~ line)

data |> filter(h == 210) |> 
  ggplot(aes(x = line, y = od, col = line))+
  geom_boxplot()+
  labs(subtitle = 'Б. 9 день інкубації' , x = 'Клітинні лінії', y = '') -> d9

ggplot(data, aes(x = h, y = od, col = line))+
  geom_smooth(alpha = 0.35)+
  labs(subtitle = 'А. Непараметрична регресія' , x = 'Час, год', y = 'OD 570') -> od


##
data |> filter(line == 'MCF7' | line == 'dTKS4') -> dt
data |> filter(line == 'MCF7' | line == 'TKS4L+') -> dt
data |> filter(line == 'MCF7' | line == 'TKS4b+') -> dt
data |> filter(line == 'dTKS4' | line == 'TKS4L+') -> dt
data |> filter(line == 'dTKS4' | line == 'TKS4b+') -> dt
data |> filter(line == 'TKS4L+' | line == 'TKS4b+') -> dt

dt$line <- as.factor(dt$line)


# Fit null model - equal a, b, and c parameters for each group
m0 <- nls(od ~ exp(a+b*h)+c, 
          start = list(a = -6, b = 0.02, c = 0.01),
          data = dt)

# Fit model with differing a, b, and c parameters by group
m1b <- nls(od ~ exp(a+h*b[line])+c, 
          start = list(a = -6, b = rep(0.02, 2), c = 0.01),
          data = dt)

m1a <- nls(od ~ exp(a[line]+h*b)+c, 
           start = list(a = rep(-6, 2), b = 0.02, c = 0.01),
           data = dt)

m1c <- nls(od ~ exp(a+h*b)+c[line], 
           start = list(a = -6, b = 0.02, c = rep(0.01, 2)),
           data = dt)

anova(m0, m1b, m1a, m1c)

###
m1 <- nls(od ~ exp(a+b[line]*h)+c, 
          start = list(a = -6, b = rep(0.02, 2), c = 0.01),
          data = dt)

m2a <- nls(od ~ exp(a[line]+h*b[line])+c, 
          start = list(a=rep(0.01, 2), b = rep(0.02, 2), c = 0.01),
          data = dt)

m2c <- nls(od ~ exp(a + h*b[line])+c[line], 
          start = list(a = 0.01, b = rep(0.02, 2), c = rep(0.01, 2)),
          data = dt)

anova(m0, m2a, m2c)

# Compare with F-test
anova(m0, m1, m2)


