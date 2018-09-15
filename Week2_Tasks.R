setwd("/media/glycosylase/EC6A2F256A2EEBD0/Users/miffka/Documents/!DataMining/Statistics3")
library(ggplot2)
library(lme4)

## Задача 1 - построить график боксплот

exp_data <- read.csv("http://www.bodowinter.com/tutorial/politeness_data.csv")

str(exp_data)
ggplot(exp_data, aes(factor(scenario), frequency, fill = attitude))+
  geom_boxplot()


## Задача 2 - построить графики плотности

ggplot(exp_data, aes(frequency, fill = subject))+
  geom_density(alpha = 0.3)+
  facet_wrap(~gender, dir = "v", strip.position = "right")


## Задача 3, 4, 5 - построить модель frequency oт attitude, случайные эффекты
# Тривиально решаются.
