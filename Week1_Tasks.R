setwd("/media/glycosylase/EC6A2F256A2EEBD0/Users/miffka/Documents/!DataMining/Statistics3")

### Задачи недели 1



## Задача 1 - тест на наличие гетероскедастичности.
# Первая колонка - зависимая переменная, остальные - независимые.
# Построить регрессионную модель, вытаскивает остатки, затем строит вспомогательную регрессию
# ЗП ~ остатки^2, и возвращает R^2 этой вспомогательной регрессии.

# Мое решение
hetero_test <- function(df) {
  fit <- lm(formula = reformulate(termlabels = colnames(df)[2:ncol(df)], 
                                  response = colnames(df)[1]), 
            data = df)
  df$resid <- (fit$residuals)^2
  addfit <- lm(formula = reformulate(termlabels = colnames(df)[2:(ncol(df)-1)],
                                     response = colnames(df)[ncol(df)]),
               data = df)
  return(summary(addfit)$r.squared)
}

hetero_test(mtcars)

# Эталон
hetero_test <-  function(test_data){
  fit <- lm(test_data[,1] ~ ., test_data[-1])
  fit_2 <- lm(fit$residuals^2 ~.,test_data[-1])
  summary(fit_2)$r.squared 
}



## Задача 2 - расчет показателя vif.
# Первая колонка - зависимая переменная, остальные - независимые.
# Для каждой независимой переменной расчитать показатель vif

set.seed(42)
test_data <- data.frame(y = rnorm(30, 5), x1 = rnorm(30, 5))
test_data$x2 <- test_data$x1^2

# Мое решение
VIF <- function(df){
  getvif <- function(i){
    f <- lm(df[,i] ~ ., df[-c(1,i)])
    return(1/(1-summary(f)$r.squared))
  }
  ans <- sapply(c(2:ncol(df)), getvif)
  names(ans) <- colnames(df)[-1]
  ans
}

VIF(mtcars)
VIF(test_data)

# Эталонное решение
VIF <- function(data){
  sapply(names(data[-1]), function(name) {
    m = lm(as.formula(paste(name, "~ .")), data[-1])
    r2 = summary(m)$r.squared
    1 / (1 - r2)})
}



## Задача 3 - последовательно исключить предикторы с vif > 10. Вернуть коэффициенты.
# Первая колонка - ЗП, остальные - НП.
# Заводим внутри функцию VIF. Пока максимальное значение вектора VIF > 10,
# исключаем из датасета предиктор с максимальным vif. Затем строим lm и возвращаем коэффиц.

# Мое решение
smart_model <- function(df){
  df <- as.data.frame(df)
  
  VIF <- function(data){
    sapply(names(data[-1]), function(name) {
      m = lm(as.formula(paste(name, "~ .")), data[-1])
      r2 = summary(m)$r.squared
      1 / (1 - r2)})
  }
  
  while (ncol(df) >= 2 && max(VIF(df)) > 10) {
    v <- VIF(df)
    df <- df[, -(which(v == max(v))[1] + 1)]
  }
  
  lm(formula = reformulate(termlabels = colnames(df)[2:ncol(df)], 
                           response = colnames(df)[1]), 
     data = df)$coefficients
}

smart_model(mtcars)
test_data <- as.data.frame(list(y = c(4.63, 5.74, 3.49, 7.23, 3.83, 4.15, 6.93, 5.24, 5.2, 4.47, 4.86, 4.66, 5.41, 6.27, 4.04, 6, 2.96, 4.87, 4.6, 4.4, 3.91, 4.58, 4.75, 4.92, 5.76, 4.7, 4.77, 4.35, 5.51, 4.28), 
                                x1 = c(6.18, 5.19, 4.49, 5.14, 5.64, 5.31, 4.53, 4.97, 4.48, 5.22, 4.33, 4.12, 4.5, 4.26, 5.19, 6.49, 6.53, 6.28, 5.65, 4.46, 4.75, 5.14, 5.05, 5.52, 5.24, 5.22, 4.19, 5.29, 4.48, 4.5), 
                                x2 = c(38.22, 26.97, 20.14, 26.46, 31.8, 28.22, 20.52, 24.72, 20.1, 27.23, 18.79, 16.97, 20.23, 18.12, 26.98, 42.11, 42.59, 39.42, 31.88, 19.85, 22.6, 26.42, 25.45, 30.46, 27.46, 27.29, 17.59, 27.94, 20.04, 20.28)))
smart_model(test_data)

# Эталонное решение
smart_model <- function(d) {
  VIF <- function(d) {
    x <- d[,-1]
    if (ncol(d) > 2) {
      sapply(names(x), function(y) 1/(1-summary(lm(as.formula(paste(y,"~.")), x))$r.sq) )
    } else 0
  }
  
  vif <- VIF(d)
  if (max(vif) > 10) { 
    smart_model(d[, -(which.max(vif) + 1)])
  } else lm(d)$coeff
}



## Задача 4 - поиск наиболее подходящей степени для трансформации независимой переменной.
# Две колонки: y и x.
# Пробегаем по значениям от -2 до 2 с шагом 0,1 и подбираем такое лямбда, чтобы
# было максимальным абсолютное значение корреляции. Вернуть x.

set.seed(42)
test_data <- data.frame(y = rnorm(10, 10, 1), x = rnorm(10, 10, 1))

# Мое решение
transform_x <- function(df, bounds = c(-2, 2), step = 0.1)
{
  vec <- round(seq(bounds[1], bounds[2], step), 1)
  names(vec) <- round(seq(bounds[1], bounds[2], step), 1)
  vc <- c()
  
  for (i in vec){
    if (isTRUE(all.equal(0, i))) a <- log(df$x)
    else a <- i/abs(i)*df$x^i
    df <- cbind(df, a)
    vc <- c(vc, abs(cor( df[,c(1, ncol(df))] ))[1,2])
  }
  
  return(df[,which(vc == max(vc))+2])
}

test_data <- as.data.frame(list(y = c(10.3, 9.33, 9.84, 11.38, 10.07, 11.59, 9.16, 9.89, 12.29, 10.66, 9.54, 8.57, 8.53, 9.1, 8.4, 9.4, 9.51, 9.45, 8.67, 9.29), 
                                x = c(14912.11, 5627.9, 9393.27, 43711.98, 11793.99, 54107.89, 4738.73, 9906.69, 108802.48, 21298.96, 6959.21, 2626.47, 2520.86, 4476.97, 2223.33, 6015.57, 6738.65, 6356.01, 2911.61, 5404.34)))
transform_x(test_data)

# Эталонное решение
transform_x = function(data)
{
  do_transform = function(x, lambd) {
    if (lambd > 0) x ^ lambd else if (lambd < 0) -(x ^ lambd) else log(x) }
  
  x_data = data[,2] 
  y_data = data[,1]
  lambdas = seq(-2, 2, 0.1)
  corrs = sapply(lambdas, 
                 function(lambd) cor(do_transform(x_data, lambd), y_data))
  lambda = lambdas[which.max(abs(corrs))]
  #print(lambda)
  do_transform(x_data, lambda)
}