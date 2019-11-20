######Consumo Livre#######

setwd ("C:\\Users\\Ana\\Documents\\PUC\\Inferencia bayesiana\\Trabalho\\Base")
#Bibliotecas#
library(Kendall)
library(openxlsx)
library(astsa)
library(geoR)
library(deseasonalize)
library(lattice)
library(forecast)

##configração dos plots#
par(mfrow=c(1,1))

##Série de ena mensal##
consumo<-read.xlsx("livre3.xlsx")

#Definindo a série de Ena como série temporal#
A <- ts(consumo[,"ALL"], start=c(2009,9), end = c(2017,4), frequency=12)
x<- ts(consumo[,"AC"], start=c(2009,9), end = c(2016,4), frequency=12)
x1<- ts(consumo[,"PRE"], start=c(2016,5), end = c(2017,4), frequency=12)

x
x1
plot(A, ylab= "MWMed", xlab= "Período")

####Transformação de BOXCOX####

l1=forecast::BoxCox.lambda(x, method=c("guerrero"),lower=-1, upper=1)
z<-forecast::BoxCox(x,0)
z<- log(x)

####sÉRIE DECOMPOSTA e FAC\FACP####
D<- decompose(A)
plot(D)

###FAC E FACP ###
x.acf <- acf(ts(z,frequency = 1), lag.max = 48)
# plot with custom x axis
plot(x.acf, xaxt = "n")
axis(1, seq(12, length(x.acf$acf), 12))

x.pacf <- pacf(ts(z,frequency = 1),lag.max = 48) 
plot(x.pacf, xaxt = "n")
axis(1, seq(12, length(x.pacf$acf), 12))

###FAC E FACP - SAZONAL ###
z1 <- diff(z, 12)
y.acf <- acf(ts(z1, frequency = 1), lag.max = 48)
plot(y.acf, xaxt = "n")
axis(1, seq(12, length(y.acf$acf), 12))

y.pacf <- pacf(ts(z1, frequency = 1), lag.max = 48)
plot(y.pacf, xaxt = "n")
axis(1, seq(12, length(y.pacf$acf), 12))

####Modelo Arima####

y<-sarima(z,1,1,0,1,1,0,12)
y
r<-y$fit$residuals
RRR<- InvBoxCox(r, lambda = 0)
Me<- mean(r)
DPe <- sd(r)
r

####Testes####

##Teste de estacionariedade##
MannKendall(x)

###Teste de normalidade dos resíduos###
tseries::jarque.bera.test (r)
shapiro.test(r)

###Teste de independencia dos resíduos
Box.test(r, lag = 1, type = "Ljung-Box")

density(A)
hist(r)

###Teste de homocedasticidade###
set.seed(2^30-1)
r1<-split(r,sample(r, size=2))
A1<-r1$`0.000642390272406265`
A2<-r1$`0.010911018359534`
y<- c(A1, A2)
group <- as.factor(c(rep(1, length(A1)), rep(2, length(A2))))
car::leveneTest(y,group)

####Sarima~ previsao####

w<- sarima.for(z,12,0,1,0,1,1,0,12)

w1<- InvBoxCox(w$pred, lambda = 0)
W2 = x+RRR
plot(x)
plot(W2)
x
W2
x1
w1
plot(w1)
MAPE <- mean(abs((x1-w1)/x1))
MAPE

####Gráficos####

#Intervalo de confiança
IC <- matrix(NA, nrow=12, ncol=2)
ICs=w$pred + 1.96 * w$se 
ICi= w$pred - 1.96*w$se
IC[,1]=InvBoxCox(ICs,0)
IC[,2]=InvBoxCox(ICi,0)

MX <- matrix(x1, ncol=1)
MW<- matrix(w1, ncol=1)
###GRÁFICO COLORIDO
plot(MX, main= " " ,type= "o", ylab= "MWMed", xlab= "Período", col= "black", ylim = c(11000,20000), lwd=2.5)
lines(IC[,1], col='dimgray', lwd=2.5, lty=1,type= "l")
lines(IC[,2], col='dimgray',lwd=2.5, lty=1, type= "l")
lines(MW, col='red',lwd=2.5, lty=3, type= "o")

legend(4, 20000, legend=c("Intervalo de confiança","Dados observados", "Previsão"),lty=c(1,1,1), lwd=c(2,2,2), col=c("dimgray","black", "red"), pt.cex=10, cex=0.8)

plot (x, xlab = "Ano", ylab = "MWMed", col = "black", xlim = c (2009, 2017), ylim = c(12000,17000))
lines (W2, col = "red")
lines (w1, col = "red", lty = 2)
legend ("topleft", legend = c ("Histórico", "Ajuste", "Previsão"), col = c ("grey","red", "red"), lty = c (1, 1, 2))
                                                                            
####Naive####

ajusteN = naive(z,12)

wN<- ajusteN$mean
Naive1 = InvBoxCox(ajusteN$mean,0)
Naive1
plot (InvBoxCox(z,0), xlab = "Ano", ylab = "MWMed", col = "grey", xlim = c (2009, 2017), ylim = c(12000,17000))
lines (InvBoxCox(ajusteN$fitted,0), col = "red")
lines (InvBoxCox(ajusteN$mean,0), col = "red", lty = 2)
legend ("topleft", legend = c ("Histórico", "Ajuste", "Previsão"), col = c ("grey",     "red", "red"), lty = c (1, 1, 2))

#Intervalo de confiança
ICN <- matrix(NA, nrow=12, ncol=2)
ICNs=ajusteN$mean + 1.96 * ajusteN$model$sd 
ICNi= ajusteN$mean - 1.96*ajusteN$model$sd
ICN[,1]=InvBoxCox(ICNs,0)
ICN[,2]=InvBoxCox(ICNi,0)

MXN <- matrix(x1, ncol=1)
MWN<- matrix(ajusteN$mean, ncol=1)
###GRÁFICO COLORIDO
plot(MXN, main= " " ,type= "o", ylab= "MWMed", xlab= "Período", col= "black", ylim = c(13000,19000), lwd=2.5)
lines(ICN[,1], col='dimgray', lwd=2.5, lty=1,type= "l")
lines(ICN[,2], col='dimgray',lwd=2.5, lty=1, type= "l")
lines(InvBoxCox(MWN,0), col='red',lwd=2.5, lty=3, type= "o")

legend(2, 19000, legend=c("Intervalo de confiança","Dados observados", "Previsão"),lty=c(1,1,1), lwd=c(2,2,2), col=c("dimgray","black", "red"), pt.cex=10, cex=0.8)


####Amortecimento Exponencial Duplo de Brown####

FITh<- HoltWinters(z, gamma = FALSE)
prevB =forecast(FITh, h=12)
FITh
plot (x, xlab = "Ano", ylab = "MWMed", col = "grey", xlim = c (2009, 2017), ylim = c(12000,18000))
lines (FITh$fitted[, 1], col = 2)
lines (prevB$mean, col = 2, lty = 2)
legend ("topleft", legend = c ("Histórico", "Ajuste", "Previsão"), col = c ("grey", 2, 2), lty = c (1, 1, 2))

####Holt winters aditivo####
Fit = HoltWinters(z)
PREVW = forecast(Fit, h=12)

plot (InvBoxCox(z,0), xlab = "Ano", ylab = "MWMed", col = "grey", xlim = c (2009, 2017), ylim = c(12000,18000))
lines (InvBoxCox(Fit$fitted[, 1], 0), col = 2)
lines (InvBoxCox(PREVW$mean, 0), col = 2, lty = 2)
legend ("topleft", legend = c ("Histórico", "Ajuste", "Previsão"), col = c ("grey", 2, 2), lty = c (1, 1, 2))

#Intervalo de confiança
ICH <- matrix(NA, nrow=12, ncol=2)
ICHs= PREVW$mean + 1.96 * PREVW$model$SSE
ICHi= PREVW$mean - 1.96 * PREVW$model$SSE
ICH[,1]=InvBoxCox(ICHs,0)
ICH[,2]=InvBoxCox(ICHi,0)

MXH <- matrix(x1, ncol=1)
MWH<- matrix(PREVW$mean, ncol=1)
###GRÁFICO COLORIDO
plot(MXH, main= " " ,type= "o", ylab= "MWMed", xlab= "Período", col= "black", ylim = c(13000,19000), lwd=2.5)
lines(ICH[,1], col='dimgray', lwd=2.5, lty=1,type= "l")
lines(ICH[,2], col='dimgray',lwd=2.5, lty=1, type= "l")
lines(InvBoxCox(MWH,0), col='red',lwd=2.5, lty=3, type= "o")

legend(2, 19000, legend=c("Intervalo de confiança","Dados observados", "Previsão"),lty=c(1,1,1), lwd=c(2,2,2), col=c("dimgray","black", "red"), pt.cex=10, cex=0.8)

####MODELOS####

MODELS = matrix(NA, nrow = 12, ncol = 4)
MODELS[,1] = x1
MODELS[,2] = w1
MODELS[,3] = InvBoxCox(ajusteN$mean,0)
MODELS[,4] = InvBoxCox(PREVW$mean,0)
MODELS

####testes####

MASE = function(y, y_pred) {
  n = length (y)
  num = (1/n) * sum ( abs (y - y_pred))
  den = (1/(n - 1)) * sum ( abs (y[-1] - lag (y, -1)[-n]))
  mase = round (num/den, 3)
  return (mase)
}
RMSE = function(y, y_pred) {
  rmse = round ( sqrt ( mean ((y - y_pred)^2)), 3)
  return (rmse)
}
MAPE = function(y, y_pred) {
  mape = round (100 * sum ( abs ((y - y_pred)/y))/ length (y), 3)
  return (mape)
}
erros = function(y, y_pred) {
  erros = rbind ( MASE (y, y_pred), MAPE (y, y_pred), RMSE (y, y_pred))
  rownames (erros) = c ("MASE", "MAPE", "RMSE")
  return (erros)
}

MASE1 = MASE(x1,w1)
MASE2 = MASE(x1,MODELS[,3])
MASE3 = MASE(x1,MODELS[,4])

MAPE1 = MAPE(x1,w1)
MAPE2 = MAPE(x1,MODELS[,3])
MAPE3 = MAPE(x1,MODELS[,4])

RMSE1 = RMSE(x1,w1)
RMSE2 = RMSE(x1,MODELS[,3])
RMSE3 = RMSE(x1,MODELS[,4])


####aNÁLISE####

#### Análise dos resultados####

am<- matrix(x, 12)#matriz com os dados por mês na linha
aa <- matrix(NA, nrow = 12)
for (i in 1:12){
  aa[i,] <- mean(am[i,])
}

ks.test(aa, w1, alternative = c("two.sided"))  # previsão SARIMA e série histórica
ks.test(aa, MODELS[,3], alternative = c("two.sided"))  # previsão NAIVE e SÉRIE HISTÓRICA
ks.test(aa, MODELS[,4], alternative = c("two.sided"))  # previsão winters e SÉRIE HISTÓRICA

ks.test(x1, w1, alternative = c("two.sided"))  # previsão SARIMA e dados previstos
ks.test(x1, MODELS[,3], alternative = c("two.sided"))  # previsão NAIVE e dados previstos
ks.test(x1, MODELS[,4], alternative = c("two.sided"))  # previsão winters e dados previstos


t.test(aa, w1, alternative = c("two.sided")) # previsão sarima e série histórica
t.test(aa, MODELS[,3], alternative = c("two.sided")) # previsão naive e série histórica
t.test(aa, MODELS[,4], alternative = c("two.sided")) # previsão winters e série histórica

t.test(x1, w1, alternative = c("two.sided")) # previsão sarima e dados previstos
t.test(x1, MODELS[,3], alternative = c("two.sided")) # previsão naive e dados previstos
t.test(x1, MODELS[,4], alternative = c("two.sided")) # previsão winters e dados previstos


LEV1<- c(aa, w1)
group1 <- as.factor(c (rep(1, length(aa)), rep(2, length(w1)))) 
car::leveneTest(LEV1,group1)   # previsão SARIMA e série histórica

LEV2<- c(aa, MODELS[,3])
group2 <- as.factor(c(rep(1, length(aa)), rep(2, length(MODELS[,3]))))
car::leveneTest(LEV2,group2) # previsão NAIVE e série histórica

LEV3<- c(aa, MODELS[,4])
group3 <- as.factor(c(rep(1, length(aa)), rep(2, length(MODELS[,4]))))
car::leveneTest(LEV3,group3) # previsão WINTERS e série histórica

LEV11<- c(x1, w1)
group11 <- as.factor(c (rep(1, length(x1)), rep(2, length(w1)))) 
car::leveneTest(LEV11,group1)   # previsão SARIMA e PREVISTO

LEV21<- c(x1, MODELS[,3])
group21 <- as.factor(c(rep(1, length(x1)), rep(2, length(MODELS[,3]))))
car::leveneTest(LEV21,group21) # previsão NAIVE e PREVISTO

LEV31<- c(x1, MODELS[,4])
group31 <- as.factor(c(rep(1, length(x1)), rep(2, length(MODELS[,4]))))
car::leveneTest(LEV31,group31) # previsão WINTERS e PREVISTO
