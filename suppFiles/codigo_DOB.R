# Ajuste DOB freeknotsplines
rm(list=ls())
#Carregar pacotes necessarios
library(readxl)
library(dplyr)
library(tidyr)
library(freeknotsplines)
library(caret)

#Ir para diretorio de trabalho
setwd("C:/Users/Tiago/Drive UNESP/Dissertacao/Vladimir_Bia")

# Leitura dos dados
dados<-read_xlsx("ddob.xlsx",sheet=1)

# Remover dados faltantes
dados<-na.omit(dados)

# Criar matrix de confusao
# Utilizado para avalar inicialmente desempenho do teste UBT
confM<-confusionMatrix(factor(dados$UBT),factor(dados$HIS_F))
confm_m<-confM$table

# Criar bancos de dados para o ajuste das curvas
dadosc<-data.frame(gather(dados,key="tempo",value="DOB",
                          '5','7.5','10','12.5','15','17.5','20','22.5','25','27.5',
                          '30','35','40','45'))

dadoscP<-filter(dadosc,HIS_F=="P")
dadoscN<-filter(dadosc,HIS_F=="N")

{
  fnseeds<-c(993,561) 
  
  fnP<-fit.search.numknots(as.numeric(dadoscP$tempo), as.numeric(dadoscP$DOB), degree=3, minknot = 0, maxknot = 10, seed = fnseeds[1])
  
  fnN<-fit.search.numknots(as.numeric(dadoscN$tempo), as.numeric(dadoscN$DOB), degree=3, minknot = 0, maxknot = 10, seed = fnseeds[2])
}

# Salvar sessao no arquivo freeknots.RData

# Ajuste DOB BASS
rm(list=ls())

# Carregar ambiente de trabalho freeknotsplines
setwd("C:/Users/Tiago/Drive UNESP/Dissertacao/Vladimir_Bia/Defesa")
load("freeknots.RData")


# Carregar pacotes necessarios
library(BASS)


# Ir para diretorio de trabalho
setwd("C:/Users/Tiago/Drive UNESP/Dissertacao/Vladimir_Bia")

# Leitura dos dados
dados<-read_xlsx("ddob.xlsx",sheet=1)

# Remover dados faltantes
dados<-na.omit(dados)

# Criar bancos de dados para o ajuste das curvas
dadosc<-data.frame(gather(dados,key="tempo",value="DOB",
                          '5','7.5','10','12.5','15','17.5','20','22.5','25','27.5',
                          '30','35','40','45'))

# Dados para curvas com exame histologico positivo
dadoscP<-filter(dadosc,HIS_F=="P")

# Dados para curvas com exame histologico negativo
dadoscN<-filter(dadosc,HIS_F=="N")

# Fixar semente para reprodutibilidade
set.seed(138)

# Rodar algoritmo BASS em 3 cadeias (Histologia Positivo)

# Especificacoes simulacao

###
# grau da spline (degree) igual a 3
# thin igual a 100
# burn-in (nburn) igual a 1000
###

{
  p<-bass(xx=data.frame(Tempo=as.numeric(dadoscP$tempo)),dadoscP$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
  p1<-bass(xx=data.frame(Tempo=as.numeric(dadoscP$tempo)),dadoscP$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
  p2<-bass(xx=data.frame(Tempo=as.numeric(dadoscP$tempo)),dadoscP$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
}

# Fixar semente para reprodutibilidade
set.seed(367)

# Rodar algoritmo BASS em 3 cadeias (Histologia Negativo)

# Especificacoes simulacao

###
# grau da spline (degree) igual a 3
# thin igual a 100
# burn-in (nburn) igual a 1000
###

{
  n<-bass(xx=data.frame(Tempo=as.numeric(dadoscN$tempo)),dadoscN$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
  n1<-bass(xx=data.frame(Tempo=as.numeric(dadoscN$tempo)),dadoscN$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
  n2<-bass(xx=data.frame(Tempo=as.numeric(dadoscN$tempo)),dadoscN$DOB,degree=3,thin=100,nmcmc=1001000,nburn = 1000)
}

setwd("C:/Users/Tiago/Drive UNESP/Dissertacao/Vladimir_Bia/Defesa")
load("bass.RData")

x11()
barplot(table(n$nbasis))
barplot(table(n1$nbasis))
barplot(table(n2$nbasis))
barplot(table(c(n$nbasis,n1$nbasis,n2$nbasis)))

barplot(table(p$nbasis))
barplot(table(p1$nbasis))
barplot(table(p2$nbasis))
barplot(table(c(p$nbasis,p1$nbasis,p2$nbasis)))

#load packages
library(readxl)
library(dplyr)
library(tidyr)
library(BASS)
library(coda)
library(RColorBrewer)

library(freeknotsplines)
library(car)

##############################################################################
#Funções criadas para esses dados

#Plotar curva ajustada e intervalos de credibilidade 95% para BASS
pred_line<-function(chain1,chain2,chain3,yourcolor="black",
                    lw=1,lt1=1,lt2=3,mcmcit=1:length(chain1$s2),
                    credibility=0.95,by=1){
  alpha<-(1+credibility)/2
  #Plotar curva
  lines(seq(5,45,by=by),colMeans(rbind(
    predict(chain1,data.frame(Tempo=seq(5,45,by=by)),mcmc.use=mcmcit),
    predict(chain2,data.frame(Tempo=seq(5,45,by=by)),mcmc.use=mcmcit),
    predict(chain3,data.frame(Tempo=seq(5,45,by=by)),mcmc.use=mcmcit)
  )
  
  )
  ,col=yourcolor,lwd=lw,lty=lt1)
  #Criar previsões médias para cada cadeia
  yn_di_ubt<-predict(chain1,seq(5,45,by=by),mcmc.use=mcmcit)
  yn_di_ubta<-predict(chain2,seq(5,45,by=by),mcmc.use=mcmcit)
  yn_di_ubtb<-predict(chain3,seq(5,45,by=by),mcmc.use=mcmcit)
  #Criar multiplicador
  mult <- qnorm(alpha) * sqrt(chain1$s2[mcmcit])
  multa <-qnorm(alpha) * sqrt(chain2$s2[mcmcit])
  multb <- qnorm(alpha) * sqrt(chain3$s2[mcmcit])
  #Criar limites inferiores
  q1 <- apply(yn_di_ubt - mult, 2, quantile, probs = 1-alpha)
  q1a <- apply(yn_di_ubta - multa, 2, quantile, probs = 1-alpha)
  q1b <- apply(yn_di_ubtb - multb, 2, quantile, probs = 1-alpha)
  #Criar limites superiores
  q2 <- apply(yn_di_ubt + mult, 2, quantile, probs = alpha)
  q2a <- apply(yn_di_ubta + multa, 2, quantile, probs = alpha)
  q2b <- apply(yn_di_ubtb + multb, 2, quantile, probs = alpha)
  #Plotar limites inferiores
  lines(seq(5,45,by=by),colMeans(rbind(q1,q1a,q1b)),lty=lt2,col=yourcolor,lwd=lw)
  #Plotar limites superiores
  lines(seq(5,45,by=by),colMeans(rbind(q2,q2a,q2b)),lty=lt2,col=yourcolor,lwd=lw)
}

#Plotar curva ajustada e intervalos de confiança 95% para freeknotsplines
fpred_line<-function(freeknotobj,dadosajuste,yourcolor="gray",lt1=1,lt2=3,lw=1,by=1){
  #Plotar curva ajustada
  lines(seq(5,45,by=by),fitted(freeknotobj,seq(5,45,by=by)),col=yourcolor,lty=lt1,lwd=lw)
  #Criar previsão édia
  fnp_fitted<-fitted(freeknotobj,seq(5,45,by=by))
  #Criar multiplicador
  mult <- 1.96 * sqrt(var(fitted(freeknotobj,as.numeric(dadosajuste$tempo))-as.vector(as.numeric(dadosajuste$DOB))))
  #Criar limite inferior
  q1 <- apply(fnp_fitted - as.numeric(mult), 1, quantile, probs = .025)
  #Criar limite superior
  q2 <- apply(fnp_fitted + as.numeric(mult), 1, quantile, probs = .975)
  #Plotar limite inferior
  lines(seq(5,45,by=by),q1,lty=lt2,col=yourcolor,lwd=lw)
  #Plotar limite superior
  lines(seq(5,45,by=by),q2,lty=lt2,col=yourcolor,lwd=lw) 
}

#Criar objeto de class BASSresiduals de resíduos para BASS com 3 cadeias
BASSresiduals<-function(chain1,chain2,chain3,x_data,y_data,mcmcit=1:length(chain1$s2)){
  #Criar vetor de y estimados para 3 cadeias
  y_est<-colMeans(rbind(
    predict(chain1,data.frame(Tempo=as.numeric(x_data)),mcmc.use=mcmcit),
    predict(chain2,data.frame(Tempo=as.numeric(x_data)),mcmc.use=mcmcit),
    predict(chain3,data.frame(Tempo=as.numeric(x_data)),mcmc.use=mcmcit)))
  #Resíduos
  residuals<-as.numeric(y_data)-as.numeric(y_est) 
  #Criar conjunto de dados de nome residuals com X, Y observado, Y estimado e resíduos
  residuals<-data.frame(matrix(cbind(as.numeric(x_data),as.numeric(y_data),y_est,residuals),ncol=4))
  colnames(residuals)<-c("X","Y Observado","Y Estimado", "Resíduos")
  #Atribuir class: BASSresiduals a objeto
  attr(residuals,'class')<-'BASSresiduals'
  return(residuals)
}

#Criar função plot para objetos da class BASSresiduals
plot.BASSresiduals<-function(obj,bottomtext=NA){
  #Dividir gráfico em 5 painéis 4 grandes em matrix 2x2 em cima
  #1 pequeno em  baixo
  layout(matrix(c(1,2,5,3,4,5),nrow=3),heights = c(1, 1,0.1))
  par(cex.lab=1.5,cex.axis=1.5,cex.main=2,mar=c(5.1,5.1,4.1,2.1))
  plot(obj$X,obj$Resíduos,xlab="Tempo",ylab="Resíduo", main="Resíduo x Tempo")
  plot(obj$'Y Observado',obj$'Y Estimado',xlab="DOB Observado",ylab="DOB Estimado", main="DOB Estimado x Observado")
  abline(a=0,b=1)
  plot(obj$'Y Estimado',obj$Resíduo,xlab="DOB Estimado",ylab="Resíduo", main="Resíduo x DOB Estimado")
  qqPlot(obj$Resíduo,xlab="Quantis da Distribuição Normal",
         ylab="Quantis da Amostra",main="Gráfico Quantil-Quantil da Distribuição Normal")
  par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  text(0,-1,bottomtext,cex=1.5)
}

#Criar objeto de class fnsresiduals de resíduos para FNS
fnsresiduals<-function(fns_model,x_data,y_data){
  #Criar vetor de y estimados para 3 cadeias
  y_est<-fitted(fns_model)
  #Resíduos
  residuals<-residuals(fns_model)
  #Criar conjunto de dados de nome residuals com X, Y observado, Y estimado e resíduos
  residuals<-data.frame(matrix(cbind(as.numeric(x_data),as.numeric(y_data),y_est,residuals),ncol=4))
  colnames(residuals)<-c("X","Y Observado","Y Estimado", "Resíduos")
  #Atribuir class: fnsresiduals a objeto
  attr(residuals,'class')<-'fnsresiduals'
  return(residuals)
}

#Criar função plot para objetos da class fnsresiduals
plot.fnsresiduals<-function(obj,bottomtext=NA){
  #Dividir gráfico em 5 painéis 4 grandes em matrix 2x2 em cima
  #1 pequeno em  baixo
  layout(matrix(c(1,2,5,3,4,5),nrow=3),heights = c(1, 1,0.1))
  par(cex.lab=1.5,cex.axis=1.5,cex.main=2,mar=c(5.1,5.1,4.1,2.1))
  plot(obj$X,obj$Resíduos,xlab="Tempo",ylab="Resíduo", main="Resíduo x Tempo")
  plot(obj$'Y Observado',obj$'Y Estimado',xlab="DOB Observado",ylab="DOB Estimado", main="DOB Estimado x Observado")
  abline(a=0,b=1)
  plot(obj$'Y Estimado',obj$Resíduos,xlab="DOB Estimado",ylab="Resíduo", main="Resíduo x DOB Estimado")
  qqPlot(obj$Resíduos,xlab="Quantis da Distribuição Normal",
         ylab="Quantis da Amostra",main="Gráfico Quantil-Quantil da Distribuição Normal")
  par(fig=c(0,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
  plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
  text(0,-1,bottomtext,cex=1.5)
}

library(convRJMCMC)

cadeias<-c(rep(1,10000),rep(2,10000),rep(3,10000))
iteracoes<-rep(seq(1:10000),3)

x11()
#Gráficos de traços negativos e positivos
traceplot(mcmc.list(mcmc(n$s2),mcmc(n1$s2),mcmc(n2$s2)),
          xlab="Iterações",ylab=expression(sigma**2))
traceplot(mcmc.list(mcmc(p$s2),mcmc(p1$s2),mcmc(p2$s2)),
          xlab="Iterações",ylab=expression(sigma**2))

#Convergência critério Castelloe e Zimmerman (2002)

#Negativos
parametro<-c(n$s2,n1$s2,n2$s2)
modelos<-c(n$nbasis,n1$nbasis,n2$nbasis)
conv_cz<-CZ_ANOVA(parametro,cadeias,modelos,iteracoes)

#Gráficos
x11()
plot(conv_cz)

#Positivos
parametro<-c(p$s2,p1$s2,p2$s2)
modelos<-c(p$nbasis,p1$nbasis,p2$nbasis)
conv_cz<-CZ_ANOVA(parametro,cadeias,modelos,iteracoes)

#Gráficos
x11()
plot(conv_cz)

#Gráficos de autocorrelação

#Negativos
x11()
acf(n$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 1"
acf(n1$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 2"
acf(n2$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 3"

#Positivos
x11()
acf(p$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 1"
acf(p1$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 2"
acf(p2$s2,main=NA,xlab="Atraso",ylab="Autocorrelação")
#"Gráfico de Autocorrelação Cadeia 3"

#Gráficos de Dispersão

#Negativos
x11()
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
plot(dadoscN$tempo,dadoscN$DOB,xlab = "Tempo após ingestão da uréia (minutos)",ylab = "DOB")
pred_line(n,n1,n2,by=0.1)
fpred_line(fnN,dadoscN,by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,1,3),col = c(1,1,"gray","gray"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Curva Ajustada freeknotsplines",
                  "Intervalo de Confiança 95% (freeknotsplines)"),ncol=2)

#versao Cairo
library(Cairo)
setwd("C:/Users/Tiago/Drive UNESP/Dissertacao/Vladimir_Bia/Defesa/DOB")
CairoPDF("Cairodispneg.pdf",width=11.7,height=8.3)
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
plot(dadoscN$tempo,dadoscN$DOB,xlab = "Tempo após ingestão da uréia (minutos)",ylab = "DOB")
pred_line(n,n1,n2,lw=2,by=0.1)
fpred_line(fnN,dadoscN,lw=2,by=0.1,)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,1,3),col = c(1,1,"gray","gray"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Curva Ajustada freeknotsplines",
                  "Intervalo de Confiança 95% (freeknotsplines)"),ncol=2)
dev.off()

#Positivos
x11()
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
ypmin<--40
ypmax<-150
plot(dadoscP$tempo,dadoscP$DOB,xlab = "Tempo",ylab = "DOB"
     ,ylim=c(ypmin,ypmax))
pred_line(p,p1,p2,by=0.1)
fpred_line(fnP,dadoscP,by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,1,3),col = c(1,1,"gray","gray"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Curva Ajustada freeknotsplines",
                  "Intervalo de Confiança 95% (freeknotsplines)"),ncol=2)

#versao Cairo
CairoPDF("Cairodisppos.pdf",width=11.7,height=8.3)
layout(matrix(c(1,2),nrow=2),heights = c(1,0.14))
par(mar=c(5.1,4.1,4.1,2.1))
ypmin<--40
ypmax<-150
plot(dadoscP$tempo,dadoscP$DOB,xlab = "Tempo",ylab = "DOB"
     ,ylim=c(ypmin,ypmax))
pred_line(p,p1,p2,lw=2,by=0.1)
fpred_line(fnP,dadoscP,lw=2,by=0.1)
abline(h=4,lty=5,col="gray")
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,5,1,3),col = c(1,1,"gray","gray","gray"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Ponto de Corte para Classificação UBT (DOB=4\211)",
                  "Curva Ajustada freeknotsplines",
                  "Intervalo de Confiança 95% (freeknotsplines)"
       ),ncol=2)
dev.off()

#Ambos BASS
x11()
CairoPDF("CairodispBASS.pdf",width=11.7,height=8.3)
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
palette(c("black","gray"))
ypmin<--40
ypmax<-150
plot(dadoscN$tempo,dadoscN$DOB,col="gray",pch=1,xlab="Tempo",
     ylab="DOB",ylim=c(ypmin,ypmax))
points(dadoscP$tempo,dadoscP$DOB,col="black",pch=4)
#"Gráfico de Dispersão DOB por Tempo Discriminado por UBT (Peso normal)"  
pred_line(n,n1,n2,2,lt1=1,lt2=2,by=0.1)
pred_line(p,p1,p2,1,lt1=1,lt2=2,by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend(-1,1,pch=c(1,4),col=c("grey","black"),legend=c("Negativo","Positivo"),title="Teste Histológico",ncol=2)
legend(-0.6,1,lty=c(1,2),col=rep("blue"),legend=c("Curva Ajustada","Intervalo de Credibilidade (0.95)"),title="Tipo de Curva",ncol=2)
dev.off()  
# #Ambos freeknotsplines
# x11()
# layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
# par(mar=c(5.1,5.1,4.1,2.1))
# palette(c("black","grey"))
# ypmin<--40
# ypmax<-150
# plot(dadosc$tempo,dadosc$DOB,col=match(dadosc$HIS_F,c("P","N")),
#      pch=ifelse(match(dadosc$HIS_F,c("P","N"))==2,2,4),
#      xlab="Tempo",ylab="DOB",ylim=c(ypmin,ypmax))
# #"Gráfico de Dispersão DOB por Tempo Discriminado por Teste Histológico"  
# fpred_line(fnN,dadoscN,yourcolor="gray")
# fpred_line(fnP,dadoscP,yourcolor="black")
# par(mar=c(rep(0,4)))
# plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
# legend(-1,1,pch=c(2,1),col=c("grey","black"),legend=c("Negativo","Positivo"),title="Teste Histológico",ncol=2)
# legend(-0.6,1,lty=c(1,2),col=rep("blue"),legend=c("Curva Ajustada","Intervalo de Credibilidade (0.95)"),title="Tipo de Curva",ncol=2)

#Ambos freeknotsplines
x11()
CairoPDF("CairodispFN.pdf",width=11.7,height=8.3)
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
palette(c("black","grey"))
ypmin<--40
ypmax<-150
plot(dadoscN$tempo,dadoscN$DOB,col="gray",pch=1,xlab="Tempo",
     ylab="DOB",ylim=c(ypmin,ypmax))
points(dadoscP$tempo,dadoscP$DOB,col="black",pch=4)
#abline(h=4,col="red",lty=5)
#"Gráfico de Dispersão DOB por Tempo Discriminado por Teste Histológico"  
fpred_line(fnN,dadoscN,yourcolor="gray",by=0.1)
fpred_line(fnP,dadoscP,yourcolor="black",by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend(-1,1,pch=c(1,4),col=c("grey","black"),legend=c("Negativo","Positivo"),title="Teste Histológico",ncol=2)
legend(-0.6,1,lty=c(1,2),col=rep("blue"),legend=c("Curva Ajustada","Intervalo de Confiança (0.95)"),title="Tipo de Curva",ncol=2)
dev.off()

#Graficos densidade preditiva 

#Negativos
x11()
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
plot(dadoscN$tempo,dadoscN$DOB,xlab = "Tempo",ylab = "DOB",
     main="Gráfico de Dispersão DOB x Tempo com Densidades Preditivas (Negativo)")
for(i in 1:10000){
  lines(seq(5,45,by=0.1),predict(n,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
  lines(seq(5,45,by=0.1),predict(n1,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
  lines(seq(5,45,by=0.1),predict(n2,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
}
pred_line(n,n1,n2,1,lt1=1,lt2=2,by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,1),col = c("black","black","grey"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Densidades Preditivas (BASS)"),ncol=3)

#Positivos
x11()
layout(matrix(c(1,2),nrow=2),heights = c(1,0.1))
par(mar=c(5.1,5.1,4.1,2.1))
ymin=-55
ymax=150
plot(dadoscP$tempo,dadoscP$DOB,xlab = "Tempo",ylab = "DOB",
     main="Gráfico de Dispersão DOB x Tempo com Densidades Preditivas (Positivo)",
     ylim=c(ymin,ymax))
for(i in 1:10000){
  lines(seq(5,45,by=0.1),predict(p,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
  lines(seq(5,45,by=0.1),predict(p1,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
  lines(seq(5,45,by=0.1),predict(p2,data.frame(Tempo=seq(5,45,by=0.1)),
                                 mcmc.use=i),col="grey")
}
pred_line(p,p1,p2,1,lt1=1,lt2=2,by=0.1)
par(mar=c(rep(0,4)))
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("center",lty=c(1,3,1),col = c("black","black","grey"),
       legend = c("Curva Ajustada BASS (Média a Posteriori)",
                  "Intervalo de Credibilidade 95% (BASS)",
                  "Densidades Preditivas (BASS)"),ncol=3)

library(car)

#Gráficos de resíduos

#Negativos
x11()
residuals_BASSmodel<-BASSresiduals(n,n1,n2,dadoscN$tempo,dadoscN$DOB)
plot(residuals_BASSmodel)

residuals_fnsmodel<-fnsresiduals(fnN,dadoscN$tempo,dadoscN$DOB)
plot(residuals_fnsmodel)

#Positivos
x11()
residuals_BASSmodel<-BASSresiduals(p,p1,p2,dadoscP$tempo,dadoscP$DOB)
plot(residuals_BASSmodel)

residuals_fnsmodel<-fnsresiduals(fnP,dadoscP$tempo,dadoscP$DOB)
plot(residuals_fnsmodel)

length(unique(dadosc$ID))
length(unique(dadoscN$ID))
length(unique(dadoscP$ID))



abline(h=4)
dev.off()
pred_line(p,p1,p2,credibility = 0.95)


# Função para verificar se o limite inferior dos positivos é maior
# que o superior dos positivos
conflimits<-function(chainp1,chainp2,chainp3,chainn1,chainn2,chainn3,yourcolor="black",
                     lw=1,lt1=1,lt2=3,mcmcit=1:length(chainp1$s2),
                     credibility=0.95,by=1){
  alpha<-(1+credibility)/2
  times<-seq(5,45,by=by)
  #Plotar curva
  #Criar previsões médias para cada cadeia dos positivos
  yn_di_ubt<-predict(chainp1,times,mcmc.use=mcmcit)
  yn_di_ubta<-predict(chainp2,times,mcmc.use=mcmcit)
  yn_di_ubtb<-predict(chainp3,times,mcmc.use=mcmcit)
  #Criar multiplicador para positivos
  mult <- qnorm(alpha) * sqrt(chainp1$s2[mcmcit])
  multa <-qnorm(alpha) * sqrt(chainp2$s2[mcmcit])
  multb <- qnorm(alpha) * sqrt(chainp3$s2[mcmcit])
  #Criar limites inferiores
  q1 <- apply(yn_di_ubt - mult, 2, quantile, probs = 1-alpha)
  q1a <- apply(yn_di_ubta - multa, 2, quantile, probs = 1-alpha)
  q1b <- apply(yn_di_ubtb - multb, 2, quantile, probs = 1-alpha)
  #Criar previsões médias para cada cadeia dos positivos
  nyn_di_ubt<-predict(chainn1,times,mcmc.use=mcmcit)
  nyn_di_ubta<-predict(chainn2,times,mcmc.use=mcmcit)
  nyn_di_ubtb<-predict(chainn3,times,mcmc.use=mcmcit)
  #Criar multiplicador para positivos
  nmult <- qnorm(alpha) * sqrt(chainn1$s2[mcmcit])
  nmulta <-qnorm(alpha) * sqrt(chainn2$s2[mcmcit])
  nmultb <- qnorm(alpha) * sqrt(chainn3$s2[mcmcit])
  #Criar limites superiores
  q2 <- apply(nyn_di_ubt + nmult, 2, quantile, probs = alpha)
  q2a <- apply(nyn_di_ubta + nmulta, 2, quantile, probs = alpha)
  q2b <- apply(nyn_di_ubtb + nmultb, 2, quantile, probs = alpha)
  #Plotar limites inferiores
  LI<-colMeans(rbind(q1,q1a,q1b))
  #Plotar limites superiores
  LS<-colMeans(rbind(q2,q2a,q2b))
  result<-data.frame(cbind(times,LI>=LS))
  colnames(result)<-c("Tempos","Maior")
  Tempos<-result$Tempos[which(result$Maior==TRUE)]
  result<-list(Tempos,credibility)
  return(result)
}
#conflimits(p,p1,p2,n,n1,n2,credibility = 0.57,by=0.5)
LimitsBASS<-lapply(seq(0.8,0.57,by=-0.005),FUN=g<-function(x){
  conflimits(p,p1,p2,n,n1,n2,credibility = x,by=0.5)
})

# Exibir resultados do LimitsBASS com os tempos em que os
# intervalos nao se cruzam para cada credibilidade
LimitsBASS

# Obter credibilidade máxima para cada tempo de modo que
# os intervalos não se cruzem
Tempos<-seq(5,45,by=0.5)

Credibilidades<-c(NULL)

for(value in seq(5,45,by=0.5)){
  for(i in seq(1,47)){
    if(is.na(match(value,LimitsBASS[[i]][[1]]))==FALSE){
      li<-LimitsBASS[[i]][[2]]
      break
    }else{
      next
    }
  }
  Credibilidades<-c(Credibilidades,li)
}
remove(li)

# Fazer grafico com a credibilidade
x11()
plot(Tempos,Credibilidades,xaxt='n')
axis(side = 1,at=c(10,16,21,30,40))
abline(v=16,lty=2)
abline(v=21,lty=2)


#Avaliacao UBT tempo 17,5
nrow(subset(dadosc,DOB>4&tempo==17.5&HIS_F=='P'))
nrow(subset(dadosc,DOB>4&tempo==17.5&HIS_F=='N'))
nrow(subset(dadosc,DOB<4&tempo==17.5&HIS_F=='N'))
nrow(subset(dadosc,DOB<4&tempo==17.5&HIS_F=='P'))

#Avaliacao UBT tempo 20
nrow(subset(dadosc,DOB>4&tempo==20&HIS_F=='P'))
nrow(subset(dadosc,DOB>4&tempo==20&HIS_F=='N'))
nrow(subset(dadosc,DOB<4&tempo==20&HIS_F=='N'))
nrow(subset(dadosc,DOB<4&tempo==20&HIS_F=='P'))

