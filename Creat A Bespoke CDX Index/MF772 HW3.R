#MF772 HW3
#Sike Yang
#a
library(credule)
#issue 1
#set risk free rate as 2% and recovery rate as 40%
CDScurve1 = bootstrapCDS(yieldcurveTenor=c(1,2,3),
                        yieldcurveRate=rep(0.02,3),
                        cdsTenors=c(1,2,3),
                        cdsSpreads=c(0.015,0.026,0.037),
                        recoveryRate=rep(0.4,3),
                        numberPremiumPerYear=1,
                        numberDefaultIntervalPerYear=1)
#issue 2
CDScurve2 = bootstrapCDS(yieldcurveTenor=c(1,2,3),
                         yieldcurveRate=rep(0.02,3),
                         cdsTenors=c(1,2,3),
                         cdsSpreads=c(0.005,0.006,0.007),
                         recoveryRate=rep(0.4,3),
                         numberPremiumPerYear=1,
                         numberDefaultIntervalPerYear=1)

#issue 3
CDScurve3 = bootstrapCDS(yieldcurveTenor=c(1,2,3),
                         yieldcurveRate=rep(0.02,3),
                         cdsTenors=c(1,2,3),
                         cdsSpreads=c(0.13,0.18,0.24),
                         recoveryRate=rep(0.4,3),
                         numberPremiumPerYear=1,
                         numberDefaultIntervalPerYear=1)

#issue 4
CDScurve4 = bootstrapCDS(yieldcurveTenor=c(1,2,3),
                         yieldcurveRate=rep(0.02,3),
                         cdsTenors=c(1,2,3),
                         cdsSpreads=c(0.00498,0.0088,0.01),
                         recoveryRate=rep(0.4,4),
                         numberPremiumPerYear=1,
                         numberDefaultIntervalPerYear=1)

#issue 5
CDScurve5 = bootstrapCDS(yieldcurveTenor=c(1,2,3),
                         yieldcurveRate=rep(0.02,3),
                         cdsTenors=c(1,2,3),
                         cdsSpreads=c(0.022,0.023,0.024),
                         recoveryRate=rep(0.4,4),
                         numberPremiumPerYear=1,
                         numberDefaultIntervalPerYear=1)

#b
lossfunction <- function(CDScurve,recovery_rate){
  loss = CDScurve$survprob
  loss[1] = (1-recovery_rate)*(1-CDScurve$survprob[1])
  for (i in 2:length(CDScurve$survprob)){
    loss[i] <- (1-recovery_rate)*(CDScurve$survprob[i-1]-CDScurve$survprob[i])
  }
  return(loss)
}
losscurve1 <- lossfunction(CDScurve1,0.4)
losscurve2 <- lossfunction(CDScurve2,0.4)
losscurve3 <- lossfunction(CDScurve3,0.4)
losscurve4 <- lossfunction(CDScurve4,0.4)
losscurve5 <- lossfunction(CDScurve5,0.4)

#c
library(Ryacas0)
s <- Sym("s")
cal_rpv01 <- function(CDScurve){
  rpv01 <- 0.5*1*(1/(1+0.02))*(1+CDScurve$survprob[1])
  for (i in 2:length(CDScurve$survprob)){
    rpv01 <- rpv01+0.5*1*(1/(1+0.02)^i)*(CDScurve$survprob[i-1]+CDScurve$survprob[i])
  }
  return(rpv01)
}
rpv01_1 <- cal_rpv01(CDScurve1)
rpv01_2 <- cal_rpv01(CDScurve2)
rpv01_3 <- cal_rpv01(CDScurve3)
rpv01_4 <- cal_rpv01(CDScurve4)
rpv01_5 <- cal_rpv01(CDScurve5)

premium_leg <- (rpv01_1+rpv01_2+rpv01_3+rpv01_4+rpv01_5)*s/5

protection <- function(CDScurve){
  leg <- (1/(1+0.02))*(1-CDScurve$survprob[1])
  for (i in 2:length(CDScurve$survprob)){
    leg <- leg+(1/(1+0.02)^i)*(CDScurve$survprob[i-1]-CDScurve$survprob[i])
  }
  leg <- leg*0.6
}

pleg1 <- protection(CDScurve1)
pleg2 <- protection(CDScurve2)
pleg3 <- protection(CDScurve3)
pleg4 <- protection(CDScurve4)
pleg5 <- protection(CDScurve5)

protection_leg <- (pleg1+pleg2+pleg3+pleg4+pleg5)/5

S <- Eval(Solve(premium_leg == protection_leg, s, 1, 0.00001))
#the fair spread is 264.85bps

#d
loss_tranches <- function(lt,k1,k2){
  R <- 0.4
  N <- 5
  k2 <- min(k2,N-R/N)
  loss = (max(lt-k1,0) - max(lt-k2,0))/(k2-k1)
  return(loss)
}

tranche_premium_leg <- function(L,K1,K2,r){
  M <- 1-((loss_tranches(L[1],K1,K2)+loss_tranches(0,K1,K2))/2)
  leg <- 0.5*(1/(1+r))*M
  for (i in 2:length(L)){
    M <- 1-((loss_tranches(L[i],K1,K2)+loss_tranches(L[i-1],K1,K2))/2)
    leg <- leg+0.5*(1/(1+r)^i)*M
  }
  return(leg)
}

tranche_protect_leg <- function(L,K1,K2,r){
  minus <- loss_tranches(L[1],K1,K2)-loss_tranches(0,K1,K2)
  leg <- (1/(1+r))*minus
  for (i in 2:length(L)){
    minus <- loss_tranches(L[i],K1,K2)-loss_tranches(L[i-1],K1,K2)
    leg <- leg + (1/((1+r)^i))*minus
  }
  return(leg)
}

tranche_spread <- function(L,K1,K2,r){
  pre_leg <- tranche_premium_leg(L,K1,K2,r)
  prot_leg <- tranche_protect_leg(L,K1,K2,r)
  s <- Sym("s")
  pre_leg = pre_leg * s
  spread = Eval(Solve(pre_leg == prot_leg, s, 1, 0.00001))
  spread = as.numeric(gsub("[^0-9.]", "", spread))
  return(spread)
}

x_crit <- function(CDScurve){
  return(qnorm(1-CDScurve$survprob))
}

pj_m <- function(xcrit, rho, m) {
  value <- (xcrit - sqrt(rho) * m) / sqrt(1 - rho)
  return(pnorm(value))
}

loss_distribution <- function(CDScurve,rho){  
  m <- 0.15
  x <- x_crit(CDScurve)
  pj <- pj_m(x,rho,m)
  bj <- c(0,0,0)
  lj <- 0.2 * 0.6
  for (i in 1:3){
    bj[i] <- 0*(1-pj[i])+lj*pj[i]
  }
  return(bj)
}

lt1 <- loss_distribution(CDScurve1,0.3)
lt2 <- loss_distribution(CDScurve2,0.3)
lt3 <- loss_distribution(CDScurve3,0.3)
lt4 <- loss_distribution(CDScurve4,0.3)
lt5 <- loss_distribution(CDScurve5,0.3)

lt <- lt1+lt2+lt3+lt4+lt5

tranche_spread(lt,0,0.05,0.02)
tranche_spread(lt,0.05,0.15,0.02)
tranche_spread(lt,0.15,0.30,0.02)
tranche_spread(lt,0.30,1,0.02)

#e
#RV01
tranche_premium_leg(lt,0,0.05,0.02)
tranche_premium_leg(lt,0.05,0.15,0.02)
tranche_premium_leg(lt,0.15,0.30,0.02)
tranche_premium_leg(lt,0.30,1,0.02)

#DV01
dv01 <- function(L,K1,K2,r){
  old_s <- tranche_spread(L,K1,K2,r)
  old_s <- as.numeric(gsub("[^0-9.]", "", old_s))
  new_s <- tranche_spread(L,K1,K2,r+1/10000)
  new_s <- as.numeric(gsub("[^0-9.]", "", new_s))
  return((new_s-old_s)/(1/10000))
}
dv01(lt,0,0.05,0.02)
dv01(lt,0.05,0.15,0.02)
dv01(lt,0.15,0.3,0.02)
dv01(lt,0.30,1,0.02)

#f
spread = c(tranche_spread(lt,0,0.05,0.02),
           tranche_spread(lt,0.05,0.15,0.02),
           tranche_spread(lt,0.15,0.30,0.02),
           tranche_spread(lt,0.30,1,0.02))
tranches <- data.frame(attachment = c(0, 0.05, 0.15, 0.30),
                       detachment = c(0.05, 0.15, 0.30, 1))
LGD <- 0.6
weight <- 0.2
JtD_matrix <- matrix(0, nrow = nrow(tranches), ncol = 5)
for (i in 1:5) {
  loss <- LGD * weight
  remaining_loss <- loss
  for (j in 1:nrow(tranches)) {
    if (remaining_loss > 0) {
      tranche_loss <- min(remaining_loss, tranches$detachment[j] - tranches$attachment[j])
      JtD_matrix[j, i] <- tranche_loss
      remaining_loss <- remaining_loss - tranche_loss
    }
  }
}
print(JtD_matrix)

#g
lt1 <- loss_distribution(CDScurve1,0.5)
lt2 <- loss_distribution(CDScurve2,0.5)
lt3 <- loss_distribution(CDScurve3,0.5)
lt4 <- loss_distribution(CDScurve4,0.5)
lt5 <- loss_distribution(CDScurve5,0.5)

lt <- lt1+lt2+lt3+lt4+lt5

tranche_spread(lt,0,0.05,0.02)
tranche_spread(lt,0.05,0.15,0.02)
tranche_spread(lt,0.15,0.30,0.02)
tranche_spread(lt,0.30,1,0.02)
#Go Long on Higher Tranches and Short Lower Tranches

#h
#Yes
#The intent is typically to minimize exposure to broad market moves 
#and isolate soread between the long and short positions.
#1.Determine the Dollar Amount to Trade
#2.Short the Equity Tranche
#3.Go Long on Senior Tranches
#4.Balancing with Other Tranches/Instruments








