rm(list=ls())

# libraries #####

library(writexl)
library(data.table)
library(deSolve)
library(tidyverse)
library(dplyr)
library(stringr) 
library(readxl)
library(fitdistrplus)
library(logspline)
library(pracma)
library(stats)
library(evd)
library(nleqslv)
library(ggpubr)

# Create Skill Distribution #####

x <- seq(0, 10, length.out = 100)

z <- dgpd(x, loc = 0, scale = 1, shape = 1, log = FALSE)

plot(x, z, type = "l", xlab = "x", ylab = "Density")

percent <- data.frame(sapply(0:9, function(i) integrate(dgpd, lower=i, upper=i+1, loc = 0, scale = 1, shape = 1, log = FALSE)$value))
names(percent)[1] <- "dis"

sum <- sum(percent$dis) # round and change a little bit such that sum is equal to 1
percent = percent %>% 
  mutate(dis3 = dis * (1-sum) + dis,
         dis4 = round(dis3, digits = 2))
sum(percent$dis4)

Kap_total = 1 # set total skill in the fishery equal to one

percent = percent %>% 
  mutate(Kap = dis4 * Kap_total)
sum(percent$Kap)

# Define Fixed Variables #####

r=0.69 #intrinsic growth rate
k=1095 #carrying capacity
a=0.5 #alpha 
p=1 #price
w=151.1 #cost parameter
q=1 #catchability coeffcient 


# Open Access #####
# According to the model equations starting at a stock level of 1 thousand metric tons
# For 102 periods of time to get to steady state

oa <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(oa) <- c('time', 'S', 'Htotal')
oa$time <- c(1:102)
oa$S[1] = 1
oa$Htotal[1] = 0

for (t in 1:101) {
  
  oa$Htotal[t] = ((q*oa$S[t])^(1/a)) * (((1-a)*(p/w))^((1-a)/a)) * Kap_total
  
  oa$S[t+1] = oa$S[t] + (r*oa$S[t]*(1-(oa$S[t]/k))) - oa$Htotal[t]
  
}

oa <- oa %>% 
  mutate(H = sapply(1:10, function (i)  ((q*S)^(1/a)) * (((1-a)*(p/w))^((1-a)/a)) *  percent$Kap[i]),
         Eff= (((1-a)*p*q*S/w)^(1/a)) * Kap_total,
         H_sum = H[,1] + H[,2]+ H[,3] + H[,4] + H[,5] + H[,6] +  H[,7] + H[,8] + H[,9] + H[,10],
         profit = sapply(1:10, function (i)
           ((((1-a)^((1-a)/a)) - ((1-a)^(1/a))) * ((p*q*S)^(1/a))* (w^((a-1)/a))) *  percent$Kap[i]),
         share1 = H[,1]/Htotal)

oa = oa %>% relocate(H_sum, .after=Htotal)

oa = oa %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate of the stock in percent
oa$gS[1]=1

oa_path <- oa[abs(oa$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state

# calculate steady state annuity profits

oa_st_st <- oa[abs(oa$gS) >= 0.001, ]
oa_st_st <- tail(oa_st_st, n=1) 
profit <- data.frame(sapply(1:10, function (i) profit <- (oa_st_st$profit[,i])))
colnames(profit) <- c("profit_oa") 

profit = profit %>%
  mutate(sum_oa = sum(profit_oa))

# Case 1: TAC_cap #####
# TAC set by regulator according to harvest control rule and distribution of quota shares based on skill distribution

taccap <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(taccap) <- c('time', 'S', 'Htotal')
taccap$time <- c(1:102)
taccap$S[1] = 1
taccap$Htotal[1] = 0

control = r/2 #linear harvest control for MSY

percent = percent %>%
  mutate(share = Kap/Kap_total) #share as skill distribution 


for (t in 1:101) {
  
  taccap$Htotal[t] = control * taccap$S[t] 
  
  taccap$S[t+1] = taccap$S[t] + (r*taccap$S[t]*(1-(taccap$S[t]/k))) - taccap$Htotal[t]
  
}


taccap <- taccap %>% 
  mutate(H = sapply(1:10, function (i) percent$share[i] * Htotal),
         Eff = sapply(1:10, function (i) (percent$share[i]*r*S/(2*q*S))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         H_sum = H[,1] + H[,2]+ H[,3] + H[,4] + H[,5] + H[,6] +  H[,7] + H[,8] + H[,9] + H[,10],
         Eff_sum = Eff[,1] + Eff[,2]+ Eff[,3] + Eff[,4] + Eff[,5] + Eff[,6] +  Eff[,7] + Eff[,8] + Eff[,9] + Eff[,10],
         profit = sapply(1:10, function (i)
           p * percent$share[i] * Htotal - w * (percent$share[i]*r/(2*q))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         share1 = H[,1]/Htotal)
taccap = taccap %>% relocate(H_sum, .after=Htotal)
taccap <- taccap %>% 
  mutate(o = p - (w*Eff_sum**a)/((1-a)*q*Kap_total**a*S))


taccap = taccap %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate of stock in percent 

taccap_path <- taccap[abs(taccap$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state

# calculate steady state annuity profits 

taccap_st_st <- taccap[abs(taccap$gS) >= 0.001, ]
taccap_st_st <- tail(taccap_st_st, n=1) 

profit_taccap <- data.frame(sapply(1:10, function (i) profit <- (taccap_st_st$profit[,i])))
colnames(profit_taccap) <- c("profit_taccap") 
profit <- bind_cols(profit, profit_taccap)

profit = profit %>%
  mutate(sum_taccap = sum(profit_taccap))



# Case 1: No Trade and Skill Allocation #####
# TAC set by regulator according to harvest control rule and distribution of quota shares based on skill distribution

tac <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(tac) <- c('time', 'S', 'Htotal')
tac$time <- c(1:102)
tac$S[1] = 1
tac$Htotal[1] = 0

control = r/2 #linear harvest control for MSY

percent = percent %>%
  mutate(share = Kap/Kap_total) #share as capital distribution 

for (t in 1:101) {
  
  tac$Htotal[t] = control * tac$S[t] 
  
  tac$S[t+1] = tac$S[t] + (r*tac$S[t]*(1-(tac$S[t]/k))) - tac$Htotal[t]
  
}


tac <- tac %>% 
  mutate(H = sapply(1:10, function (i) percent$share[i] * Htotal),
         Eff = sapply(1:10, function (i) (percent$share[i]*r*S/(2*q*S))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         H_sum = H[,1] + H[,2]+ H[,3] + H[,4] + H[,5] + H[,6] +  H[,7] + H[,8] + H[,9] + H[,10],
         Eff_sum = Eff[,1] + Eff[,2]+ Eff[,3] + Eff[,4] + Eff[,5] + Eff[,6] +  Eff[,7] + Eff[,8] + Eff[,9] + Eff[,10],
         profit = sapply(1:10, function (i)
           p * percent$share[i] * Htotal - w * (percent$share[i]*r/(2*q))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         share1 = H[,1]/Htotal)
tac = tac %>% relocate(H_sum, .after=Htotal)
tac <- tac %>% 
  mutate(o = p - (w*Eff_sum**a)/((1-a)*q*Kap_total**a*S))

tac = tac %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate of stock in percent
tac$gS[1]=1

tac_path <- tac[abs(tac$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state

# calculate steady state annuity profits 

tac_st_st <- tac[abs(tac$gS) >= 0.001, ]
tac_st_st <- tail(tac_st_st, n=1) 

profit_tac <- data.frame(sapply(1:10, function (i) profit <- (tac_st_st$profit[,i])))
colnames(profit_tac) <- c("profit_tac") 
profit <- bind_cols(profit, profit_tac)

profit = profit %>%
  mutate(sum_tac = sum(profit_tac))


# Case 2: Full Trade and Skill Allocation #####

percent = percent %>%
  mutate(share = Kap/Kap_total) #share as capital distribution 

m = 0 #no transaction cost

trade_no_cost <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(trade_no_cost) <- c('time', 'S', 'Htotal')
trade_no_cost$time <- c(1:102)

trade_no_cost$S[1] = 1 
trade_no_cost$Htotal[1] = 0

# calculate willingness to pay based on demand function as explained in the model for limited transferability: 


for (t in 1:101) {
  
  trade_no_cost$Htotal[t] = control * trade_no_cost$S[t] 
  
  trade_no_cost$S[t+1] = trade_no_cost$S[t] + (r*trade_no_cost$S[t]*(1-(trade_no_cost$S[t]/k))) - trade_no_cost$Htotal[t]
  
}

trade_no_cost <- trade_no_cost %>% 
  mutate(a = (p + m * percent$share[1] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[1]**(-1)),
         b = (p + m * percent$share[2] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[2]**(-1)),
         c = (p + m * percent$share[3] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[3]**(-1)),
         d = (p + m * percent$share[4] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[4]**(-1)),
         e = (p + m * percent$share[5] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[5]**(-1)),
         f = (p + m * percent$share[6] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[6]**(-1)),
         g = (p + m * percent$share[7] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[7]**(-1)),
         h = (p + m * percent$share[8] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[8]**(-1)),
         i = (p + m * percent$share[9] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[9]**(-1)),
         j = (p + m * percent$share[10] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[10]**(-1)),
         sum1 = a+b+c+d+e+f+g+h+i+j,
         a1 = (1 / (m + 2*w*S**(-2)*percent$Kap[1]**(-1))),
         a2 = (1 / (m + 2*w*S**(-2)*percent$Kap[2]**(-1))),
         a3 = (1 / (m + 2*w*S**(-2)*percent$Kap[3]**(-1))),
         a4 = (1 / (m + 2*w*S**(-2)*percent$Kap[4]**(-1))),
         a5 = (1 / (m + 2*w*S**(-2)*percent$Kap[5]**(-1))),
         a6 = (1 / (m + 2*w*S**(-2)*percent$Kap[6]**(-1))),
         a7 = (1 / (m + 2*w*S**(-2)*percent$Kap[7]**(-1))),
         a8 = (1 / (m + 2*w*S**(-2)*percent$Kap[8]**(-1))),
         a9 = (1 / (m + 2*w*S**(-2)*percent$Kap[9]**(-1))),
         a10 = (1 / (m + 2*w*S**(-2)*percent$Kap[10]**(-1))),
         sum2 = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10,
         o = (sum1 - Htotal) / sum2
  )

trade_no_cost <- trade_no_cost %>% 
  dplyr::select(time, S, Htotal, o)

# calculate Harvest: 

for (t in 1:101) {
  
  #fisher 1
  n = (1/(1-a))*w*(1/(q*percent$Kap[1]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[1]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H1[t] = harvest$x[1]
  
  # fisher 2
  n = (1/(1-a))*w*(1/(q*percent$Kap[2]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[2]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H2[t] = harvest$x[1]
  
  # fisher 3
  n = (1/(1-a))*w*(1/(q*percent$Kap[3]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[3]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H3[t] = harvest$x[1]
  
  # fisher 4
  n = (1/(1-a))*w*(1/(q*percent$Kap[4]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[4]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H4[t] = harvest$x[1]
  
  # fisher 5
  n = (1/(1-a))*w*(1/(q*percent$Kap[5]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[5]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H5[t] = harvest$x[1]
  
  # fisher 6
  n = (1/(1-a))*w*(1/(q*percent$Kap[6]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[6]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H6[t] = harvest$x[1]
  
  # fisher 7
  n = (1/(1-a))*w*(1/(q*percent$Kap[7]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[7]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H7[t] = harvest$x[1]
  
  # fisher 8
  n = (1/(1-a))*w*(1/(q*percent$Kap[8]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[8]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H8[t] = harvest$x[1]
  
  # fisher 9
  n = (1/(1-a))*w*(1/(q*percent$Kap[9]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[9]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H9[t] = harvest$x[1]
  
  # fisher 10
  n = (1/(1-a))*w*(1/(q*percent$Kap[10]**a*trade_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_no_cost$o[t] - n*x**b - m*x + m*percent$share[10]*trade_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_no_cost$H10[t] = harvest$x[1]
  
  
}

# calculate indiviudal profits 

trade_no_cost <- trade_no_cost %>% 
  mutate(H_sum = H1 + H2+ H3 + H4 + H5 + H6 +  H7 + H8 + H9 + H10,
         profit1 =
           p*H1 - w * (H1/(q*percent$Kap[1]**a*S))**(1/(1-a))
         - o * (H1 - percent$share[1] * Htotal) 
         - m/2 * (H1 - percent$share[1] * Htotal)**2,
         profit2 =
           p*H2 - w * (H2/(q*percent$Kap[2]**a*S))**(1/(1-a)) 
         - o * (H2 - percent$share[2] * Htotal) 
         - m/2 * (H2 - percent$share[2] * Htotal)**2,
         profit3 =
           p*H3 - w * (H3/(q*percent$Kap[3]**a*S))**(1/(1-a))
         - o * (H3 - percent$share[3] * Htotal) 
         - m/2 * (H3 - percent$share[3] * Htotal)**2,
         profit4 =
           p*H4 - w * (H4/(q*percent$Kap[4]**a*S))**(1/(1-a)) 
         - o * (H4 - percent$share[4] * Htotal) 
         - m/2 * (H4 - percent$share[4] * Htotal)**2,
         profit5 =
           p*H5 - w * (H5/(q*percent$Kap[5]**a*S))**(1/(1-a))
         - o * (H5 - percent$share[5] * Htotal) 
         - m/2 * (H5 - percent$share[5] * Htotal)**2,
         profit6 =
           p*H6 - w * (H6/(q*percent$Kap[6]**a*S))**(1/(1-a))
         - o * (H6 - percent$share[6] * Htotal) 
         - m/2 * (H6 - percent$share[6] * Htotal)**2,
         profit7 =
           p*H7 - w * (H7/(q*percent$Kap[7]**a*S))**(1/(1-a))
         - o * (H7 - percent$share[7] * Htotal) 
         - m/2 * (H7 - percent$share[7] * Htotal)**2,
         profit8 =
           p*H8 - w * (H8/(q*percent$Kap[8]**a*S))**(1/(1-a))
         - o * (H8 - percent$share[8] * Htotal) 
         - m/2 * (H8 - percent$share[8] * Htotal)**2,
         profit9 =
           p*H9 - w * (H9/(q*percent$Kap[9]**a*S))**(1/(1-a))
         - o * (H9 - percent$share[9] * Htotal) 
         - m/2 * (H9 - percent$share[9] * Htotal)**2,
         profit10 =
           p*H10 - w * (H10/(q*percent$Kap[10]**a*S))**(1/(1-a))
         - o * (H10 - percent$share[10] * Htotal) 
         - m/2 * (H10 - percent$share[10] * Htotal)**2,
         mcost1 = 
           m / 2 * (H1-percent$share[1]*Htotal)**2,
         mcost2 = 
           m / 2 * (H2-percent$share[2]*Htotal)**2,
         mcost3 =
           m / 2 * (H3-percent$share[3]*Htotal)**2,
         mcost4 = 
           m / 2 * (H4-percent$share[4]*Htotal)**2,
         mcost5 = 
           m / 2 * (H5-percent$share[5]*Htotal)**2,
         mcost6 = 
           m / 2 * (H6-percent$share[6]*Htotal)**2,
         mcost7 =
           m / 2 * (H7-percent$share[7]*Htotal)**2,
         mcost8 = 
           m / 2 * (H8-percent$share[8]*Htotal)**2,
         mcost9 = 
           m / 2 * (H9-percent$share[9]*Htotal)**2,
         mcost10 = 
           m / 2 * (H10-percent$share[10]*Htotal)**2,
         msum = mcost1 + mcost2 + mcost3 +mcost4 +mcost5+mcost6+mcost7+mcost8+mcost9+mcost10
  )


trade_no_cost = trade_no_cost %>% relocate(H_sum, .after=Htotal)

trade_no_cost = trade_no_cost %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate in percent of stock 
trade_no_cost$gS[1] = 1

trade_no_cost_path <- trade_no_cost[abs(trade_no_cost$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state


trade_no_cost_st_st <- trade_no_cost[abs(trade_no_cost$gS) >= 0.001, ]
trade_no_cost_st_st <- tail(trade_no_cost_st_st, n=1) 

# calculate steady state profits

profit_no_cost <- c(trade_no_cost_st_st$profit1[1], trade_no_cost_st_st$profit2[1], trade_no_cost_st_st$profit3[1], trade_no_cost_st_st$profit4[1], 
                    trade_no_cost_st_st$profit5[1], trade_no_cost_st_st$profit6[1], trade_no_cost_st_st$profit7[1], trade_no_cost_st_st$profit8[1],
                    trade_no_cost_st_st$profit9[1], trade_no_cost_st_st$profit10[1])
profit_no_cost <- data.frame(profit_no_cost)
profit <- bind_cols(profit, profit_no_cost)

profit = profit %>%
  mutate(sum_no_cost = sum(profit_no_cost))

# Case 3: Limited Trade and Skill Allocation #####


percent = percent %>%
  mutate(share = Kap/Kap_total) #share as capital distribution 

m = 0.1 #transaction cost parameter = 0.1

trade_cost <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(trade_cost) <- c('time', 'S', 'Htotal')
trade_cost$time <- c(1:102)

trade_cost$S[1] = 1 
trade_cost$Htotal[1] = 0


for (t in 1:101) {
  
  trade_cost$Htotal[t] = control * trade_cost$S[t] 
  
  trade_cost$S[t+1] = trade_cost$S[t] + (r*trade_cost$S[t]*(1-(trade_cost$S[t]/k))) - trade_cost$Htotal[t]
  
}

# calculate willingness to pay based on demand function as explained in the model for limited transferability: 

trade_cost <- trade_cost %>% 
  mutate(a = (p + m * percent$share[1] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[1]**(-1)),
         b = (p + m * percent$share[2] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[2]**(-1)),
         c = (p + m * percent$share[3] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[3]**(-1)),
         d = (p + m * percent$share[4] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[4]**(-1)),
         e = (p + m * percent$share[5] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[5]**(-1)),
         f = (p + m * percent$share[6] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[6]**(-1)),
         g = (p + m * percent$share[7] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[7]**(-1)),
         h = (p + m * percent$share[8] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[8]**(-1)),
         i = (p + m * percent$share[9] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[9]**(-1)),
         j = (p + m * percent$share[10] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[10]**(-1)),
         sum1 = a+b+c+d+e+f+g+h+i+j,
         a1 = (1 / (m + 2*w*S**(-2)*percent$Kap[1]**(-1))),
         a2 = (1 / (m + 2*w*S**(-2)*percent$Kap[2]**(-1))),
         a3 = (1 / (m + 2*w*S**(-2)*percent$Kap[3]**(-1))),
         a4 = (1 / (m + 2*w*S**(-2)*percent$Kap[4]**(-1))),
         a5 = (1 / (m + 2*w*S**(-2)*percent$Kap[5]**(-1))),
         a6 = (1 / (m + 2*w*S**(-2)*percent$Kap[6]**(-1))),
         a7 = (1 / (m + 2*w*S**(-2)*percent$Kap[7]**(-1))),
         a8 = (1 / (m + 2*w*S**(-2)*percent$Kap[8]**(-1))),
         a9 = (1 / (m + 2*w*S**(-2)*percent$Kap[9]**(-1))),
         a10 = (1 / (m + 2*w*S**(-2)*percent$Kap[10]**(-1))),
         sum2 = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10,
         o = (sum1 - Htotal) / sum2
  )

trade_cost <- trade_cost %>% 
  dplyr::select(time, S, Htotal, o)

# calculate Harvest: 

for (t in 1:101) {
  
  #fisher 1
  n = (1/(1-a))*w*(1/(q*percent$Kap[1]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[1]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H1[t] = harvest$x[1]
  
  # fisher 2
  n = (1/(1-a))*w*(1/(q*percent$Kap[2]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[2]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H2[t] = harvest$x[1]
  
  # fisher 3
  n = (1/(1-a))*w*(1/(q*percent$Kap[3]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[3]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H3[t] = harvest$x[1]
  
  # fisher 4
  n = (1/(1-a))*w*(1/(q*percent$Kap[4]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[4]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H4[t] = harvest$x[1]
  
  # fisher 5
  n = (1/(1-a))*w*(1/(q*percent$Kap[5]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[5]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H5[t] = harvest$x[1]
  
  # fisher 6
  n = (1/(1-a))*w*(1/(q*percent$Kap[6]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[6]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H6[t] = harvest$x[1]
  
  # fisher 7
  n = (1/(1-a))*w*(1/(q*percent$Kap[7]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[7]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H7[t] = harvest$x[1]
  
  # fisher 8
  n = (1/(1-a))*w*(1/(q*percent$Kap[8]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[8]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H8[t] = harvest$x[1]
  
  # fisher 9
  n = (1/(1-a))*w*(1/(q*percent$Kap[9]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[9]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H9[t] = harvest$x[1]
  
  # fisher 10
  n = (1/(1-a))*w*(1/(q*percent$Kap[10]**a*trade_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_cost$o[t] - n*x**b - m*x + m*percent$share[10]*trade_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_cost$H10[t] = harvest$x[1]
  
  
}

# calculate profits


trade_cost <- trade_cost %>% 
  mutate(H_sum = H1 + H2+ H3 + H4 + H5 + H6 +  H7 + H8 + H9 + H10,
         profit1 =
           p*H1 - w * (H1/(q*percent$Kap[1]**a*S))**(1/(1-a))
         - o * (H1 - percent$share[1] * Htotal) 
         - m/2 * (H1 - percent$share[1] * Htotal)**2,
         profit2 =
           p*H2 - w * (H2/(q*percent$Kap[2]**a*S))**(1/(1-a)) 
         - o * (H2 - percent$share[2] * Htotal) 
         - m/2 * (H2 - percent$share[2] * Htotal)**2,
         profit3 =
           p*H3 - w * (H3/(q*percent$Kap[3]**a*S))**(1/(1-a))
         - o * (H3 - percent$share[3] * Htotal) 
         - m/2 * (H3 - percent$share[3] * Htotal)**2,
         profit4 =
           p*H4 - w * (H4/(q*percent$Kap[4]**a*S))**(1/(1-a)) 
         - o * (H4 - percent$share[4] * Htotal) 
         - m/2 * (H4 - percent$share[4] * Htotal)**2,
         profit5 =
           p*H5 - w * (H5/(q*percent$Kap[5]**a*S))**(1/(1-a))
         - o * (H5 - percent$share[5] * Htotal) 
         - m/2 * (H5 - percent$share[5] * Htotal)**2,
         profit6 =
           p*H6 - w * (H6/(q*percent$Kap[6]**a*S))**(1/(1-a))
         - o * (H6 - percent$share[6] * Htotal) 
         - m/2 * (H6 - percent$share[6] * Htotal)**2,
         profit7 =
           p*H7 - w * (H7/(q*percent$Kap[7]**a*S))**(1/(1-a))
         - o * (H7 - percent$share[7] * Htotal) 
         - m/2 * (H7 - percent$share[7] * Htotal)**2,
         profit8 =
           p*H8 - w * (H8/(q*percent$Kap[8]**a*S))**(1/(1-a))
         - o * (H8 - percent$share[8] * Htotal) 
         - m/2 * (H8 - percent$share[8] * Htotal)**2,
         profit9 =
           p*H9 - w * (H9/(q*percent$Kap[9]**a*S))**(1/(1-a))
         - o * (H9 - percent$share[9] * Htotal) 
         - m/2 * (H9 - percent$share[9] * Htotal)**2,
         profit10 =
           p*H10 - w * (H10/(q*percent$Kap[10]**a*S))**(1/(1-a))
         - o * (H10 - percent$share[10] * Htotal) 
         - m/2 * (H10 - percent$share[10] * Htotal)**2,
         mcost1 = 
           m / 2 * (H1-percent$share[1]*Htotal)**2,
         mcost2 = 
           m / 2 * (H2-percent$share[2]*Htotal)**2,
         mcost3 =
           m / 2 * (H3-percent$share[3]*Htotal)**2,
         mcost4 = 
           m / 2 * (H4-percent$share[4]*Htotal)**2,
         mcost5 = 
           m / 2 * (H5-percent$share[5]*Htotal)**2,
         mcost6 = 
           m / 2 * (H6-percent$share[6]*Htotal)**2,
         mcost7 =
           m / 2 * (H7-percent$share[7]*Htotal)**2,
         mcost8 = 
           m / 2 * (H8-percent$share[8]*Htotal)**2,
         mcost9 = 
           m / 2 * (H9-percent$share[9]*Htotal)**2,
         mcost10 = 
           m / 2 * (H10-percent$share[10]*Htotal)**2,
         msum = mcost1 + mcost2 + mcost3 +mcost4 +mcost5+mcost6+mcost7+mcost8+mcost9+mcost10
  )


trade_cost = trade_cost %>% relocate(H_sum, .after=Htotal)

trade_cost = trade_cost %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate in percent of stock 
trade_cost$gS[1] = 1

trade_cost_path <- trade_cost[abs(trade_cost$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state


trade_cost_st_st <- trade_cost[abs(trade_cost$gS) >= 0.001, ]
trade_cost_st_st <- tail(trade_cost_st_st, n=1) 

# steady state profits

profit_cost <- c(trade_cost_st_st$profit1[1], trade_cost_st_st$profit2[1], trade_cost_st_st$profit3[1], trade_cost_st_st$profit4[1], 
                 trade_cost_st_st$profit5[1], trade_cost_st_st$profit6[1], trade_cost_st_st$profit7[1], trade_cost_st_st$profit8[1],
                 trade_cost_st_st$profit9[1], trade_cost_st_st$profit10[1])
profit_cost <- data.frame(profit_cost)
profit <- bind_cols(profit, profit_cost)

profit = profit %>%
  mutate(sum_cost = sum(profit_cost))

# Case 4: No Trade and Equal Allocation #####

tac_eq <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(tac_eq) <- c('time', 'S', 'Htotal')
tac_eq$time <- c(1:102)
tac_eq$S[1] = 1
tac_eq$Htotal[1] = 0

control = r/2 #linear harvest control for MSY

percent$share = 0.1 #allocation now equal 

for (t in 1:101) {
  
  tac_eq$Htotal[t] = control * tac_eq$S[t] 
  
  tac_eq$S[t+1] = tac_eq$S[t] + (r*tac_eq$S[t]*(1-(tac_eq$S[t]/k))) - tac_eq$Htotal[t]
  
}


tac_eq <- tac_eq %>% 
  mutate(H = sapply(1:10, function (i) percent$share[i] * Htotal),
         Eff = sapply(1:10, function (i) (percent$share[i]*r*S/(2*q*S))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         H_sum = H[,1] + H[,2]+ H[,3] + H[,4] + H[,5] + H[,6] +  H[,7] + H[,8] + H[,9] + H[,10],
         Eff_sum = Eff[,1] + Eff[,2]+ Eff[,3] + Eff[,4] + Eff[,5] + Eff[,6] +  Eff[,7] + Eff[,8] + Eff[,9] + Eff[,10],
         profit = sapply(1:10, function (i)
           p * percent$share[i] * Htotal - w * (percent$share[i]*r/(2*q))^(1/(1-a)) * percent$Kap[i]^(a/(a-1))),
         share1 = H[,1]/Htotal)
tac_eq = tac_eq %>% relocate(H_sum, .after=Htotal)
tac_eq <- tac_eq %>% 
  mutate(o = p - (w*Eff_sum**a)/((1-a)*q*Kap_total**a*S))


tac_eq = tac_eq %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate in percent of stock 
tac_eq$gS[1]=1

tac_eq_path <- tac_eq[abs(tac_eq$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state

## steady state profits

tac_eq_st_st <- tac_eq[abs(tac_eq$gS) >= 0.001, ]
tac_eq_st_st <- tail(tac_eq_st_st, n=1) 

profit_tac_eq <- data.frame(sapply(1:10, function (i) profit <- (tac_eq_st_st$profit[,i])))
colnames(profit_tac_eq) <- c("profit_tac_eq") 
profit <- bind_cols(profit, profit_tac_eq)

profit = profit %>%
  mutate(sum_tac_eq = sum(profit_tac_eq))


# Case 5: Full Trade and Equal Allocation ##### 

percent$share = 0.1 # equal allocation

m = 0 # no transaction cost

trade_eq_no_cost <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(trade_eq_no_cost) <- c('time', 'S', 'Htotal')
trade_eq_no_cost$time <- c(1:102)

trade_eq_no_cost$S[1] = 1 #
trade_eq_no_cost$Htotal[1] = 0

# calculate willingness to pay


for (t in 1:101) {
  
  trade_eq_no_cost$Htotal[t] = control * trade_eq_no_cost$S[t] #tac in each t = total harvest (assumption)
  
  trade_eq_no_cost$S[t+1] = trade_eq_no_cost$S[t] + (r*trade_eq_no_cost$S[t]*(1-(trade_eq_no_cost$S[t]/k))) - trade_eq_no_cost$Htotal[t]
  
}

trade_eq_no_cost <- trade_eq_no_cost %>% 
  mutate(a = (p + m * percent$share[1] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[1]**(-1)),
         b = (p + m * percent$share[2] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[2]**(-1)),
         c = (p + m * percent$share[3] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[3]**(-1)),
         d = (p + m * percent$share[4] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[4]**(-1)),
         e = (p + m * percent$share[5] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[5]**(-1)),
         f = (p + m * percent$share[6] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[6]**(-1)),
         g = (p + m * percent$share[7] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[7]**(-1)),
         h = (p + m * percent$share[8] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[8]**(-1)),
         i = (p + m * percent$share[9] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[9]**(-1)),
         j = (p + m * percent$share[10] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[10]**(-1)),
         sum1 = a+b+c+d+e+f+g+h+i+j,
         a1 = (1 / (m + 2*w*S**(-2)*percent$Kap[1]**(-1))),
         a2 = (1 / (m + 2*w*S**(-2)*percent$Kap[2]**(-1))),
         a3 = (1 / (m + 2*w*S**(-2)*percent$Kap[3]**(-1))),
         a4 = (1 / (m + 2*w*S**(-2)*percent$Kap[4]**(-1))),
         a5 = (1 / (m + 2*w*S**(-2)*percent$Kap[5]**(-1))),
         a6 = (1 / (m + 2*w*S**(-2)*percent$Kap[6]**(-1))),
         a7 = (1 / (m + 2*w*S**(-2)*percent$Kap[7]**(-1))),
         a8 = (1 / (m + 2*w*S**(-2)*percent$Kap[8]**(-1))),
         a9 = (1 / (m + 2*w*S**(-2)*percent$Kap[9]**(-1))),
         a10 = (1 / (m + 2*w*S**(-2)*percent$Kap[10]**(-1))),
         sum2 = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10,
         o = (sum1 - Htotal) / sum2
  )

trade_eq_no_cost <- trade_eq_no_cost %>% 
  dplyr::select(time, S, Htotal, o)

# calcultae harvest

for (t in 1:101) {
  
  #fisher 1
  n = (1/(1-a))*w*(1/(q*percent$Kap[1]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[1]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H1[t] = harvest$x[1]
  
  # fisher 2
  n = (1/(1-a))*w*(1/(q*percent$Kap[2]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[2]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H2[t] = harvest$x[1]
  
  # fisher 3
  n = (1/(1-a))*w*(1/(q*percent$Kap[3]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[3]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H3[t] = harvest$x[1]
  
  # fisher 4
  n = (1/(1-a))*w*(1/(q*percent$Kap[4]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[4]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H4[t] = harvest$x[1]
  
  # fisher 5
  n = (1/(1-a))*w*(1/(q*percent$Kap[5]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[5]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H5[t] = harvest$x[1]
  
  # fisher 6
  n = (1/(1-a))*w*(1/(q*percent$Kap[6]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[6]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H6[t] = harvest$x[1]
  
  # fisher 7
  n = (1/(1-a))*w*(1/(q*percent$Kap[7]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[7]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H7[t] = harvest$x[1]
  
  # fisher 8
  n = (1/(1-a))*w*(1/(q*percent$Kap[8]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[8]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H8[t] = harvest$x[1]
  
  # fisher 9
  n = (1/(1-a))*w*(1/(q*percent$Kap[9]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[9]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H9[t] = harvest$x[1]
  
  # fisher 10
  n = (1/(1-a))*w*(1/(q*percent$Kap[10]**a*trade_eq_no_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_no_cost$o[t] - n*x**b - m*x + m*percent$share[10]*trade_eq_no_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_no_cost$H10[t] = harvest$x[1]
  
  
}

# calculate profits


trade_eq_no_cost <- trade_eq_no_cost %>% 
  mutate(H_sum = H1 + H2+ H3 + H4 + H5 + H6 +  H7 + H8 + H9 + H10,
         profit1 =
           p*H1 - w * (H1/(q*percent$Kap[1]**a*S))**(1/(1-a))
         - o * (H1 - percent$share[1] * Htotal) 
         - m/2 * (H1 - percent$share[1] * Htotal)**2,
         profit2 =
           p*H2 - w * (H2/(q*percent$Kap[2]**a*S))**(1/(1-a)) 
         - o * (H2 - percent$share[2] * Htotal) 
         - m/2 * (H2 - percent$share[2] * Htotal)**2,
         profit3 =
           p*H3 - w * (H3/(q*percent$Kap[3]**a*S))**(1/(1-a))
         - o * (H3 - percent$share[3] * Htotal) 
         - m/2 * (H3 - percent$share[3] * Htotal)**2,
         profit4 =
           p*H4 - w * (H4/(q*percent$Kap[4]**a*S))**(1/(1-a)) 
         - o * (H4 - percent$share[4] * Htotal) 
         - m/2 * (H4 - percent$share[4] * Htotal)**2,
         profit5 =
           p*H5 - w * (H5/(q*percent$Kap[5]**a*S))**(1/(1-a))
         - o * (H5 - percent$share[5] * Htotal) 
         - m/2 * (H5 - percent$share[5] * Htotal)**2,
         profit6 =
           p*H6 - w * (H6/(q*percent$Kap[6]**a*S))**(1/(1-a))
         - o * (H6 - percent$share[6] * Htotal) 
         - m/2 * (H6 - percent$share[6] * Htotal)**2,
         profit7 =
           p*H7 - w * (H7/(q*percent$Kap[7]**a*S))**(1/(1-a))
         - o * (H7 - percent$share[7] * Htotal) 
         - m/2 * (H7 - percent$share[7] * Htotal)**2,
         profit8 =
           p*H8 - w * (H8/(q*percent$Kap[8]**a*S))**(1/(1-a))
         - o * (H8 - percent$share[8] * Htotal) 
         - m/2 * (H8 - percent$share[8] * Htotal)**2,
         profit9 =
           p*H9 - w * (H9/(q*percent$Kap[9]**a*S))**(1/(1-a))
         - o * (H9 - percent$share[9] * Htotal) 
         - m/2 * (H9 - percent$share[9] * Htotal)**2,
         profit10 =
           p*H10 - w * (H10/(q*percent$Kap[10]**a*S))**(1/(1-a))
         - o * (H10 - percent$share[10] * Htotal) 
         - m/2 * (H10 - percent$share[10] * Htotal)**2,
         mcost1 = 
           m / 2 * (H1-percent$share[1]*Htotal)**2,
         mcost2 = 
           m / 2 * (H2-percent$share[2]*Htotal)**2,
         mcost3 =
           m / 2 * (H3-percent$share[3]*Htotal)**2,
         mcost4 = 
           m / 2 * (H4-percent$share[4]*Htotal)**2,
         mcost5 = 
           m / 2 * (H5-percent$share[5]*Htotal)**2,
         mcost6 = 
           m / 2 * (H6-percent$share[6]*Htotal)**2,
         mcost7 =
           m / 2 * (H7-percent$share[7]*Htotal)**2,
         mcost8 = 
           m / 2 * (H8-percent$share[8]*Htotal)**2,
         mcost9 = 
           m / 2 * (H9-percent$share[9]*Htotal)**2,
         mcost10 = 
           m / 2 * (H10-percent$share[10]*Htotal)**2,
         msum = mcost1 + mcost2 + mcost3 +mcost4 +mcost5+mcost6+mcost7+mcost8+mcost9+mcost10
  )


trade_eq_no_cost = trade_eq_no_cost %>% relocate(H_sum, .after=Htotal)

trade_eq_no_cost = trade_eq_no_cost %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate in percent of stock 
trade_eq_no_cost$gS[1] = 1

trade_eq_no_cost_path <- trade_eq_no_cost[abs(trade_eq_no_cost$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state


trade_eq_no_cost_st_st <- trade_eq_no_cost[abs(trade_eq_no_cost$gS) >= 0.001, ]
trade_eq_no_cost_st_st <- tail(trade_eq_no_cost_st_st, n=1) 

# steady state profits

profit_eq_no_cost <- c(trade_eq_no_cost_st_st$profit1[1], trade_eq_no_cost_st_st$profit2[1], trade_eq_no_cost_st_st$profit3[1], trade_eq_no_cost_st_st$profit4[1], 
                       trade_eq_no_cost_st_st$profit5[1], trade_eq_no_cost_st_st$profit6[1], trade_eq_no_cost_st_st$profit7[1], trade_eq_no_cost_st_st$profit8[1],
                       trade_eq_no_cost_st_st$profit9[1], trade_eq_no_cost_st_st$profit10[1])
profit_eq_no_cost <- data.frame(profit_eq_no_cost)
profit <- bind_cols(profit, profit_eq_no_cost)

profit = profit %>%
  mutate(sum_eq_no_cost = sum(profit_eq_no_cost))


# Case 6: Limited Trade and Equal Allocation #####

percent$share = 0.1 # equal allocation

m = 0.01 # transaction cost parameter

trade_eq_cost <- data.frame(matrix(NA, nrow = 102, ncol = 3))
names(trade_eq_cost) <- c('time', 'S', 'Htotal')
trade_eq_cost$time <- c(1:102)

trade_eq_cost$S[1] = 1 
trade_eq_cost$Htotal[1] = 0

# calculate willingness to pay


for (t in 1:101) {
  
  trade_eq_cost$Htotal[t] = control * trade_eq_cost$S[t] 
  
  trade_eq_cost$S[t+1] = trade_eq_cost$S[t] + (r*trade_eq_cost$S[t]*(1-(trade_eq_cost$S[t]/k))) - trade_eq_cost$Htotal[t]
  
}

trade_eq_cost <- trade_eq_cost %>% 
  mutate(a = (p + m * percent$share[1] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[1]**(-1)),
         b = (p + m * percent$share[2] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[2]**(-1)),
         c = (p + m * percent$share[3] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[3]**(-1)),
         d = (p + m * percent$share[4] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[4]**(-1)),
         e = (p + m * percent$share[5] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[5]**(-1)),
         f = (p + m * percent$share[6] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[6]**(-1)),
         g = (p + m * percent$share[7] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[7]**(-1)),
         h = (p + m * percent$share[8] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[8]**(-1)),
         i = (p + m * percent$share[9] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[9]**(-1)),
         j = (p + m * percent$share[10] * Htotal) / (m + 2 * w * S**(-2) * percent$Kap[10]**(-1)),
         sum1 = a+b+c+d+e+f+g+h+i+j,
         a1 = (1 / (m + 2*w*S**(-2)*percent$Kap[1]**(-1))),
         a2 = (1 / (m + 2*w*S**(-2)*percent$Kap[2]**(-1))),
         a3 = (1 / (m + 2*w*S**(-2)*percent$Kap[3]**(-1))),
         a4 = (1 / (m + 2*w*S**(-2)*percent$Kap[4]**(-1))),
         a5 = (1 / (m + 2*w*S**(-2)*percent$Kap[5]**(-1))),
         a6 = (1 / (m + 2*w*S**(-2)*percent$Kap[6]**(-1))),
         a7 = (1 / (m + 2*w*S**(-2)*percent$Kap[7]**(-1))),
         a8 = (1 / (m + 2*w*S**(-2)*percent$Kap[8]**(-1))),
         a9 = (1 / (m + 2*w*S**(-2)*percent$Kap[9]**(-1))),
         a10 = (1 / (m + 2*w*S**(-2)*percent$Kap[10]**(-1))),
         sum2 = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10,
         o = (sum1 - Htotal) / sum2
  )

trade_eq_cost <- trade_eq_cost %>% 
  dplyr::select(time, S, Htotal, o)

# calculate harvest

for (t in 1:101) {
  
  #fisher 1
  n = (1/(1-a))*w*(1/(q*percent$Kap[1]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[1]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H1[t] = harvest$x[1]
  
  # fisher 2
  n = (1/(1-a))*w*(1/(q*percent$Kap[2]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[2]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H2[t] = harvest$x[1]
  
  # fisher 3
  n = (1/(1-a))*w*(1/(q*percent$Kap[3]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[3]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H3[t] = harvest$x[1]
  
  # fisher 4
  n = (1/(1-a))*w*(1/(q*percent$Kap[4]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[4]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H4[t] = harvest$x[1]
  
  # fisher 5
  n = (1/(1-a))*w*(1/(q*percent$Kap[5]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[5]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H5[t] = harvest$x[1]
  
  # fisher 6
  n = (1/(1-a))*w*(1/(q*percent$Kap[6]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[6]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H6[t] = harvest$x[1]
  
  # fisher 7
  n = (1/(1-a))*w*(1/(q*percent$Kap[7]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[7]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H7[t] = harvest$x[1]
  
  # fisher 8
  n = (1/(1-a))*w*(1/(q*percent$Kap[8]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[8]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H8[t] = harvest$x[1]
  
  # fisher 9
  n = (1/(1-a))*w*(1/(q*percent$Kap[9]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[9]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H9[t] = harvest$x[1]
  
  # fisher 10
  n = (1/(1-a))*w*(1/(q*percent$Kap[10]**a*trade_eq_cost$S[t]))**(1/(1-a))
  b = (a/(1-a))
  
  fn <- function (x) {p - trade_eq_cost$o[t] - n*x**b - m*x + m*percent$share[10]*trade_eq_cost$Htotal[t]}
  
  harvest <- nleqslv(1,fn)
  
  trade_eq_cost$H10[t] = harvest$x[1]
  
  
}

# calculate profits


trade_eq_cost <- trade_eq_cost %>% 
  mutate(H_sum = H1 + H2+ H3 + H4 + H5 + H6 +  H7 + H8 + H9 + H10,
         profit1 =
           p*H1 - w * (H1/(q*percent$Kap[1]**a*S))**(1/(1-a))
         - o * (H1 - percent$share[1] * Htotal) 
         - m/2 * (H1 - percent$share[1] * Htotal)**2,
         profit2 =
           p*H2 - w * (H2/(q*percent$Kap[2]**a*S))**(1/(1-a)) 
         - o * (H2 - percent$share[2] * Htotal) 
         - m/2 * (H2 - percent$share[2] * Htotal)**2,
         profit3 =
           p*H3 - w * (H3/(q*percent$Kap[3]**a*S))**(1/(1-a))
         - o * (H3 - percent$share[3] * Htotal) 
         - m/2 * (H3 - percent$share[3] * Htotal)**2,
         profit4 =
           p*H4 - w * (H4/(q*percent$Kap[4]**a*S))**(1/(1-a)) 
         - o * (H4 - percent$share[4] * Htotal) 
         - m/2 * (H4 - percent$share[4] * Htotal)**2,
         profit5 =
           p*H5 - w * (H5/(q*percent$Kap[5]**a*S))**(1/(1-a))
         - o * (H5 - percent$share[5] * Htotal) 
         - m/2 * (H5 - percent$share[5] * Htotal)**2,
         profit6 =
           p*H6 - w * (H6/(q*percent$Kap[6]**a*S))**(1/(1-a))
         - o * (H6 - percent$share[6] * Htotal) 
         - m/2 * (H6 - percent$share[6] * Htotal)**2,
         profit7 =
           p*H7 - w * (H7/(q*percent$Kap[7]**a*S))**(1/(1-a))
         - o * (H7 - percent$share[7] * Htotal) 
         - m/2 * (H7 - percent$share[7] * Htotal)**2,
         profit8 =
           p*H8 - w * (H8/(q*percent$Kap[8]**a*S))**(1/(1-a))
         - o * (H8 - percent$share[8] * Htotal) 
         - m/2 * (H8 - percent$share[8] * Htotal)**2,
         profit9 =
           p*H9 - w * (H9/(q*percent$Kap[9]**a*S))**(1/(1-a))
         - o * (H9 - percent$share[9] * Htotal) 
         - m/2 * (H9 - percent$share[9] * Htotal)**2,
         profit10 =
           p*H10 - w * (H10/(q*percent$Kap[10]**a*S))**(1/(1-a))
         - o * (H10 - percent$share[10] * Htotal) 
         - m/2 * (H10 - percent$share[10] * Htotal)**2,
         mcost1 = 
           m / 2 * (H1-percent$share[1]*Htotal)**2,
         mcost2 = 
           m / 2 * (H2-percent$share[2]*Htotal)**2,
         mcost3 =
           m / 2 * (H3-percent$share[3]*Htotal)**2,
         mcost4 = 
           m / 2 * (H4-percent$share[4]*Htotal)**2,
         mcost5 = 
           m / 2 * (H5-percent$share[5]*Htotal)**2,
         mcost6 = 
           m / 2 * (H6-percent$share[6]*Htotal)**2,
         mcost7 =
           m / 2 * (H7-percent$share[7]*Htotal)**2,
         mcost8 = 
           m / 2 * (H8-percent$share[8]*Htotal)**2,
         mcost9 = 
           m / 2 * (H9-percent$share[9]*Htotal)**2,
         mcost10 = 
           m / 2 * (H10-percent$share[10]*Htotal)**2,
         msum = mcost1 + mcost2 + mcost3 +mcost4 +mcost5+mcost6+mcost7+mcost8+mcost9+mcost10
  )


trade_eq_cost = trade_eq_cost %>% relocate(H_sum, .after=Htotal)

trade_eq_cost = trade_eq_cost %>%
  mutate(gS = ((S - lag(S))/S * 100)) # calculate growth rate in percent of stock 
trade_eq_cost$gS[1] = 1

trade_eq_cost_path <- trade_eq_cost[abs(trade_eq_cost$gS) >= 0.001, ] # delete when growth rate of stock is smaller then 1% -> nearly steady state


trade_eq_cost_st_st <- trade_eq_cost[abs(trade_eq_cost$gS) >= 0.001, ]
trade_eq_cost_st_st <- tail(trade_eq_cost_st_st, n=1) 

# steady state profits

profit_eq_cost <- c(trade_eq_cost_st_st$profit1[1], trade_eq_cost_st_st$profit2[1], trade_eq_cost_st_st$profit3[1], trade_eq_cost_st_st$profit4[1], 
                    trade_eq_cost_st_st$profit5[1], trade_eq_cost_st_st$profit6[1], trade_eq_cost_st_st$profit7[1], trade_eq_cost_st_st$profit8[1],
                    trade_eq_cost_st_st$profit9[1], trade_eq_cost_st_st$profit10[1])
profit_eq_cost <- data.frame(profit_eq_cost)
profit <- bind_cols(profit, profit_eq_cost)

profit = profit %>%
  mutate(sum_eq_cost = sum(profit_eq_cost))

# Create Equity Values Based on (inverse) Simpson-Index ##### 

profit_comp = profit %>%
  mutate(ratio_oa = profit_oa/sum_oa, .after = sum_oa,
         ratio_tac = profit_tac/sum_tac,
         ratio_no_cost = profit_no_cost/sum_no_cost,
         ratio_cost = profit_cost/sum_cost,
  )

profit_comp = profit_comp %>% relocate(ratio_tac, .after=sum_tac)
profit_comp = profit_comp %>% relocate(ratio_no_cost, .after=sum_no_cost)
profit_comp = profit_comp %>% relocate(ratio_cost, .after=sum_cost)

profit_comp = profit_comp %>%
  mutate(in_no_cost_trade = (ratio_no_cost - ratio_tac) / ratio_tac * 100,
         in_cost_nocost = (ratio_cost - ratio_no_cost) / ratio_cost * 100
  )

equity_tac = c()
equity_tac = profit %>%
  mutate(p = profit_tac / sum_tac,
         p_2 = p**2,
         sum_p2 = sum(p_2),
         D_tac = 1/sum_p2)
equity_tac <- subset(equity_tac, select= c("D_tac"))

equity_tac_eq = c()
equity_tac_eq = profit %>%
  mutate(p = profit_tac_eq / sum_tac_eq,
         p_2 = p**2,
         sum_p2 = sum(p_2),
         D_tac_eq = 1/sum_p2)
equity_tac_eq <- subset(equity_tac_eq, select= c("D_tac_eq"))

equity_cost = c()
equity_cost = profit %>%
  mutate(p = profit_eq_cost / sum_eq_cost,
         p_2 = p**2,
         sum_p2 = sum(p_2),
         D_eq_cost = 1/sum_p2)
equity_cost <- subset(equity_cost, select= c("D_eq_cost"))

equity_no_cost = c()
equity_no_cost = profit %>%
  mutate(p = profit_eq_no_cost / sum_eq_no_cost,
         p_2 = p**2,
         sum_p2 = sum(p_2),
         D_eq_no_cost = 1/sum_p2)
equity_no_cost <- subset(equity_no_cost, select= c("D_eq_no_cost"))


equity <- bind_cols(equity_tac, equity_tac_eq, equity_cost, equity_no_cost)

# Equity Plot#####

smooth_equity <- ggplot()+
  scale_x_continuous(expand = c(0,0), breaks=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), limits = c(0.9,10.1) ) +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_tac ,color="a"), se=F, size=1) +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_tac_eq, color="b"), se=F, size=1) +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_eq_cost, color="c"), se=F, size=1) +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_eq_no_cost, color="d"), se=F, size=1) +
  ylab("Profit in million Euro") +
  xlab("Fisher") +
  scale_color_manual(values=c("darkseagreen","slategray", "slategray3", "slategray1"),
                     name="Transferability",
                     label=c("All - skill","INTQs - equal", "ITQs cost - equal", "ITQs no cost - equal"))
smooth_equity

smooth_equity <- smooth_equity + theme_bw()
smooth_equity

smooth_equity <- smooth_equity + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                       legend.position = c(0.85, 0.85),
                                       axis.text=element_text(size=14),
                                       axis.title=element_text(size=15,face="bold"),
                                       legend.title = element_text(size=15, face="bold"),
                                       legend.text = element_text(size=15))
smooth_equity
ggsave("smooth_equity.eps", dpi = "print", width = 9, height = 6)

# Comparison Profit Plot Skill Allocation#####

bar_tac <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z", y=profit_tac, fill="darkseagreen"), color="black", stat="identity") +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_tac), se=F, color="black", size=0.5) +
  ylab("Individual profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("darkseagreen"),
                    name="Transferability",
                    label=c("INTQs, ITQs no cost, ITQs cost"))+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10"))
bar_tac

bar_tac <- bar_tac + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                           legend.position = c(0.65, 0.91),
                           panel.background=element_blank(),
                           axis.text=element_text(size=11),
                           axis.title=element_text(size=12,face="bold"),
                           legend.title = element_text(size=10, face="bold"),
                           legend.text = element_text(size=10))
bar_tac

#ggsave("bar_cap.eps", dpi = "print", width = 7, height = 5)

bar_sum_profit_steady <- ggplot() + 
  geom_bar(data=profit[1, ], aes(x="1. no trade", y=sum_tac), color="black", fill="darkseagreen4", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="2. trade with cost", y=sum_cost), color="black", fill="darkseagreen2", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="3. trade without cost", y=sum_no_cost), color="black", fill="honeydew2", stat="identity") +
  ylab("Total profit in million Euro") +
  xlab("Transferability") +
  scale_x_discrete(labels=c("INTQs", "ITQs no cost", "ITQs cost"))
bar_sum_profit_steady

bar_sum_profit_steady <- bar_sum_profit_steady + theme(panel.background=element_blank(), axis.text=element_text(size=11),
                                                       axis.title=element_text(size=12,face="bold"),
                                                       axis.title.x = element_text(vjust=-1))
bar_sum_profit_steady

# Comparison Profit Plot Equal Allocation #####


bar_tac_eq2 <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="1 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  ylab("Individual profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("darkseagreen","slategray"),
                    name="Transferability",
                    label=c("All - skill","INTQs"))+
  scale_x_discrete(labels=c("1","1","2", "2","3", "3","4", "4","5", "5",
                            "6", "6","7", "7","8", "8","9", "9","10", "10"))
bar_tac_eq2

bar_tac_eq2 <- bar_tac_eq2 + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                   legend.position = c(0.85, 0.91),
                                   panel.background=element_blank(),
                                   axis.text=element_text(size=11),
                                   axis.title=element_text(size=12,face="bold"),
                                   legend.title = element_text(size=10, face="bold"),
                                   legend.text = element_text(size=10))
bar_tac_eq2


bar_tac_eq <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (a)", y=profit_tac_eq, fill="slategray"), color="black", stat="identity") +
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_tac_eq), se=F, color="black", size=0.5) +
  ylab("Profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("slategray"),
                    name="Transferability",
                    label=c("INTQ"))+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10"))
bar_tac_eq

bar_tac_eq <- bar_tac_eq + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                 legend.position = c(0.85, 0.91),
                                 panel.background=element_blank(),
                                 axis.text=element_text(size=11),
                                 axis.title=element_text(size=12,face="bold"),
                                 legend.title = element_text(size=10, face="bold"),
                                 legend.text = element_text(size=10))
bar_tac_eq



bar_eq_cost <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") + 
  geom_bar(data=profit[6, ], aes(x="6 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (b)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") + 
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_eq_cost), se=F, color="black", size=0.5) +
  ylab("Profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("slategray3"),
                    name="Transferability",
                    label=c("ITQ with cost"))+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10"))
bar_eq_cost

bar_eq_cost <- bar_eq_cost + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                   legend.position = c(0.85, 0.91),
                                   panel.background=element_blank(),
                                   axis.text=element_text(size=11),
                                   axis.title=element_text(size=12,face="bold"),
                                   legend.title = element_text(size=10, face="bold"),
                                   legend.text = element_text(size=10))
bar_eq_cost

bar_eq_cost2 <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="1 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (a)", y=profit_eq_cost, fill="slategray3"), color="black", stat="identity") +
  ylab("Individual profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("darkseagreen","slategray3"),
                    name="Transferability",
                    label=c("All - skill","ITQs with cost"))+
  scale_x_discrete(labels=c("1","1","2", "2","3", "3","4", "4","5", "5",
                            "6", "6","7", "7","8", "8","9", "9","10", "10"))
bar_eq_cost2

bar_eq_cost2 <- bar_eq_cost2 + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                     legend.position = c(0.85, 0.91),
                                     panel.background=element_blank(),
                                     axis.text=element_text(size=11),
                                     axis.title=element_text(size=12,face="bold"),
                                     legend.title = element_text(size=10, face="bold"),
                                     legend.text = element_text(size=10))
bar_eq_cost2


bar_eq_no_cost <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") + 
  geom_bar(data=profit[6, ], aes(x="6 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (c)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") + 
  geom_smooth(data=profit, aes(x=seq(1:10), y=profit_eq_no_cost), se=F, color="black", size=0.5) +
  ylab("Profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("slategray1"),
                    name="Transferability",
                    label=c("ITQ no cost"))+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10"))
bar_eq_no_cost

bar_eq_no_cost <- bar_eq_no_cost + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                         legend.position = c(0.85, 0.91),
                                         panel.background=element_blank(),
                                         axis.text=element_text(size=11),
                                         axis.title=element_text(size=12,face="bold"),
                                         legend.title = element_text(size=10, face="bold"),
                                         legend.text = element_text(size=10))
bar_eq_no_cost


bar_eq_no_cost2 <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="1 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (a)", y=profit_eq_no_cost, fill="slategray1"), color="black", stat="identity") +
  ylab("Individual profit in million Euro") +
  xlab("Fisher") +
  scale_fill_manual(values=c("darkseagreen","slategray1"),
                    name="Transferability",
                    label=c("All - skill","ITQs no cost"))+
  scale_x_discrete(labels=c("1","1","2", "2","3", "3","4", "4","5", "5",
                            "6", "6","7", "7","8", "8","9", "9","10", "10"))
bar_eq_no_cost2

bar_eq_no_cost2 <- bar_eq_no_cost2 + theme(axis.title.x = element_text(vjust=-1, hjust=0.5),
                                           legend.position = c(0.85, 0.91),
                                           panel.background=element_blank(),
                                           axis.text=element_text(size=11),
                                           axis.title=element_text(size=12,face="bold"),
                                           legend.title = element_text(size=10, face="bold"),
                                           legend.text = element_text(size=10))
bar_eq_no_cost2


bar_sum_profit_steady_equal <- ggplot() + 
  geom_bar(data=profit[1, ], aes(x="a", y=sum_tac_eq), fill="slategray", color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="b", y=sum_eq_cost), fill="slategray3", color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="c", y=sum_eq_no_cost), fill="slategray1", color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="d", y=sum_no_cost), fill="darkseagreen", color="black", stat="identity") +
  ylab("Total profit in million Euro") +
  xlab("Transferability") +
  scale_x_discrete(labels=c(str_wrap("INTQs equal",width = 5), str_wrap("ITQs(cost) equal", width = 5), 
                            str_wrap("ITQs(nocost) equal", width = 5), str_wrap("All skill", width = 5)))
bar_sum_profit_steady_equal

bar_sum_profit_steady_equal <- bar_sum_profit_steady_equal + theme(panel.background=element_blank(),
                                                                   axis.text=element_text(size=11),
                                                                   axis.title=element_text(size=12,face="bold"),
                                                                   axis.title.x = element_text(vjust=-1))
bar_sum_profit_steady_equal


# Comparison Profit Plot Skill and Equal Allocations #####

bar_steady_profit_equal2 <- ggplot() +
  geom_bar(data=profit[1, ], aes(x="1 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (0)", y=profit_tac, fill="a"), color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="1 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[6, ], aes(x="6 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (a)", y=profit_tac_eq, fill="b"), color="black", stat="identity") +
  geom_bar(data=profit[1, ], aes(x="1 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") + 
  geom_bar(data=profit[6, ], aes(x="6 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (b)", y=profit_eq_cost, fill="c"), color="black", stat="identity") + 
  geom_bar(data=profit[1, ], aes(x="1 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[2, ], aes(x="2 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[3, ], aes(x="3 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[4, ], aes(x="4 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[5, ], aes(x="5 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") + 
  geom_bar(data=profit[6, ], aes(x="6 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[7, ], aes(x="7 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[8, ], aes(x="8 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[9, ], aes(x="9 (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") +
  geom_bar(data=profit[10, ], aes(x="9z (c)", y=profit_eq_no_cost, fill="d"), color="black", stat="identity") + 
  ylab("profit in million USD") +
  xlab("fisher") +
  scale_fill_manual(values=c("darkseagreen", "slategray",  "slategray3", "slategray1"),
                    name="Transferability",
                    label=c("comparison to equal", "no trade (a)", "cost (b)", "no cost (c)"))
bar_steady_profit_equal2

bar_steady_profit_equal2 <- bar_steady_profit_equal2 + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                             axis.title.x = element_text(vjust=+2, hjust=0.5),
                                                             legend.position = c(0.8, 0.8),
                                                             panel.background=element_blank())
bar_steady_profit_equal2


# Harvest Plots #####

harvest_tac <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=tac_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=tac_path, aes(x=S, y=H[,1], color="b")) +
  geom_line(data=tac_path, aes(x=S, y=H[,1], color="b2", size="b")) +
  geom_line(data=tac_path, aes(x=S, y=H[,2], color="c", size="c")) +
  geom_line(data=tac_path, aes(x=S, y=H[,3], color="d", size="d")) +
  geom_line(data=tac_path, aes(x=S, y=H[,4], color="e")) +
  geom_line(data=tac_path, aes(x=S, y=H[,5], color="f")) +
  geom_line(data=tac_path, aes(x=S, y=H[,6], color="g")) +
  geom_line(data=tac_path, aes(x=S, y=H[,7], color="h")) +
  geom_line(data=tac_path, aes(x=S, y=H[,10], color="k")) +
  ylab("Growth/ Harvest (thousand metric tons)") +
  xlab("Biomass (thousand tons)")+
  title("INTQ & ITQ - skill") +
  scale_color_manual(values=c("red","yellow", "blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "black"),
                                   name="Legend", 
                     label=c("Total", "Fisher 1-10", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_tac

harvest_tac <- harvest_tac + theme_bw()

harvest_tac <- harvest_tac + theme(axis.text=element_text(size=12),
                                   axis.title=element_text(size=13,face="bold"),
                                   legend.key.size = unit(1.5,"line"), 
                                   legend.title = element_text(size=14,face="bold"),
                                   legend.text = element_text(size=14),
                                   legend.position = "bottom")
harvest_tac
harvest_tac <- harvest_tac + guides(color=guide_legend(nrow=2))
harvest_tac

harvest_tac_open <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=tac_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=tac_path, aes(x=S, y=H[,1], color="b", size="b")) +
  geom_line(data=tac_path, aes(x=S, y=H[,2], color="c", size="c")) +
  geom_line(data=tac_path, aes(x=S, y=H[,3], color="d", size="d")) +
  geom_line(data=tac_path, aes(x=S, y=H[,4], color="e")) +
  geom_line(data=tac_path, aes(x=S, y=H[,5], color="f")) +
  geom_line(data=tac_path, aes(x=S, y=H[,6], color="g")) +
  geom_line(data=tac_path, aes(x=S, y=H[,7], color="h")) +
  geom_line(data=tac_path, aes(x=S, y=H[,10], color="k")) +
  geom_line(data=oa_path, aes(x=S, y=Htotal, color="x", size="x"), linetype = "dashed") +
  ylab("Growth/ Harvest (thousand metric tons)") +
  xlab("Biomass (thousand tons)")+
  title("INTQs & ITQs - skill") +
  scale_color_manual(values=c("red", "blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "green2", "black"),
                                   name="Legend", 
                     label=c("MSY/ Total", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Open Access", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_tac_open

harvest_tac_open <- harvest_tac_open + theme_bw()

harvest_tac_open <- harvest_tac_open + theme(axis.text=element_text(size=12),
                                             axis.title=element_text(size=13,face="bold"),
                                             legend.key.size = unit(1.5,"line"), 
                                             legend.title = element_text(size=14,face="bold"),
                                             legend.text = element_text(size=14))
harvest_tac_open
ggsave("harvest_tac_open.eps", dpi = "print", width = 7, height = 5)


harvest_cost_all <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=trade_no_cost_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H1, color="b", size="b")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H2, color="c", size="c")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H3, color="d", size="d")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H4, color="e")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H5, color="f")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H6, color="g")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H7, color="h")) +
  geom_line(data=trade_no_cost_path, aes(x=S, y=H10, color="k")) +
  ylab("Growth/ Harvest (thousand tons)") +
  xlab("Biomass (thousand tons)")+
  scale_color_manual(values=c("red","blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "black"),
                                   name="Harvest", 
                     label=c("Total", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_cost_all

harvest_cost_all <- harvest_cost_all + theme_bw()

harvest_cost_all + theme(axis.text=element_text(size=12),
                         axis.title=element_text(size=13,face="bold"),
                         legend.key.size = unit(1.5,"line"), 
                         legend.title = element_text(size=12,face="bold"),
                         legend.text = element_text(size=10))

harvest_tac_trade_no_cost <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=trade_cost_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H1, color="b", size="b")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H2, color="c", size="c")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H3, color="d", size="d")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H4, color="e")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H5, color="f")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H6, color="g")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H7, color="h")) +
  geom_line(data=trade_cost_path, aes(x=S, y=H10, color="k")) +
  ylab("Growth/ Harvest (thousand tons)") +
  xlab("Biomass (thousand tons)")+
  scale_color_manual(values=c("red","blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "black"),
                                   name="Harvest", 
                     label=c("Total", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_tac_trade_no_cost

harvest_tac_trade_no_cost <- harvest_tac_trade_no_cost + theme_bw()

harvest_tac_trade_no_cost + theme(axis.text=element_text(size=12),
                                  axis.title=element_text(size=13,face="bold"),
                                  legend.key.size = unit(1.5,"line"), 
                                  legend.title = element_text(size=12,face="bold"),
                                  legend.text = element_text(size=10))

harvest_tac_eq <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=tac_eq_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,1], color="b", size="b")) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,2])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,3])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,4])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,5])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,6])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,7])) +
  geom_line(data=tac_eq_path, aes(x=S, y=H[,10])) +
  ylab("Growth/ Harvest (thousand metric tons)") +
  xlab("Biomass (thousand tons)")+
  title("INTQ - equal") +
  scale_color_manual(values=c("red","yellow", "black"),
                     name="Harvest", 
                     label=c("Total", "All Fishers", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_tac_eq

harvest_tac_eq <- harvest_tac_eq + theme_bw()

harvest_tac_eq <- harvest_tac_eq + theme(axis.text=element_text(size=12),
                                         axis.title=element_text(size=13,face="bold"),
                                         legend.key.size = unit(1.5,"line"), 
                                         legend.title = element_text(size=12,face="bold"),
                                         legend.text = element_text(size=10))

harvest_eq_no_cost <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H1, color="b", size="b")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H2, color="c", size="c")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H3, color="d", size="d")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H4, color="e")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H5, color="f")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H6, color="g")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H7, color="h")) +
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=H10, color="k")) +
  ylab("Growth/ Harvest (thousand tons)") +
  xlab("Biomass (thousand tons)")+
  title("ITQ - no cost - equal") +
  scale_color_manual(values=c("red","blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "black"),
                                   name="Harvest", 
                     label=c("Total", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_eq_no_cost

harvest_eq_no_cost <- harvest_eq_no_cost + theme_bw()

harvest_eq_no_cost <- harvest_eq_no_cost + theme(axis.text=element_text(size=12),
                                                 axis.title=element_text(size=13,face="bold"),
                                                 legend.key.size = unit(1.5,"line"), 
                                                 legend.title = element_text(size=12,face="bold"),
                                                 legend.text = element_text(size=10))

harvest_eq_cost <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0), breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=trade_eq_cost_path, aes(x=S, y=Htotal, color="a", size="a")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H1, color="b", size="b")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H2, color="c", size="c")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H3, color="d", size="d")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H4, color="e")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H5, color="f")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H6, color="g")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H7, color="h")) +
  geom_line(data=trade_eq_cost_path, aes(x=S, y=H10, color="k")) +
  ylab("Growth/ Harvest (thousand metric tons)") +
  xlab("Biomass (thousand tons)")+
  title("ITQ - cost - skill") +
  scale_color_manual(values=c("red","blue","blueviolet", "burlywood", "darkgray", 
                                   "darkgoldenrod", "coral", "cyan4", "chartreuse4", "black"),
                                   name="Harvest", 
                     label=c("Total", "Fisher 1", "Fisher 2", "Fisher 3", "Fisher 4", "Fisher 5"
                             , "Fisher 6", "Fisher 7 & 8", "Fisher 9 & 10", "Growth Function")) +
  scale_size_manual(values = c(1.1, 1, 1, 1),guide = "none") +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_eq_cost

harvest_eq_cost <- harvest_eq_cost + theme_bw()

harvest_eq_cost <- harvest_eq_cost + theme(axis.text=element_text(size=12),
                                           axis.title=element_text(size=13,face="bold"),
                                           legend.key.size = unit(1.5,"line"), 
                                           legend.title = element_text(size=12,face="bold"),
                                           legend.text = element_text(size=10))

harvest_total <- ggplot () +
  geom_function(fun=function(x) r*x*(1-(x/k)), aes(color="z")) +
  scale_x_continuous( expand = c(0,0) , breaks=c(0, 142, 300, 450, 600, 750, 900, 1050), limits = c(0,1150) ) +
  scale_y_continuous( expand = c(0,0), limits = c(0,250) ) + 
  geom_line(data=trade_eq_no_cost_path, aes(x=S, y=Htotal, color="a")) +
  geom_line(data=oa_path, aes(x=S, y=Htotal, color="b")) +
  ylab("Growth/ Harvest (thousand tons)") +
  xlab("Biomass (thousand tons)")+
  scale_color_manual(values=c("red", "blue","black"),
                     name="Legend", 
                     label=c("Harvest MSY", "Harvest Open Access", "Growth Function")) +
  geom_vline(xintercept = 142, linetype = "dashed", color = "black") 
harvest_total

harvest_total <- harvest_total + theme_bw()

harvest_total + theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=13,face="bold"),
                      legend.key.size = unit(1.5,"line"), 
                      legend.title = element_text(size=12,face="bold"),
                      legend.text = element_text(size=10))


#ggsave("harvest_total.eps", dpi = "print", width = 7, height = 5)

# Put plots together ######

# comparison eq harvest
harvest_comp_eq <- ggarrange(plotlist= list(harvest_tac_eq, harvest_eq_cost, harvest_tac), labels = c("A", "B", "C", "D"), label.x = 0.5, label.y = 0.99, legend = "none", common.legend = TRUE)
harvest_comp_eq #uses the legend from the first enteres plot (harvest_tac)

#ggsave("harvest_comp_eq_work.eps", dpi = "print", width = 10, height = 9)

#bar equal without skill
bartac <- ggarrange(plotlist = list(bar_tac, bar_sum_profit_steady), labels=c("A", "B"))
bartac

#ggsave("bar_tac.eps", dpi = "print", width = 8, height = 4)

# bar equal compared with skill
bareq2 <- ggarrange(plotlist= list(bar_tac_eq2, bar_eq_cost2, bar_eq_no_cost2, bar_sum_profit_steady_equal), 
                    labels = c("A", "B", "C", "D"))
bareq2

#ggsave("bar_eq2.eps", dpi = "print", width = 10, height = 9)



