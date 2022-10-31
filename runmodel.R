library(stringr)
library(data.table)
library(dplyr)
library(brms)
library(betareg)

#setwd('C:/Users/stany/Dropbox/polyu2/research/ovitrap')
dat = fread('mergedat.csv')
dat$oindex2 = (dat$oindex * (dat$nodevice-1)+0.5)/dat$nodevice
dat$monthf = as.factor(dat$month)


#fix(dat)
unique(dat$Year)
length(unique(dat$Year))
unique(dat$month)
length(unique(dat$month))
unique(dat$areano)
length(unique(dat$areano))
ldf = data.frame(Year=rep(2010:2021,each=64*12),month=rep(rep(1:12,each=64),12),areano = rep(1:64,12*12))
ldf  = ldf %>% left_join(dat)

datfit = ldf
datfit$obs = as.numeric((!is.na(datfit$device)) & (!is.na(datfit$Templag1)) & 
                          (!is.na(datfit$Templag2)) & (!is.na(datfit$Templag3)) & 
                          (!is.na(datfit$Templag4)) & (!is.na(datfit$nino4lag1)) & 
                          (!is.na(datfit$Rainfalllag2)) & (!is.na(datfit$Rainfalllag3)) & 
                          (!is.na(datfit$Rainfalllag4)) & (!is.na(datfit$Rainfalllag1)) &
                          (!is.na(datfit$Ozonelag2)) & (!is.na(datfit$Ozonelag3)) & 
                          (!is.na(datfit$Ozonelag4)) & (!is.na(datfit$Ozonelag1)) &
                          (!is.na(datfit$oindex)) )

datfit$Rainfalllag1_hveg = as.numeric(datfit$ndvi700m>0.16) * datfit$Rainfalllag1 
datfit$Templag1_hveg = as.numeric(datfit$ndvi700m>0.16) * datfit$Templag1 


datfit = filter(datfit,obs>0)

timevec = list()
m = 1
for (i in unique(datfit$Year)) {
  for (j in unique(datfit$month)) {
    timevec[[m]] = data.frame(Year=i,month=j,t=m) 
    m = m + 1
  }
}
allareatime = Reduce('rbind',timevec)
datfit = datfit %>% left_join(allareatime)
areanos = apply(allareatime,1,function(x) dplyr::filter(datfit,Year==x['Year'],month==x['month'],obs>0)$areano)
datfit = filter(datfit,!(t %in% allareatime[which(sapply(areanos,length)==0),]$t))
areanos = apply(filter(allareatime,t %in% datfit$t),1,function(x) dplyr::filter(datfit,Year==x['Year'],month==x['month'],obs>0)$areano)
allaratime = filter(allareatime,t %in% datfit$t)
iiobs = sapply(areanos,function(x) as.numeric(1:64 %in% x))
datfit$t = datfit$t - length(setdiff(1:max(datfit$t),datfit$t))
# delete the first one

b0 = betareg(oindex2 ~ ndvi700m + monthf + distw + Templag1 + nino4lag1,data=datfit)
b = betareg(oindex2 ~ ndvi700m + monthf + distw + Templag1 + nino4lag1,data=datfit)
b = betareg(oindex2 ~ ndvi700m + distw + Templag1 + Rainfalllag1 + Templag4 + Rainfalllag4 + Ozonelag1 + nino4lag1 + monthf,data=datfit)
AIC(b)
AIC(b0)
summary(b)
summary(b0)


datfit$
  library(rstan)
f1 = formula(oindex2 ~ ndvi700m + t + distw + Templag1 + Templag2 + Templag3 + Templag4 + Rainfalllag1 + Rainfalllag2 + Rainfalllag3 + Rainfalllag4 + nino4lag1 + monthf)
f1 = formula(oindex2 ~ ndvi700m + t + distw + Templag1 + Rainfalllag1 + Ozonelag1 + Ozonelag2 + Rainfalllag1_hveg + Templag1_hveg + device + monthf)
modelmatrix = model.matrix(f1,datfit)

standata = list(
  N = dim(filter(datfit,obs>0))[1],
  K = dim(modelmatrix)[2],
  X = modelmatrix,
  y = filter(datfit,obs>0)$oindex2
)

filter(datfit,obs>0,is.na(oindex2))

datfit2 = datfit
datfit2$oindex2[datfit2$obs]=1
datfit2$monthf = as.factor(datfit2$month)
datfit2[8365,]
for (var in c('oindex2',setdiff(str_split(as.character(f1)[3],'\\s')[[1]],c('+','monthf')))) {
  datfit2[var][is.na(datfit2[var])] = 0
}
modelmatrix2 = model.matrix(f1,datfit2)
#dim(modelmatrix)
#select(datfit2,c('oindex2',setdiff(str_split(as.character(f1)[3],'\\s')[[1]],'+')))

distmat = fread('distmat.csv')
phi = -log(0.05)/15
cormat = exp(-phi * as.matrix(distmat))

standata4 = list(
  N = dim(filter(datfit,obs>0))[1],
  K = dim(modelmatrix)[2],
  Narea = length(unique(filter(datfit,obs>0)$areano)),
  T = length(unique(filter(datfit,obs>0)$t)),
  X = modelmatrix,
  y = filter(datfit,obs>0)$oindex2,
  cormat = cormat,
  month = filter(datfit,obs>0)$month,
  areano = filter(datfit,obs>0)$areano
)

f1 = formula(oindex2 ~ ndvi700m + t + distw + Templag1 + Rainfalllag1 + Ozonelag1 + Ozonelag2 + Rainfalllag1_hveg + Templag1_hveg + device + nino4lag1 + nino12lag1+ monthf)
modelmatrix = model.matrix(f1,datfit)

distmat = fread('distmat.csv')
phi = -log(0.05)/15
cormat = exp(-phi * as.matrix(distmat))

standata6 = list(
  N = dim(filter(datfit,obs>0))[1],
  K = dim(modelmatrix)[2],
  Narea = length(unique(filter(datfit,obs>0)$areano)),
  T = length(unique(filter(datfit,obs>0)$t)),
  X = modelmatrix,
  y = filter(datfit,obs>0)$oindex2,
  cormat = cormat,
  month = filter(datfit,obs>0)$month,
  areano = filter(datfit,obs>0)$areano
)

f1 = formula(oindex2 ~ ndvi700m + t + distw + Templag1 + Rainfalllag1 + Ozonelag1 + Ozonelag2 + Rainfalllag1_hveg + Templag1_hveg + device + nino4lag1 + nino12lag1)
modelmatrix = model.matrix(f1,datfit)

standata6nos = list(
  N = dim(filter(datfit,obs>0))[1],
  K = dim(modelmatrix)[2],
  Narea = length(unique(filter(datfit,obs>0)$areano)),
  T = length(unique(filter(datfit,obs>0)$t)),
  X = modelmatrix,
  y = filter(datfit,obs>0)$oindex2,
  cormat = cormat,
  month = filter(datfit,obs>0)$month,
  areano = filter(datfit,obs>0)$areano
)

library(rstan)

fit1 <- stan(
  file = 'betareg.stan',
  data = standata6nos,
  chains = 4,
  warmup = 1000,
  iter = 3000,
  cores = 4,
  init = '0')
loo1 = loo(fit1)
loo1

fitg <- stan(
  file = 'gaussreg.stan',
  data = standata6nos,
  chains = 4,
  warmup = 1000,
  iter = 3000,
  cores = 4,
  init = '0')

loog = loo(fitg)
loog


fit1n <- stan(
  file = 'betareg.stan',
  data = standata6,
  chains = 4,
  warmup = 1000,
  iter = 3000,
  cores = 4,
  init = '0')

#fit1val = extract(fit1)
#loo1 = loo(fit1)

#fit4 <- stan(
#  file = 'betareg4.stan',
#  data = standata4,
#  chains = 4,
#  warmup = 1000,
#  iter = 2000,
#  cores = 1,
#  init = '0')

#fit5 <- stan(
#  file = 'betareg5.stan',
#  data = standata4,
#  chains = 4,
#  warmup = 1000,
#  iter = 2000,
#  cores = 1,
#  init = '0')

fit6 <- stan(
  file = 'betareg4.stan',
  data = standata6,
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 4,
  init = '0')

loo6 = loo(fit6)
loo6

#fitp <- stan(
#  file = 'poisreg.stan',
#  data = standata6,
#  chains = 4,
#  warmup = 1000,
##  iter = 2000,
#  cores = 1,
#  init = '0')


loo1 = loo(fit1)
loo4 = loo(fit4)
loo5 = loo(fit5)
loo6 = loo(fit6)
loop = loo(fitp)

loo::compare(loo1,loog,loo6,loo1n)



###
##oindex2 ~ ndvi700m + t + distw + Templag1 + Templag2 + Templag3 + 
Templag4 + Rainfalllag1 + Rainfalllag2 + Rainfalllag3 + Rainfalllag4 + 
  nino4lag1 + monthf


Computed from 4000 by 6443 log-likelihood matrix

Estimate    SE
elpd_loo  14672.1 103.8
p_loo        27.4   0.9
looic    -29344.2 207.5
------
  Monte Carlo SE of elpd_loo is 0.1.

All Pareto k estimates are good (k < 0.5).
See help('pareto-k-diagnostic') for details.
> loo4

Computed from 4000 by 6443 log-likelihood matrix

Estimate    SE
elpd_loo  14850.3 104.8
p_loo       103.6   4.7
looic    -29700.7 209.6
------
  Monte Carlo SE of elpd_loo is 0.2.

Pareto k diagnostic values:
  Count Pct.    Min. n_eff
(-Inf, 0.5]   (good)     6441  100.0%  816       
(0.5, 0.7]   (ok)          2    0.0%  145       
(0.7, 1]   (bad)         0    0.0%  <NA>      
  (1, Inf)   (very bad)    0    0.0%  <NA>      
  
  All Pareto k estimates are ok (k < 0.7).
See help('pareto-k-diagnostic') for details.

> loo5

Computed from 4000 by 6443 log-likelihood matrix

Estimate    SE
elpd_loo  18460.3  74.3
p_loo      4024.9  34.8
looic    -36920.6 148.6
------
  Monte Carlo SE of elpd_loo is NA.

Pareto k diagnostic values:
  Count Pct.    Min. n_eff
(-Inf, 0.5]   (good)      691  10.7%   150       
(0.5, 0.7]   (ok)       1766  27.4%   52        
(0.7, 1]   (bad)      3294  51.1%   7         
(1, Inf)   (very bad)  692  10.7%   2         
See help('pareto-k-diagnostic') for details.

b0 = betareg(f1,data=datfit)


fit2val = extract(fit2)
loo2 = loo(fit2)

brmn = brms::brm(
  formula = bf(
    f1
  ),
  family = beta,
  data = datfit,
  inits = '0'
)


brm0 = brms::brm(
  formula = bf(
    oindex ~ 1 + monthf,
    zi ~ 1 + monthf
  ),
  family = zero_inflated_beta,
  data = datfit,
  inits = '0'
)

make_stancode(
  oindex2 ~ 1 + monthf,
  family = Beta,
  data = datfit)




brm1 = brms::brm(
  formula = bf(
    oindex ~ 1 + monthf + distw + ndvi + Templag1 + Rainfalllag1,
    zi ~ 1 + monthf
  ),
  family = zero_inflated_beta,
  data = datfit,
  inits = '0'
)

brm2 = brms::brm(
  formula = bf(
    oindex ~ 1 + monthf + distw + ndvi + Templag1 + Rainfalllag1 + Templag4 + Rainfalllag4,
    zi ~ 1 + monthf
  ),
  family = zero_inflated_beta,
  data = datfit,
  inits = '0'
)


summary(brm0)
loo0 = loo(brm0)
loo0
summary(brm1)
loo1 = loo(brm1)
