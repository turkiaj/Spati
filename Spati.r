
library(nlme)
library(lmfor)
data(spati)

setwd("~/projects/Biostats/Spati")

# inlude only observations with measured growth
spati<-spati[spati$id2>0,]
spati$plot <- with(spati,factor(plot))

####### LME

spati_lme <- lme(id1~id2,random=~id2|plot,data=spati,
                 correlation = corAR1())
summary(spati_lme)

#Linear mixed-effects model fit by REML
#Data: spati 
#AIC     BIC    logLik
#26267.65 26312.9 -13126.82
#
#Random effects:
#  Formula: ~id2 | plot
#Structure: General positive-definite, Log-Cholesky parametrization
#StdDev    Corr  
#(Intercept) 5.9402229 (Intr)
#id2         0.1753516 -0.623
#Residual    3.7233015       
#
#Correlation Structure: AR(1)
#Formula: ~1 | plot 
#Parameter estimate(s):
#  Phi 
#0.1199893 
#Fixed effects: id1 ~ id2 
#Value Std.Error   DF   t-value p-value
#(Intercept) 2.3678008 0.7958357 4687  2.975238  0.0029
#id2         0.7842088 0.0260120 4687 30.147903  0.0000
#Correlation: 
#  (Intr)
#id2 -0.627

####### BRMS - Tarjoaa LME4-tyyppisen rajapinnan Stanin p??lle

library(brms)
rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ())

spati_brm <- brm(formula = id1  ~ id2 + (1 + id2|plot),
                  data = spati, family = gaussian(),
                  warmup = 1000, iter = 2000, chains = 4,
                  control = list(adapt_delta = 0.95),
                  autocor = cor_ar())

stancode(spati_brm)
spati_brm  

#Family: gaussian 
#Links: mu = identity; sigma = identity 
#Formula: id1 ~ id2 + (1 + id2 | plot) 
#Data: spati (Number of observations: 4747) 
#Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
#total post-warmup samples = 4000

#Correlation Structures:
#  Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#ar[1]     0.11      0.02     0.08     0.14       4000 1.00

#Group-Level Effects: 
#  ~plot (Number of levels: 59) 
#Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#sd(Intercept)          6.01      0.61     4.96     7.30       1244 1.00
#sd(id2)                0.18      0.02     0.14     0.22       1523 1.00
#cor(Intercept,id2)    -0.60      0.09    -0.76    -0.40       2212 1.00
#
#Population-Level Effects: 
#  Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#Intercept     2.21      0.82     0.58     3.81        558 1.00
#id2           0.79      0.03     0.74     0.84       1227 1.00
#
#Family Specific Parameters: 
#  Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
#sigma     3.70      0.04     3.63     3.77       4000 1.00

# D matriisi vastaa my?s aika hyvin. T?ss? summaryss? on vaan keskihajonta eik? varianssi.

# T?st? saa generoidun Stan koodin ulos
stancode(spati_brm)

# Kaikki estimoinnit ovat jakaumia, my?s nuo ryhm?kohtaiset varianssit
plot(spati_brm)

plot(marginal_effects(spati_brm), points = TRUE)

D_brm <- VarCorr(spati_brm)
D_brm

###### RSTAN - T?ll? voi ajaa Stan-malleja R:st?

library(rstan)

# Jos koneessa on useampi prossoriydin niin t?m? nopeuttaa
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores()) 

params <- within(list(),
                 {
                   N <- nrow(spati)
                   Y <- as.vector(spati$id1)
                   X <- cbind(1,as.matrix(spati$id2))   
                   Z <- X                               # Z is indexed by 'group' so it's not block diagonal matrix
                   p <- 2                               # fixed-effects
                   k <- 2                               # random-effects
                   J <- length(levels(spati$plot))
                   group <- as.integer(spati$plot)      # group index for data rows
                 })

#spati_model <- stan_model(file = "spati.stan")

spati_fit <- stan(file="spati.stan", data=params, warmup=1000, iter=2000, chains=4, control = list(adapt_delta = 0.95))

saveRDS(spati_fit, file="models/spati_BLMM.rds")

spati_fit <- readRDS("models/spati_BLMM.rds")

print(spati_fit, pars = c("beta[1]", "beta_Intercept", "sigma_e", "sigma_b[1]", "sigma_b[2]", "ar1"))
summary(spati_lme)
spati_fit

#                mean se_mean   sd 2.5%  25%  50%  75% 97.5% n_eff Rhat
# beta[1]        0.81    0.00 0.03 0.75 0.79 0.81 0.82  0.86  1023    1
# beta_Intercept 1.93    0.04 0.79 0.42 1.39 1.92 2.44  3.53   495    1
# sigma_e        3.71    0.00 0.04 3.64 3.68 3.71 3.74  3.79  4000    1
# sigma_b[1]     5.73    0.01 0.56 4.73 5.33 5.70 6.10  6.94  1542    1
# sigma_b[2]     0.18    0.00 0.02 0.14 0.16 0.18 0.19  0.22  1296    1

# Vertaillaan random-efektejä -- nämä ovat samoja
ranef(spati_lme)
print(spati_fit, pars = c("b"))

b_posterior <- extract(spati_fit, pars = c("b"))$b
ranef_blmm <- colMeans(b_posterior)

# Stanilla voi my?s approksimoida jakauman variational bayes-menetelm?ll?

spati_model <- stan_model(file = "spati.stan")
spati_fit_vb <- vb(spati_model, data=params, output_samples=3000, iter=4000, seed=123)

# T?m? on huomattavasti nopeampi kuin HMC-samplays, mutta estimaatit ovatkin sitten sinne p?in..
print(spati_fit_vb, pars = c("beta[1]", "beta_Intercept", "sigma_e", "sigma_b[1]", "sigma_b[2]"))

#                mean   sd 2.5%  25%  50%  75% 97.5%
# beta[1]        0.70 0.01 0.68 0.69 0.70 0.70  0.71
# beta_Intercept 1.25 0.12 1.02 1.17 1.25 1.33  1.49
# sigma_e        4.22 0.04 4.14 4.19 4.22 4.25  4.30
# sigma_b[1]     0.01 0.01 0.00 0.01 0.01 0.02  0.05
# sigma_b[2]     0.17 0.00 0.17 0.17 0.17 0.17  0.18

#### dens overlay

library(bayesplot)

posterior <- extract(spati_fit, pars = c("Y_rep"))
posterior_y_50 <- posterior$Y_rep[1:50,]

plot(ppc_dens_overlay(params$Y, posterior_y_50) + 
       coord_cartesian(xlim = c(-1, 50)))

### TESTATAAN ENNUSTEITA SPATI-DATALLA

source("../mebn/MEBN.r")

localsummary <- mebn.localsummary(spati_fit)
err <- mebn.ranef_BLUP_vars(params, localsummary)

localsummary

# - pistokoe

err$C1 # t?m?n pit?isi olla var(beta.hat) eli sama kuin ...?

#[,1]         [,2]
#[1,]  0.376100103 -0.009222113
#[2,] -0.009222113  0.000552146

# eli betan varianssi on [0.37610, 0.00055] int beta

beta_var <- diag(err$C1)
beta_sd <- sqrt(beta_var)

#[1] 0.61327001 0.02349779

# Pit?isi olla
# 0.79  0.03  (stan)
# 0.80  0.03  (brms, pop-level eff. est. error)
# 0.75  0.025 (nlme)

# DEBUG---
#mebn.ranef_BLUP_vars <- function(localparams, localsummary)
#{

# Number of fixed and random effects 
p <- params$p
k <- params$k

n <- params$N  # used to calculate R
groups <- length(unique(params$group))
group_obs <- n/groups

# Here X and Z are design matrices for whole data, not just one group
X <- params$X     # fixef design matrix 
Z <- mebn.ranef_designmatrix(X,groups,group_obs)
Y <- as.vector(spati$id1)

dim(Z)

group_obs

# We also need D for whole data (D on ok!!)
Di <- as.matrix(localsummary$D)
D <- diag(1,groups) %x% Di # block diagonal matrix for every group with Kronecker product

sigma<-localsummary$std_error

R <- diag(sigma^2,nrow=n,ncol=n) # residual var-cov matrix 

# Henderson's Mixed-Effect Equations

H11 <- t(X)%*%solve(R)%*%X
H12 <- t(X)%*%solve(R)%*%Z
H22 <- t(Z)%*%solve(R)%*%Z+solve(D)

H <- cbind(rbind(H11,t(H12)),rbind(H12,H22))

r<-rbind(t(X)%*%solve(R)%*%Y,t(Z)%*%solve(R)%*%Y)
solve(H)%*%r

r

# pit?isi tulla 

#beta_int
#beta
#ranef_int
#ranef

# The inverse of H provides variance-covariance matrix of the estimation errors
Hinv <- solve(H)

# Submatrices of Hinv
C1 <- Hinv[1:p,1:p]                   # var(beta.hat) 
C2 <- Hinv[(p+1):(p+k),(p+1):(p+k)]   # var(b.tilde-b), var(BLUP)
C12 <- Hinv[(p+1):(p+k),1:p]          # cov(beta.hat, b.tilde-b)

beta_var <- diag(C1)
beta_sd <- sqrt(beta_var)
beta_sd

# [1] 0.61327001 0.02349779 --- ei sinnep?ink??n? pit?isi olla

# 0.79  0.03  (stan)

