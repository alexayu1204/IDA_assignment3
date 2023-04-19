install.packages("mice")
install.packages("ggplot2")
library(mice)
library(ggplot2)
library(JointAI)
require(reshape2)
require(RColorBrewer)
require(ggplot2)
require(mitml)

#Q1a
cc(nhanes)
nrow(cc(nhanes))
incom_percent <- (nrow(nhanes)-nrow(cc(nhanes)))/nrow(nhanes)
incom_percent

#Q1b
imps<- mice(nhanes, printFlag = FALSE, seed = 1)
imps
#predict bmi from age, hyp, and chl by the normal linear regression model
fits<- with(imps, lm(bmi ~ age + hyp + chl))
class(fits)
#result
pool_est <- pool(fits)
pool_est

#Q1c
#Repeat the analysis for seed2,3,4,5,6
est_2<-pool(with(mice(nhanes, printFlag = FALSE, seed = 2),
                  lm(bmi ~ age + hyp + chl)))
est_2
est_3<-pool(with(mice(nhanes, printFlag = FALSE, seed = 3),
                  lm(bmi ~ age + hyp + chl)))
est_3
est_4<-pool(with(mice(nhanes, printFlag = FALSE, seed = 4),
                  lm(bmi ~ age + hyp + chl)))
est_4
est_5<-pool(with(mice(nhanes, printFlag = FALSE, seed = 5),
                  lm(bmi ~ age + hyp + chl)))
est_5
est_6<-pool(with(mice(nhanes, printFlag = FALSE, seed = 6),
                  lm(bmi ~ age + hyp + chl)))
est_6

#Q1d
#Repeat the analysis with M = 100 with the same seeds
ests_1<-pool(with(mice(nhanes, printFlag = FALSE, seed = 1, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_1
ests_2<-pool(with(mice(nhanes, printFlag = FALSE, seed = 2, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_2
ests_3<-pool(with(mice(nhanes, printFlag = FALSE, seed = 3, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_3
ests_4<-pool(with(mice(nhanes, printFlag = FALSE, seed = 4, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_4
ests_5<-pool(with(mice(nhanes, printFlag = FALSE, seed = 5, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_5
ests_6<-pool(with(mice(nhanes, printFlag = FALSE, seed = 6, m = 100),
                  lm(bmi ~ age + hyp + chl)))
ests_6



#Q2
#Load data
load("dataex2.Rdata")

#Initialize counters
counter_sri <- 0
counter_bootstrap <- 0

#Number of datasets
n_datasets <- dim(dataex2)[3]

#Loop through datasets
for (i in 1:n_datasets) {
  # Stochastic Regression Imputation
  sri_imputations <- mice(dataex2[,,i], method = "norm.nob",
                          m = 20, seed = 1, printFlag = FALSE)
  sri_fits <- with(sri_imputations, lm(Y~X))
  sri_estimates <- pool(sri_fits)
  sri_summary <- summary(sri_estimates, conf.int = TRUE)
  if (sri_summary[2, c(7)] <= 3 & sri_summary[2, c(8)] >= 3) {
    counter_sri <- counter_sri + 1
  }
  
  # Bootstrap-based Version
  bootstrap_imputations <- mice(dataex2[,,i], method = "norm.boot",
                                m = 20, seed = 1, printFlag = FALSE)
  bootstrap_fits <- with(bootstrap_imputations, lm(Y~X))
  bootstrap_estimates <- pool(bootstrap_fits)
  bootstrap_summary <- summary(bootstrap_estimates, conf.int = TRUE)
  if (bootstrap_summary[2, c(7)] <= 3 & bootstrap_summary[2, c(8)] >= 3) {
    counter_bootstrap <- counter_bootstrap + 1
  }
}

#Calculate empirical coverage probabilities
ecp_sri <- counter_sri/n_datasets
ecp_bootstrap <- counter_bootstrap/n_datasets

#Print results
cat("Empirical coverage probability for Stochastic Regression Imputation: ", 
    ecp_sri, "\n")
cat("Empirical coverage probability for Bootstrap-based Version: ", 
    ecp_bootstrap, "\n")


#Q4a
load('dataex4.Rdata') # Load data
q4_imp0 = mice(dataex4, m=50, seed=1, printFlag=FALSE, maxit=0)
#preventing x2 from being imputed
q4_imp0$predictorMatrix["x2",] = 0 
#Use MICE with new predictor Matrix
q4_imps = mice(dataex4, m=50, seed=1, printFlag=FALSE,
                predictorMatrix=q4_imp0$predictorMatrix) 
q4_ests = with(q4_imps, lm(y ~ x1 + x2 + x1*x2)) 
# Extracting relevant statistics
summary(pool(q4_ests), conf.int=TRUE)[-1,c(1,2,7,8)] 

#Q4b
#Create the interaction variable (z=x1*x2)
dataex4$z = dataex4$x1 * dataex4$x2
#MICE setup
imp_0_pass = mice(dataex4, maxit = 0, seed=1, m=50)
#Specify that z is derived from x1 and x3
passmeth = imp_0_pass$method
passmeth["z"] = "~I(x1*x2)"
passpred = imp_0_pass$predictorMatrix
#won't use z to impute x1 and x2
passpred[c("x1", "x2"), "z"] = 0
#won't use y to impute z
passpred["z", "y"] = 0 
#then we start imputation
imp_pass = mice(dataex4, 
                method = passmeth, 
                predictorMatrix = passpred, 
                m = 50, 
                seed = 1, 
                printFlag = FALSE)
imp_pass_ests = with(imp_pass, lm(y ~ x1 + x2 + z))
summary(pool(imp_pass_ests), conf.int=TRUE)[-1,c(1,2,7,8)] 

#4c
imps_4c = mice(dataex4, 
               m=50, 
               seed=1, 
               printFlag=FALSE)
ests_4c = with(imps_4c, 
               lm(y ~ x1 + x2 + z))
summary(pool(ests_4c), 
        conf.int=TRUE)[-1,c(1,2,7,8)]

#Q5
#Load necessary libraries and data
load("NHANES2.Rdata")
require(JointAI)
library(JointAI)
require(devtools)
require(reshape2)
require(RColorBrewer)
require(ggplot2)
source_url("https://gist.githubusercontent.com/NErler/0d00375da460dd33839b98faeee2fdab/raw/c6f537ecf80eddcefd94992ec7926aa57d454536/propplot.R")

#get some insights from the data set and missing data
dim(NHANES2)
str(NHANES2)
summary(NHANES2)
mdpat_mice <- md.pattern(NHANES2)
mdpat_mice
md_pattern(NHANES2, pattern = FALSE, color = c('#34111b', '#e30f41'))

par(mar = c(3,3,2,1), mgp = c(2, 0.6, 0))
plot_all(NHANES2, breaks = 30, ncol = 4)

#impute using mice(defalut)
imp0 <- mice(NHANES2, maxit = 0)
imp0

meth <- imp0$method
meth["hgt"] <- "norm"
meth

post <- imp0$post
post["hgt"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i], c(0, 2.5))"

#impute using mice with maxit = 20, m = 20
imp <- mice(NHANES2, method = meth,
            maxit = 20, m = 20, seed = 1, printFlag = FALSE)
imp$loggedEvents

plot(imp, layout = c(4,4))
densityplot(imp)
xyplot(imp, hgt ~ wgt | gender, pch = c(1, 20))


#imputation analysis
fit<- with(imp, lm(wgt~gender+age+hgt+WC))
summary(fit$analyses[[1]])
comp1 <- complete(imp, 1)

#create an residual plot
plot(fit$analyses[[1]]$fitted.values, residuals(fit$analyses[[1]]),
     xlab = "Fitted values", ylab = "Residuals")

#plot wgt against every variable
plot(comp1$wgt ~ comp1$gender, xlab = "gender", ylab = "wgt")
plot(comp1$wgt ~ comp1$age, xlab = "age", ylab = "wgt")
plot(comp1$wgt ~ comp1$hgt, xlab = "Height in metres", ylab = "wgt")
plot(comp1$wgt ~ comp1$WC, xlab = "Waist circumference in cm", ylab = "wgt")

#check the normal QQ plot
qqnorm(rstandard(fit$analyses[[1]]), xlim = c(-4, 4), ylim = c(-6, 6))
qqline(rstandard(fit$analyses[[1]]), col = 2)

pooled_ests <- pool(fit)
summary(pooled_ests, conf.int = TRUE)

pool.r.squared(pooled_ests, adjusted = TRUE)

fit_no_gender <- with(imp, lm(wgt ~ age + hgt + WC))
D1(fit, fit_no_gender)

fit_no_age<- with(imp, lm(wgt ~ gender + hgt + WC))
D1(fit, fit_no_age)

fit_no_hgt<- with(imp, lm(wgt ~ age+ gender + WC))
D1(fit, fit_no_hgt)

fit_no_WC<- with(imp, lm(wgt ~ age + gender + hgt))
D1(fit, fit_no_WC)