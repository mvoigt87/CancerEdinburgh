###
### A file for creating a mask for the final analysis
### - Part 1 - Time Series Analysis for stage at diagnosis -
###

### 0.1 - Load Packages
library(tidyverse)
library(forcats)
library(data.table)
library(broom)
library(survival)
require(MASS)
require(Hmisc)

### 0.2 - Load Data (load any data set)
# set path to data
setwd("C:/Users/.../RCodeData")
load('BREAST.Rdata')

# give the data set a short name
d <- breastcure

# 1. Descriptives

table(d$sex1)
# low case number will probably lead to the exclusion of male cancer patients
d <- d %>% filter(sex1=="female")

### 1.1. variable of main interest - stage at diagnosis
table(d$dxstate1)
round(prop.table(table(d$dxstate1)),3)

### 1.1.2 Age at diagnosis/stage at diagnosis
  # probably we will need to create meaningful age groups - using a age group variable (ageGr)

table(d$ageGr)
round(prop.table(table(d$ageGr)),3)
round(prop.table(table(d$dxstate1, d$ageGr),2),3)


### 2. General Variables for Analysis
#####################################

# 0-1 Variable for stage at diagnosis
d <- d %>% mutate(diag2 = ifelse(DxState<=1,"early","late"))
class(d$DxState) # in case it is not a factor = change
d$diag2 <- as.factor(d$diag2)

# Year by which to split the time series
d1 <- d %>% filter(DxYear<=2000)
d2 <- d %>% filter(DxYear>2000)


### 3. Missing Value treatment
##############################
d %>% count(is.na(d$Age))          # data.table package might be better 

d$Age[is.na(d$Age)] <- mean(d$Age,na.rm=T) # set to mean age


### Change reference categories for covariate analysis
### --------------------------------------------------

# Age at diagnosis
# d <- d %>% mutate(age = (DxYear+(DxMon/12)) - BirthYear)

# Environmental variables

# SES (Education)

# marital status variable
class(d$marstat1)
d <- within(d, marstat1 <- relevel(marstat1, ref = "married")) 



### 4. Logistic regression - (early vs. late diagnosis)
#######################################################

model_1 <- glm(diag2 ~ marstat1 + ageGr ,           # add variables as it goes
               family=binomial(link='logit'),
               data=d) 

summary(model_1)

# Obtain the odds ratio
coeff <- exp(model_1$coefficients)

# 4.2 Model test
anova(model_1, test="Chisq")


# 4.3 Comparing different times

model_1A <- glm(diag2 ~ marstat1 + ageGr ,           
               family=binomial(link='logit'),
               data=d1)                         # before 2000

summary(model_1A)

model_1B <- glm(diag2 ~ marstat1 + ageGr ,           
                family=binomial(link='logit'),
                data=d2)                        # after 2000

summary(model_1B)


### 5. Ordered logistic regression 
###################################
  # (outcome can be a an factor variable with ordered categories from little to much)

# prepare the outcome
class(d$dxstate1)

## fit ordered logit model and store results 'm'
m <- polr(dxstate1 ~ marstat1 + ageGr, data = d, Hess=TRUE)

## view a summary of the model
summary(m)

# Help for interpretation from: https://stats.idre.ucla.edu/r/dae/ordinal-logistic-regression/ 
# Next we see the usual regression output coefficient table including the value of each coefficient, standard errors, and t value, which is simply the ratio of the coefficient to its standard error. There is no significance test by default.
# Next we see the estimates for the two intercepts, which are sometimes called cutpoints. The intercepts indicate where the latent variable is cut to make the three groups that we observe in our data. Note that this latent variable is continuous. In general, these are not used in the interpretation of the results. The cutpoints are closely related to thresholds, which are reported by other statistical packages.
# Finally, we see the residual deviance, -2 * Log Likelihood of the model as well as the AIC. Both the deviance and AIC are useful for model comparison.

## Odds ratio, coefficients, and p-values
exp(coef(m))

(ctable <- coef(summary(m)))

p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2

## combined table
(ctable <- cbind(ctable, "p value" = p))

## normally distributed CIs
confint.default(m)

## OR and CI
(ci <- confint(m))
exp(cbind(OR = coef(m), ci))


### 6. Time Trend Analysis
##########################
X <- aggregate(d$diag2=="late", by=list(Year=d$DxYear, AgeGR=d$ageGr), FUN=sum)
X.2 <- aggregate(d$diag2=="early", by=list(Year=d$DxYear, AgeGR=d$ageGr), FUN=sum)

X <- X %>% inner_join(X.2, by= c("Year", "AgeGR")) %>% mutate(perc = x.x/(x.x+x.y))

X %>% ggplot(aes(x=Year, y=perc, color=AgeGR)) +
  geom_line() +                                # looks a little noisy with many age groups and few cases
  scale_y_continuous(name = "% late diagnosis")

# plot(X$Year,X$x, type = "l")
