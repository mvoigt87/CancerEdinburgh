###
### A file for creating a mask for the final analysis
### - Test with SEER data -
###

# 0.1 - Load Packages
library(tidyverse)
library(forcats)
library(data.table)
library(broom)
library(survival)

# 0.2 - Load Data (load any data set)
setwd("C:/Users/Mathias/Documents/LongPop_Madrid/Secondment/Edinburgh/RCodeData")
load('BREAST.Rdata')

d <- breastcure # fill in data set name for breastcure

# 1. Checking for the variables of interest

# 1.1. Time/Age information

# time between diagonosis and censorship/event

d <- d %>% mutate(time.d=Survmonths)  
    # (depending on the way the information is coded it might be necessary to transform the time information)

# age at censorship/death (for an analysis with left-truncated data)
d <- d %>% mutate(age.d = (DxYear+(DxMon/12)+(Survmonths/12))-BirthYear) %>% 
 # this approach requires to have an entry age as well (call it age.e - find age variable in the data set)
  mutate(age.e = Age)

# quick check
summary(d$age.d)
hist(d$age.d)


# d %>% mutate(event = as.factor(ifelse(event==1,"dead","alive"))) %>% ggplot(aes(d$age.d, fill=event)) +  
#   geom_histogram(binwidth = 0.1) + 
#   scale_fill_discrete(name=" ") + 
#   scale_x_continuous(name=" ") +
#   scale_y_continuous(name=" ") +
#   theme_bw()

# Creating different survival objects
# -----------------------------------

# time since diagnosis
# ....................

# Surv(time = d$time.d, event = d$event)

# time since diagnosis
# ....................

# Surv(time = d$age.e, time2 = d$age.d, event = d$event)




# 2. Exploratory Kaplan Meier Plots
# ---------------------------------

# 2.1 No covariates + time since diagnosis

KM1 <- survfit(coxph(Surv(time = d$time.d, event = d$event)~1, data = d),type="kaplan-meier")
# KM1
# summary(KM1)

# plot
tidy(KM1) %>% ggplot(aes(x=time,y=estimate)) +
  geom_step() +
  scale_x_continuous(name = "Time Since Diagnosis in Months") +
  scale_y_continuous(name = "Estimated Survival Probability") +
  theme_bw()


# 2.2 No covariates - Age at death as time dimension (left truncated data)

KM2 <- survfit(coxph(Surv(time = d$age.e, time2 = d$age.d, event = d$event)~1, data = d),type="kaplan-meier")
# KM2
# summary(KM2)

# plot - age patterns with high mortality in younger ages
tidy(KM2) %>% ggplot(aes(x=time,y=estimate)) +
  geom_step() +
  scale_x_continuous(name = "Time Since Diagnosis in Months") +
  scale_y_continuous(name = "Estimated Survival Probability") +
  theme_bw()


# 2.3 Example Kaplan-Meier (time since diagnosis) with covariate effects (example - stage at diagnosis)

# extract different levels
table(d$dxstate1)

# KMEs
KM3_1 <- survfit(Surv(time = time.d,
                      event = event) ~ 1, data = subset(d,dxstate1=="in situ"),
                       type="kaplan-meier")
KM3_2 <- survfit(Surv(time = time.d,
                      event = event) ~ 1, data = subset(d,dxstate1=="localized"),
                 type="kaplan-meier")
KM3_3 <- survfit(Surv(time = time.d,
                      event = event) ~ 1, data = subset(d,dxstate1=="regional"),
                 type="kaplan-meier")
KM3_4 <- survfit(Surv(time = time.d,
                      event = event) ~ 1, data = subset(d,dxstate1=="distant"),
                 type="kaplan-meier")

# Cleaning for ggplot (change variable of interest)
     
KM3_A <- tidy(KM3_1) %>% dplyr::select(time,estimate) %>% mutate(dxstate = "in situ")
KM3_B <- tidy(KM3_2) %>% dplyr::select(time,estimate) %>% mutate(dxstate = "localized")
KM3_C <- tidy(KM3_3) %>% dplyr::select(time,estimate) %>% mutate(dxstate = "regional")
KM3_D <- tidy(KM3_4) %>% dplyr::select(time,estimate) %>% mutate(dxstate = "distant")

KM3 <- union(KM3_A, KM3_B) %>% union(KM3_C) %>% union(KM3_D)

# plot
plot_KM3 <- KM3 %>% ggplot(aes(x=time,y=estimate, color=dxstate)) +
  geom_step() +
  scale_x_continuous(name = "Time Since Diagnosis in Months") +
  scale_y_continuous(name = "Estimated Survival Probability") +
  scale_color_manual(values=c("green2", "goldenrod2", "dodgerblue1", "orangered"), name=" ",
                     breaks=c("in situ","localized","regional","distant")) +                  # change colors and values
  theme_bw()

# move the legend
plot_KM3 <- plot_KM3 + theme(legend.position = c(0.85, 0.25)) +      # move legend position around by x,y coordinate changes
  scale_shape_discrete(guide=FALSE)
plot_KM3



### Change reference categories for covariate analysis
### --------------------------------------------------

# sex variable
class(d$sex1)
d <- within(d, sex1 <- relevel(sex1, ref = "female"))  
# surgery variable
class(d$surg1)
d <- within(d, surg1 <- relevel(surg1, ref = "no surgery")) 
# diagnosis variable
class(d$dxstate1)
d <- within(d, dxstate1 <- relevel(dxstate1, ref = "in situ")) 

# 3. Survival Models
# ------------------


# 3.1. Cox Model

# 3.1.1 Model with covariate of interest (i.e. DxState)

fit_1 <- coxph(Surv(time = time.d,
                    event = event) ~ dxstate1, data=d)

summary(fit_1)
# test for proportional hazards
cox.zph(fit_1)

# 3.1.1 Model with more covariates

fit_2 <- coxph(Surv(time = time.d,
                    event = event) ~ dxstate1 + surg1 + ageGR + marstat1, data=d)

summary(fit_2)
# test for proportional hazards
cox.zph(fit_2)




