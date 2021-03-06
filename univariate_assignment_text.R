

trees <- read.csv('https://raw.githubusercontent.com/dmcglinn/quant_methods/gh-pages/data/treedata_subset.csv')
trees
# 1. Carry out an exploratory analysis using the tree dataset. Metadata for the tree study can be found here. Specifically, I would like you to develop and compare models for species cover for a habitat generalist Acer rubrum (Red maple) and a habitat specialist Abies fraseri (Frasier fir). Because this dataset includes both continuous and discrete explanatory variables use the function Anova in the packages car as such
library(car)
# Anova(my_mod, type=3)
# Compare the p-values you observe using the function Anova to those generated using summary.

# For each species address the following additional questions:

# how well does the exploratory model appear to explain cover?
# which explanatory variables are the most important?
# do model diagnostics indicate any problems with violations of OLS assumptions?... answered later in this document.
# are you able to explain variance in one species better than another, why might this be the case?
# Prior to addressing the above questions you will want to restructure and subset the data using the # # following R code:

# we wish to model species cover across all sampled plots
# create site x sp matrix for two species 
sp_cov = with(trees, tapply(cover, list(plotID, spcode), 
                            function(x) round(mean(x))))
sp_cov = ifelse(is.na(sp_cov), 0, sp_cov)
sp_cov = data.frame(plotID = row.names(sp_cov), sp_cov)
# create environmental matrix
cols_to_select = c('elev', 'tci', 'streamdist', 'disturb', 'beers')
env = aggregate(trees[ , cols_to_select], by = list(trees$plotID), 
                function(x) x[1])
names(env)[1] = 'plotID'
# merge species and environmental matrices
site_dat = merge(sp_cov, env, by='plotID')
# subset species of interest
abies = site_dat[ , c('ABIEFRA', cols_to_select)]
acer  = site_dat[ , c('ACERRUB', cols_to_select)]
names(abies)[1] = 'cover'
names(acer)[1] = 'cover'

# loading more libraries
library(ggplot2)
library(gridExtra)
library(scatterplot3d)
library(MASS)

# viewing the data
abies
acer
names(abies)
names(acer)

# boxplots for quantitative----
ggplot(data = abies) + 
  geom_boxplot(mapping = aes(x = disturb, y = cover)) +  
  labs(x = 'Disturbance', y = 'Cover', title = 'abies disturbance versus cover') 

ggplot(data = acer) + 
  geom_boxplot(mapping = aes(x = disturb, y = cover)) +  
  labs(x = 'Disturbance', y = 'Cover', title = 'acer disturbance versus cover') 


quantile(abies$cover[abies$disturb == 'VIRGIN'])
quantile(abies$cover[abies$disturb == 'CORPLOG'])
quantile(abies$cover[abies$disturb == 'SETTLE'])
quantile(abies$cover[abies$disturb == 'LT-SEL'])

quantile(acer$cover[acer$disturb == 'VIRGIN'])
quantile(acer$cover[acer$disturb == 'CORPLOG'])
quantile(acer$cover[acer$disturb == 'SETTLE'])
quantile(acer$cover[acer$disturb == 'LT-SEL'])

#ABIES original plots----
plot(cover ~ elev, data = abies, xlab = 'Elevation (m)',
     ylab = 'Cover', main= 'abies elevation versus cover')

plot(cover ~ tci, data = abies, xlab = 'Topographic Coverage Index',
     ylab = 'Cover', main= 'abies tci versus cover')

plot(cover ~ streamdist, data = abies, xlab = 'Stream Distance (m)',
     ylab = 'Cover', main= 'abies stream distance versus cover')

plot(cover ~ beers, data = abies, xlab = 'Beers',
     ylab = 'Cover', main= 'abies beers versus cover')

#ACER original plots----
plot(cover ~ elev, data = acer, xlab = 'Elevation (m)',
     ylab = 'Cover', main= 'acer elevation versus cover')

plot(cover ~ tci, data = acer, xlab = 'Topographic Coverage Index',
     ylab = 'Cover', main= 'acer tci versus cover')

plot(cover ~ streamdist, data = acer, xlab = 'Stream Distance (m)',
     ylab = 'Cover', main= 'acer stream distance versus cover')

plot(cover ~ beers, data = acer, xlab = 'Beers',
     ylab = 'Cover', main= 'acer beers versus cover')

# Intercept only models for Abies and Acer, cover as y variable----
null_mod_abies = lm(cover ~ 1, data = abies)
null_mod_acer = lm(cover ~ 1, data = acer)     
null_mod_abies
null_mod_acer

mean(abies$cover)
mean(acer$cover)

plot(cover ~ 1, data = abies)
abline(null_mod_abies, lwd = 2)
abline(h = mean(abies$cover), col = 'red', lty = 2, lwd = 2)

plot(cover ~ 1, data = acer)
abline(null_mod_acer, lwd = 2)
abline(h = mean(acer$cover), col = 'red', lty = 2, lwd = 2)

# create single variable main effect models for each----
# abies
elev_mod_abies = lm(cover ~ elev, data = abies)
tci_mod_abies = lm(cover ~ tci, data = abies)
streamdist_mod_abies = lm(cover ~ streamdist, data = abies)
disturb_mod_abies = lm(cover ~ disturb, data = abies)
beers_mod_abies = lm(cover ~ beers, data = abies)

# acer 
elev_mod_acer = lm(cover ~ elev, data = acer)
tci_mod_acer = lm(cover ~ tci, data = acer)
streamdist_mod_acer = lm(cover ~ streamdist, data = acer)
disturb_mod_acer = lm(cover ~ disturb, data = acer)
beers_mod_acer = lm(cover ~ beers, data = acer)

# abies single model main effect summaries ----
summary(elev_mod_abies)
Anova(elev_mod_abies, type= 3)
summary(tci_mod_abies)
summary(streamdist_mod_abies)
summary(disturb_mod_abies)
summary(beers_mod_abies)

# significant p-values = elev, disturb (b/c of VIRGIN), streamdist

# acer single model main effect summaries----
summary(elev_mod_acer)
summary(tci_mod_acer)
summary(streamdist_mod_acer)
summary(disturb_mod_acer)
summary(beers_mod_acer)

# significant p-values = all EXCEPT tci

# beers has a significant effect on acers, but not abies in these models

# all main effects for both----
all_mod_abies = lm(cover ~ elev + tci + streamdist + disturb + beers, data = abies) 
all_mod_acer = lm(cover ~ elev + tci + streamdist + disturb + beers, data = acer)

summary(all_mod_abies)
summary(all_mod_acer)

# identifying outlers and modifying acer because of them----
# outliers in this plot tci vs acer cover
# used this code to identify the outliers... identify(acer$cover ~ acer$tci, n = 2)
# remove tci outliers
acer_subset = acer[-c(121, 318), ]
plot(cover ~ tci, data= acer_subset, main="acer subset tci versus cover", xlab='tci', ylab= 'cover')

dim(acer)
dim(acer_subset)

# plots redone here for convenience... can take out in final probably----
plot(cover ~ elev, data = acer, xlab = 'Elevation (m)',
     ylab = 'Cover', main= 'acer elevation versus cover')

plot(cover ~ tci, data = acer, xlab = 'Topographic Coverage Index',
     ylab = 'Cover', main= 'acer tci versus cover')

plot(cover ~ streamdist, data = acer, xlab = 'Stream Distance (m)',
     ylab = 'Cover', main= 'acer stream distance versus cover')

plot(cover ~ beers, data = acer, xlab = 'Beers',
     ylab = 'Cover', main= 'acer beers versus cover')

# remove possible stream distance acer outlier
# identify(acer$cover ~ acer$streamdist, n=1)
acer_subset = acer[-c(121, 318, 297,56,187), ]
dim(acer_subset)

plot(cover ~ streamdist, data = acer_subset, xlab = 'Stream Distance (m)',
     ylab = 'Cover', main= 'acer stream distance versus cover subset')

# new single variable models with the acersubset----
elev_mod_acersubset = lm(cover ~ elev, data = acer_subset)
tci_mod_acersubset = lm(cover ~ tci, data = acer_subset)
streamdist_mod_acersubset = lm(cover ~ streamdist, data = acer_subset)
disturb_mod_acersubset = lm(cover ~ disturb, data = acer_subset)
beers_mod_acersubset = lm(cover ~ beers, data = acer_subset)

# comparing the original to the subset without outliers----
summary(elev_mod_acer)
summary(elev_mod_acersubset)
summary(tci_mod_acer)
summary(tci_mod_acersubset)
summary(streamdist_mod_acer)
summary(streamdist_mod_acersubset)
summary(disturb_mod_acer)
summary(disturb_mod_acersubset)
summary(beers_mod_acer)
summary(beers_mod_acersubset)

# comparing the original main model with the subset main model----
all_mod_acersubset = lm(cover ~ elev + tci + streamdist + disturb + beers, data = acer_subset)
summary(all_mod_acersubset)
summary(all_mod_acer)

# removing the outliers did slightly improve the R-squared and adjusted R- squared, it also greatly increased the significance of tci( lowered the p value), so the subset without the outliers will be used from now on

# model with elev, tci, and beers----
etb_mod_acersubset = lm(cover ~ elev + tci + beers , data = acer_subset)
summary(etb_mod_acersubset)

etb_interaction_mod_acersubset = lm(cover ~ elev + tci + beers + tci * elev + elev * beers + beers * tci + beers * tci * elev , data = acer_subset)
summary(etb_interaction_mod_acersubset)

# interaction effect with elev:tci and elev:beers
AIC(etb_interaction_mod_acersubset)
#[1] 3413.76
AIC(etb_mod_acersubset)
#[1] 3419.997
## the above model just proved my suspicions on how important elevation was to the model, it will not be used over other models however because the r squared was lower

#adding elevation interaction to main model-------------
elev_interaction_mod_acersubset = lm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers + elev * streamdist + elev * disturb, data = acer_subset)
summary(elev_interaction_mod_acersubset)

elev_interaction_mod_acersubset = lm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers + elev * streamdist + elev * disturb, data = acer_subset)

# adding interaction to whole model----
full_mod_acersubset = update(all_mod_acersubset, ~ . + elev * disturb * tci * streamdist * beers)
all_mod_acersubset = lm(cover ~ elev + tci + streamdist + disturb + beers, data = acer_subset)

summary(all_mod_acersubset)
summary(elev_interaction_mod_acersubset)
summary(full_mod_acersubset) 
# this makes it absurdly long, but it does increase the r-squared...


anova(all_mod_acersubset, full_mod_acersubset)

AIC(full_mod_acersubset)
AIC(all_mod_acersubset)
AIC(elev_interaction_mod_acersubset)
# ele

# remove abies outlier in tci and elevation  
# identify(abies$cover ~ abies$tci, n=2)

abies_subset = abies[-c(121,33,109,56,53), ]
dim(abies_subset)

# comparing the main models of original and subset
all_mod_abiessubset = lm(cover ~ elev + tci + streamdist + disturb + beers, data = abies_subset)
summary(all_mod_abiessubset)
summary(all_mod_abies)


# abies with interaction, full and elevation----
full_mod_abiessubset = update(all_mod_abiessubset, ~ . + elev * disturb * tci * streamdist * beers)
elev_interaction_mod_abiessubset = lm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers + elev * streamdist + elev * disturb, data = abies_subset)

summary(all_mod_abiessubset)
summary(elev_interaction_mod_abiessubset)
summary(full_mod_abiessubset) 

# interaction had an even larger effect on the abies distributions


# Poisson calculation start----
# assume all poisson calculations are done using the subsets with outliers removed
summary(elev_interaction_mod_abiessubset)
summary(elev_interaction_mod_acersubset)

# poisson individual variables acer
elev_mod_acer_poi = glm(cover ~ elev, data = acer_subset, family = "poisson")
tci_mod_acer_poi = glm(cover ~ tci, data = acer_subset, family = "poisson")
streamdist_mod_acer_poi = glm(cover ~ streamdist, data = acer_subset, family = "poisson")
disturb_mod_acer_poi = glm(cover ~ disturb, data = acer_subset, family = "poisson")
beers_mod_acer_poi = glm(cover ~ beers, data = acer_subset, family = "poisson")

# poisson individual variables abies 
elev_mod_abies_poi = glm(cover ~ elev, data = abies_subset, family = "poisson")
tci_mod_abies_poi = glm(cover ~ tci, data = abies_subset, family = "poisson")
streamdist_mod_abies_poi = glm(cover ~ streamdist, data = abies_subset, family = "poisson")
disturb_mod_abies_poi = glm(cover ~ disturb, data = abies_subset, family = "poisson")
beers_mod_abies_poi = glm(cover ~ beers, data = abies_subset, family = "poisson")

# summaries of individual poisson models
# abies
summary(elev_mod_abies_poi)
summary(tci_mod_abies_poi)
summary(streamdist_mod_abies_poi)
summary(disturb_mod_abies_poi)
summary(beers_mod_abies_poi)

# acer
summary(elev_mod_acer_poi)
summary(tci_mod_acer_poi)
summary(streamdist_mod_acer_poi)
summary(disturb_mod_acer_poi)
summary(beers_mod_acer_poi)


# poisson with elevation interaction only----
abies_poi = glm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers + elev * streamdist + elev * disturb, data = abies_subset, 
                family='poisson')
acer_poi = glm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers + elev * streamdist + elev * disturb, data = acer_subset, 
               family='poisson')

summary(abies_poi)
summary(acer_poi)

# pseudo r ^2
pseudo_r2 = function(glm_mod) {
  1 -  glm_mod$deviance / glm_mod$null.deviance
}
pseudo_r2(abies_poi)
pseudo_r2(acer_poi)

# comparing poisson with equivalent OLS----
summary(abies_poi)
summary(elev_interaction_mod_abiessubset)
pseudo_r2(abies_poi)

summary(acer_poi)
summary(elev_interaction_mod_acersubset)
pseudo_r2(acer_poi)

# comparing with anova
Anova(abies_poi, type=3)
Anova(elev_interaction_mod_abiessubset, type=3)
Anova(acer_poi, type=3)
Anova(elev_interaction_mod_acersubset, type=3)

# poisson might be better at estimating abies, but not acer based on the pseudo
# r squared values

reduced <- stepAIC(full_mod_abiessubset)

# model diagnostics with plot of model
# observed vs predicted plot lines 
# change order of levels to choose which is the intercept
# use contrast to visualize level code

# checking plots of the elevation interaction model and removing more outliers---
par(mfrow = c(2,2))
plot(elev_interaction_mod_abiessubset)
# removed a few more outliers
# take out 56, 109,53  from abies subset

plot(elev_interaction_mod_acersubset)
# take out 56 , 187 from acer

# semifinal models to use:----
# elev_interaction_mod_acersubset, elev_interaction_mod_abiessubset
# abies_poi (only includes elevation interaction)
# acer_poi (only includes elevation interaction)

# use stepAIC----
# to check the model and decide if any variables can be removed----
stepAIC(elev_interaction_mod_abiessubset)
# shows that it may make sense to remove elev:beers from the abies model
stepAIC(elev_interaction_mod_acersubset)
# shows that elevation and stream distance interaction may not be neccesarry

# retry with those taken out----
#OLS
abies_test_mod = lm(cover ~ elev + tci + streamdist + disturb + elev * tci + elev * streamdist + elev * disturb, data = abies_subset)
acer_test_mod = lm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers +  elev * disturb, data = acer_subset)
summary(abies_test_mod)
summary(acer_test_mod)
#POISSON
abies_test_poi = glm(cover ~ elev + tci + streamdist + disturb + elev * tci + elev * streamdist + elev * disturb, data = abies_subset, family = "poisson")
acer_test_poi = glm(cover ~ elev + tci + streamdist + disturb + beers + elev * tci + elev * beers +  elev * disturb, data = acer_subset, family = "poisson")
summary(abies_test_poi)
summary(acer_test_poi)
summary(elev_interaction_mod_abiessubset)

# the adjusted r squares show this did help
# after stepAIC abies
AIC(abies_test_mod)
AIC(abies_test_poi)
# original abies
AIC(elev_interaction_mod_abiessubset)
AIC(abies_poi)
# after stepAIC acer
AIC(acer_test_mod)
AIC(acer_test_poi)
# original acer
AIC(elev_interaction_mod_acersubset)
AIC(acer_poi)

# the AIC's also show this was helpful, so we can rename the test models as our
# final models to present... keep in mind these final models are all using the subsetted data without outliers
# the name was just shortened for easier access and readability. 

# final models----
final_abies_model <- abies_test_mod
final_acer_model <- acer_test_mod
final_acer_poisson <- acer_test_poi
final_abies_poisson <- abies_test_poi

# 2) and----
# 1) part: "do model diagnostics indicate any problems with violations of OLS assumptions?"----
# OLS versus poisson for each species
# new model to prove this point with the abies species.
abies_elev_poi = glm(cover ~ elev, data= abies_subset, family="poisson")
# also using elev_interaction_mod_abiessubset = lm(cover ~ elev, data= abies_subset)

# When comparing these plots it confirms suspicion that the Poisson distribution is a better fit
# for the abies species, and this is most likely due to the highly skewed elevation data. Since it can
# be seen in the original (cover ~ elevation) plot for abies, the data is highly concentrated around 0 
# coverage until an elevation of 1500, where the trend turns visually turns into a positive correlation between 
# elevation and abies cover. 

par(mfrow=c(1,1))
plot(cover ~ elev, data= abies_subset)

# This can be explained by the fact that the abies species are fraisir firs and they 
# are highly selective on their habitat condition, they must be further north, or high in elevation, in order to 
# grow. So, when the model lm(cover ~ elev, data= abies[_subset]), or any other model including elevation as a 
# variable describing cover will not be taking into account this highly skewed data as it is working under a normal 
# distribution. This can be seen with the hotpink link on the graph below. This is the residual line for cover ~ elev
# for the abies species under OLS conditions. 

plot(cover ~ elev, data= abies_subset)
abline(elev_interaction_mod_abiessubset, col="hot pink")

# The line trends downwards, which we know to be false as the data trend
# actually increases as elevation increases, after a specific point. However, it can be seen through the blue line on the 
# graph that when the model glm(cover ~ elev, data= abies[_subset], family="poisson") is plotted the regression line is 
# much more in line with what we would expect to predict the data. 

plot(cover ~ elev, data= abies_subset)
abline(elev_interaction_mod_abiessubset, col="hot pink")
abline(abies_elev_poi, col= "light blue")

# The intercept begins around the time the data points begin to increase, and it continues with an upward trend until 
# the end of the graph, just like the data points.

# However, for the acer species,
# when comparing these plots it confirms suspicion that the OLS distribution is a better fit
# for the acer species, the residuals vs leverage graph and scale- location are much more
# better distributed along the line, and more evenly distributed in space during the OLS regression
# rather than the poisson regression. 

par(mfrow= c(2,2))
plot(final_acer_model)
plot(final_acer_poisson)

# more proof that OLS is better for acer, and poisson is better for abies
pseudo_r2(final_acer_poisson)
pseudo_r2(final_abies_poisson)

# Compare these to the output adjusted r squared from normal 
summary(final_acer_model)
summary(final_abies_model)
# this further proves that OLS is better fitted for the Acer species as the pseudo r^2 value
# is less than the normal r^2 value, whereas the pseudo R^2 value for abies is higher than the OLS value

# So, it can be concluded that the best models to represent the cover for species Acer and Abies are:
# the OLS regression of Acer (final_acer_model) and the poisson regression of Abies (final_abies_poisson),
# as the skewed data in the abies species makes it less accurate when using OLS regression.

#' 4. (optional) Examine the behavior of the function stepAIC() using the exploratory 
#' models developed above. This is a very simple and not very robust machine learning 
#' stepwise algorithm that uses AIC to select a best model. By default it does a backward selection routine.

# repeat of the stepAIC used in the model selection above to answer question four.
stepAIC(elev_interaction_mod_acersubset)

# using the output from the stepAIC(elev_interaction_mod_acersubset) it can be used to describe the functionality of stepAIC().
# when used the stepwise function takes a backward approach by first looking at all of the variables included in a model taken in the function.
# By running through the function the first time it is identified that the overall AIC is 1310.9, however after the second output when the function
# is modified it is identified to have an AIC of 1309.36. It can be seen that this is equivalent to the AIC Value corresponding to the elev:streamdist 
# ,it can be concluded that the stepAIC function identifies the AIC value of the model IF the variable it is describing is dropped. From this output it 
# can be further concluded that these are the variables which can be removed to to improve the model in a simple way.So, using the Acer subset, the variables 
# elev:streamdist is the onlt variables that should be considered to be removed. However, in the stepAIC(elev_interaction_mod_abiessubset), shown below, it can be 
# determined that removing elev:beers and beers will both improve the proposed model. 

stepAIC(elev_interaction_mod_abiessubset)



