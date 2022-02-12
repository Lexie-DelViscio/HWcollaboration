## Multivariate Assignment 
library(vegan) 
data(dune)
data(dune.env)
head(dune)
head(dune.env)
names(dune) # species 
names(dune.env) # variable

par(mfrow=(c(1,2)))
dune_mds = metaMDS(dune)

plot(dune_mds, type='n')
text(dune_mds, 'sp', cex=.5)
# generate vector of colors 
color_vect = rev(terrain.colors(6))[-1]
points(dune_mds, 'sites', pch=19, 
       col=color_vect[dune.env$Moisture])
legend('topright', paste("Moisture =", 1:5, sep=''), 
       col=color_vect, pch=19)
dune_mds
str(dune_mds)

# axes kinda uninformative
# the first axis is broader than the second
# first variable matches up with the first axis pretty well

# Species Eleopalu and Ranuflam are usually found close to one another in an area with alot of higher moisture points,
# whereas Rumeacet is concentrated around lots of low moisture points

# it seems like species on the left tend to line up more with lower levels of moisture, 
# species towards the right of the plot line up more with higher levels pf moisture.
# showing there may be a strong correlation between moisture and where each specific species is found 

# it computes between every sample a distance matrix, ignores where both are missing since it is a non euclidean method 

# Carry out a direct ordination using CCA in order to test any potential hypotheses that 
# you developed after examining the MDS plot. Specifically, carry out a test of the entire model 
# (i.e., including all constrained axes) and also carry out tests at the scale of individual 
# explanatory variables you included in your model if you included more than one variable. 
# Plot your results.

# Canonical Correspondence Analysis (CCA)
# Let's carry out a Canonical Correspondence Analysis (CCA) as well. CCA is appropriate 
# for modeling unimodal or hump-shaped responses to explanatory variables (rather than linear as with RDA).

cca_dune <- cca(dune ~ ., data=dune.env)
RsquareAdj(cca_dune, 100)
anova(cca_dune, permutations = 999)
anova(cca_dune, by='margin', permutations = 999)

plot(cca_dune, type='n', scaling=1)
orditorp(cca_dune, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune, display='cn', col='red')


# just moisture
cca_dune_moisture <- cca(dune ~ dune.env$Moisture)
RsquareAdj(cca_dune_moisture, 100)
anova(cca_dune_moisture, permutations = 999)
anova(cca_dune_moisture, by='margin', permutations = 999)

plot(cca_dune_moisture, type='n', scaling=1)
orditorp(cca_dune_moisture, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_moisture, display='cn', col='red')
dune.env$Moisture
levels(dune.env$Moisture)



#Just A1
cca_dune_A1 <- cca(dune ~ dune.env$A1)
RsquareAdj(cca_dune_A1, 100)
anova(cca_dune_A1, permutations = 999)
anova(cca_dune_A1, by='margin', permutations = 999)

plot(cca_dune_A1, type='n', scaling=1)
orditorp(cca_dune_A1, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_A1, display='bp', col='red')

#Just manure
cca_dune_Manure <- cca(dune ~ dune.env$Manure)
RsquareAdj(cca_dune_Manure, 100)
anova(cca_dune_Manure, permutations = 999)
anova(cca_dune_Manure, by='margin', permutations = 999)

plot(cca_dune_Manure, type='n', scaling=1)
orditorp(cca_dune_Manure, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_Manure, display='cn', col='red')

#Just Management
cca_dune_Management <- cca(dune ~ dune.env$Management)
RsquareAdj(cca_dune_Management, 100)
anova(cca_dune_Management, permutations = 999)
anova(cca_dune_Management, by='margin', permutations = 999)

plot(cca_dune_Management, type='n', scaling=1)
orditorp(cca_dune_Management, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_Management, display='cn', col='red')

#Just Use
cca_dune_Use <- cca(dune ~ dune.env$Use)
RsquareAdj(cca_dune_Use, 100)
anova(cca_dune_Use, permutations = 999)
anova(cca_dune_Use, by='margin', permutations = 999)

plot(cca_dune_Use, type='n', scaling=1)
orditorp(cca_dune_Use, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_Use, display='cn', col='red')


#vector moisture 
vector_moisture <- as.numeric(as.character(dune.env$Moisture))
vector_moisture
dune.env$Moisture
vector_moisture_cca <- cca(dune~vector_moisture)
RsquareAdj(vector_moisture_cca, 100)
anova(vector_moisture_cca, permutations = 999)
anova(vector_moisture_cca, by='margin', permutations = 999)

plot(vector_moisture_cca, type='n', scaling=1)
orditorp(vector_moisture_cca, display='sp', cex=0.5, scaling=1, col='blue')
text(vector_moisture_cca, display='bp', col='red')

# vector manure
vector_manure <- as.numeric(as.character(dune.env$Manure))
vector_manure_cca <- cca(dune~vector_manure)
RsquareAdj(vector_manure_cca, 100)
anova(vector_manure_cca, permutations = 999)
anova(vector_manure_cca, by='margin', permutations = 999)

plot(vector_manure_cca, type='n', scaling=1)
orditorp(vector_manure_cca, display='sp', cex=0.5, scaling=1, col='blue')
text(vector_manure_cca, display='bp', col='red')

#new dune cca
cca_dune_new <- cca(dune ~ A1 + Use + Management + vector_moisture + vector_manure, data=dune.env)
RsquareAdj(cca_dune_new, 100)
anova(cca_dune_new, permutations = 999)
anova(cca_dune_new, by='margin', permutations = 999)

plot(cca_dune_new, type='n', scaling=1)
orditorp(cca_dune_new, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_new, display='cn', col='red')


par(mfrow=(c(1,2)))
plot(dune_mds, type='n')
text(dune_mds, 'sp', cex=.5)
# generate vector of colors 
color_vect <- rev(terrain.colors(6))[-1]
points(dune_mds, 'sites', pch=19, 
       col=color_vect[dune.env$Moisture])
legend('topright', paste("Moisture =", 1:5, sep=''), 
       col=color_vect, pch=19)
dune_mds
str(dune_mds)

plot(cca_dune_new, type='n', scaling=1)
orditorp(cca_dune_new, display='sp', cex=0.5, scaling=1, col='blue')
text(cca_dune_new, display='cn', col='red')



RsquareAdj(cca_dune_new, 100)
RsquareAdj(cca_dune_A1, 100)
RsquareAdj(cca_dune_Management, 100)
RsquareAdj(cca_dune_Use, 100)
RsquareAdj(cca_dune_moisture, 100)
RsquareAdj(cca_dune_Manure, 100)
RsquareAdj(vector_manure_cca, 100)
RsquareAdj(vector_moisture_cca, 100)

# should try removing use from model, only 1.6 rsquared adj
# chnaging to vectors, a;though better to visualize does decrease the r square, but increase adj r squre,
# so overall probably a better model even though individuallly the vecotrs have lower r squres than their counterparts
# they help to decrease the overfitting of the model due to an excess of variables

without_use_cca <- update(cca_dune_new, . ~ . - Use)
par(mfrow=c(1,1))
RsquareAdj(without_use_cca, 100)
anova(without_use_cca, permutations = 999)
anova(without_use_cca, by='margin', permutations = 999)

plot(without_use_cca, type='n', scaling=1)
orditorp(without_use_cca, display='sp', cex=0.5, scaling=1, col='blue')
text(without_use_cca, display='cn', col='red')

# does increase the adj r but only by 1%, overall a very similar looking model this 


# manure plot like moisture nmds
plot(dune_mds, type='n')
text(dune_mds, 'sp', cex=.5)
# generate vector of colors 
color_vect = rev(heat.colors(6))[-1]
points(dune_mds, 'sites', pch=19, 
       col=color_vect[dune.env$Manure])
legend('topright', paste("Manure =", 1:5, sep=''), 
       col=color_vect, pch=19)

# this further shows the correlation between the cca and nmds ordination results. 

dune.env$Manure
?terrain.colors
