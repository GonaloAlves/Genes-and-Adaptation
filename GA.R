#Libarys

library(readr)
library(fitdistrplus)
library(glmmTMB)
library(ggplot2)
library(car)
library(lme4)
library(lmerTest)
library(emmeans)

#Import and read table
gadata = read.table("Genes e Adaptação - Folha1.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",", na.strings =  "NA")
gadata

#Summary
summary(gadata) #still missing data
str(gadata) #Locations some are coded in different ways


#distribution of data

hist(gadata$Dead_F) #Mostly in 0 but half of 0 are 1
hist(gadata$Total_F) #2 and 10 are common as expect
hist(gadata$Nr_Eggs) #so many 0s
hist(gadata$Nr_Juv) #same as eggs
hist(gadata$HatchingRate) #not normal distrubuted
hist(gadata$Nr_Female) #not much either
hist(gadata$Nr_Male) #even less than female
hist(gadata$SR) #not normal distributed to


#var independentes/fix factors = popstruct e temp
#var dependentes = sr e hr
#random factors = block, Petridish, location? and dead_F? 

#models

#check right distrubution

hist(gadata$HatchingRate)
hist(gadata$SR)
  

#we add as.numeric and na.omit to omit the NA's 
descdist(as.numeric(na.omit(gadata$HatchingRate)), boot = 100, discrete = TRUE) #cannot tell difference in plot
descdist(as.numeric(na.omit(ola)), boot = 100, discrete = TRUE) #cannot tell difference in plot

descdist(as.numeric(na.omit(gadata$Nr_Eggs)), boot = 100, discrete = TRUE) #cannot tell difference in plot



lm_HR = glmmTMB(cbind(Nr_Juv, Nr_Eggs) ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = "binomial")

summary(lm_HR)

boxplot(gadata$HatchingRate ~ gadata$PopStruct*gadata$Temp)

Anova(lm_HR)
emmeans(lm_HR,specs =pairwise~Temp:PopStruct ,type= "response")


#lm_HR_nb <- glmmTMB(HatchingRate ~ PopStruct*Temp + (1|PetriDish) + (1|Block) + (1|Group), data = gadata, family=nbinom1, ziformula = ~1) #NA cant figure that out

#if we do the poisson in the same way as nb its INF


ola = cbind(gadata$Nr_Female, gadata$Nr_Male)
ola

#Model for SexRation
  
  lm_SR = glmmTMB(cbind(Nr_Female, Nr_Male) ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = "binomial") #NA cant figure that out 

summary(lm_SR)

Anova(lm_SR)
emmeans(lm_SR,specs =pairwise~Temp:PopStruct ,type= "response")


#Model for HatchingRate

lm_SR = glmmTMB(cbind(Nr_Eggs, Nr_Juv) ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = "binomial") #NA cant figure that out 

summary(lm_SR)

Anova(lm_SR)
emmeans(lm_SR,specs =pairwise~Temp:PopStruct ,type= "response")

#model for fecundity 

lm_F = glmmTMB(Nr_Eggs ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = nbinom1)

summary(lm_F)

Anova(lm_F)
emmeans(lm_F,specs =pairwise~Temp:PopStruct ,type= "response")


#plots

#fecundity
ggplot(data = gadata, aes(x = factor(Temp), y = Nr_Eggs, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
 # facet_wrap(~ factor(Temp))+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Number of Eggs")+
  xlab("Temperature") +
  ggtitle("Fecundity in different Population") 

#hatching rate

boxplot(gadata$HatchingRate ~ gadata$PopStruct*gadata$Temp) #between temps is clearly different, but not in the Structure

ggplot(data = gadata, aes(x = factor(Temp), y = HatchingRate, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
 # facet_wrap(~ PopStruct)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Hatching Rate")+
  xlab("Temperature") +
  ggtitle("Hatching Rate in different Population") 

#sexratio
ggplot(data = gadata, aes(x = factor(Temp), y = SR, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
 # facet_wrap(~ PopStruct)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Sex Ratio")+
  xlab("Temperature") +
  ggtitle("Sex Ratio in different Population") 



boxplot(gadata$SR ~ gadata$PopStruct*gadata$Temp) #nothing is relevant

meansr = mean(gadata$SR)

boxplot(gadata$Nr_Eggs ~ gadata$PopStruct*gadata$Temp) #com addpoints

#sr = cbind binominal

#hr = cbind binominal

#fecundity feito


