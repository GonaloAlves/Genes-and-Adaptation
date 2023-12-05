#Libraries that were used in script

library(readr)
library(fitdistrplus)
library(glmmTMB)
library(ggplot2)
library(car)
library(emmeans)



#Import and read table (Portuguese decimal format)
gadata = read.table("Genes e Adaptação - Folha1.csv", header = T, sep = ";" ,stringsAsFactors = TRUE, dec = ",", na.strings =  "NA")
gadata

#Summary of the data
summary(gadata) 
str(gadata) 



#var independent/fix factors = Popstruct e Temp
#var dependent = SR, HR and Fecundity
#random factors = block and Group


#check right distrubution (not conclusive)

#we add as.numeric and na.omit to omit the NA's 
descdist(as.numeric(na.omit(gadata$HatchingRate)), boot = 100, discrete = TRUE) #cannot tell difference in plot
descdist(as.numeric(na.omit(ola)), boot = 100, discrete = TRUE) #cannot tell difference in plot
descdist(as.numeric(na.omit(gadata$Nr_Eggs)), boot = 100, discrete = TRUE) #cannot tell difference in plot


#Model for SexRatio
  
lm_SR = glmmTMB(cbind(Nr_Female, Nr_Male) ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = "binomial") #NA cant figure that out 

summary(lm_SR)

Anova(lm_SR)
emmeans(lm_SR,specs =pairwise~Temp:PopStruct ,type= "response")


#Model for HatchingRate

lm_HR = glmmTMB(cbind(Nr_Eggs, Nr_Juv) ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = "binomial") #NA cant figure that out 

summary(lm_HR)

Anova(lm_HR)
emmeans(lm_HR,specs =pairwise~Temp:PopStruct ,type= "response")

#model for Fecundity 

lm_FP = glmmTMB(Fecundity ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = poisson)

lm_F1 = glmmTMB(Fecundity ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = nbinom1)

lm_F2 = glmmTMB(Fecundity ~ PopStruct*Temp + (1|Block) + (1|Group), data = gadata, family = nbinom2)

summary(lm_F1)
AIC(lm_FP,lm_F1,lm_F2) #Comparing all models for se who fits better (lm_F1 - negative binominal 1)

Anova(lm_F1) 
emmeans(lm_F1,specs =pairwise~Temp:PopStruct ,type= "response")


#Plots

#fecundity
ggplot(data = gadata, aes(x = factor(Temp), y = Fecundity, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
 # facet_wrap(~ factor(Temp))+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Fecundity")+
  xlab("Temperature") +
  ggtitle("Fecundity in different populations and treatments") 


#hatching rate
ggplot(data = gadata, aes(x = factor(Temp), y = HatchingRate, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
 # facet_wrap(~ PopStruct)+
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))+
  ylab("Hatching Rate")+
  xlab("Temperature") +
  ggtitle("Hatching Rate in different populations and treatments") 


#sexratio
ggplot(data = gadata, aes(x = factor(Temp), y = SR, fill = PopStruct)) +
  geom_boxplot() +
  stat_summary(fun.y = "mean", geom = "point", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  stat_summary(fun.y = "mean", geom = "line", aes(group = PopStruct), position = position_dodge(width = 0.75)) +
  #facet_wrap(~ PopStruct) +
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  ylab("Sex Ratio") +
  xlab("Temperature") +
  ggtitle("Sex Ratio in different populations and treatments") 


#Try of new graph (Not used in paper)
ggplot(data = gadata, aes(x = factor(Total_F), y = Dead_F)) +
  geom_point(aes(color = Dead_F), position = position_jitter(width = 0.2, height = 0.2)) +
  stat_summary(fun.y = "mean", geom = "point", color = "red", shape = 18, size = 3, position = position_dodge(width = 0.75)) +
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~ PopStruct) +
  theme(legend.position = "right", axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10, face = "bold"), axis.title.y = element_text(size = 10, face = "bold"), plot.title = element_text(hjust = 0.5, size = 10, face = "bold")) +
  ylab("Dead Females") +
  xlab("Total Females") +
  ggtitle("Females and Dead")



