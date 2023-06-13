# Analysis script for Small Hive Beet control via anthranillic diamide
# Script by Dr Lewis J Bartlett & Ethan Hackmeyer
# lewis.bartlett@uga.edu


######################################################################

###############
getwd()
###############

library(afex)
library(emmeans)

###############

# # Larval SHB Mortality Analysis - Dousing

SHBLarvae.douse.mort <- read.csv(file='SHBLarvaeLawnDrenchDousing.csv', header=T, stringsAsFactors = F)

M1 <- glm(SHBLarvae.douse.mort$TD ~ SHBLarvae.douse.mort $Treatment,
          family = 'binomial')

anova(M1, test='Chisq')
# p = 0.00129. Significant whether a plate is treated or untreated.

SHBLarvae.douse.mort$Comb <- cbind(SHBLarvae.douse.mort$ND,SHBLarvae.douse.mort$NS)

M2 <- glm(SHBLarvae.douse.mort$Comb ~ SHBLarvae.douse.mort$Treatment,
          family = 'binomial')

anova(M2, test='Chisq')
# p < 0.001. Significant for an individual larvae if it was in a treated plate or not.

################

# # SHB Adult Exposure

SHBAdult.douse.mort <- read.csv(file='AdultSHBDousing.csv', header=T, stringsAsFactors = F)

# No analysis undertaken as all survived (see manuscript)

################

# # Adult SHB Mortality & Reproduction Analysis

AdultSHB.mort <- read.csv(file='SHBAdultFeedMortalityRepro.csv', header = T, stringsAsFactors = F)

AdultSHB.mort$CombD4 <- cbind(AdultSHB.mort$D4D, AdultSHB.mort$D4A)
D4Model <- mixed(CombD4 ~ Dose + (1 | CageID/Dose), family = 'binomial', method = "LRT", data = AdultSHB.mort)
nice(D4Model)

AdultSHB.mort$CombD6 <- cbind(AdultSHB.mort$D6D, AdultSHB.mort$D6A)
D6Model <- mixed(CombD6 ~ Dose + (1 | CageID/Dose), family = 'binomial', method = "LRT", data = AdultSHB.mort)
nice(D6Model)

AdultSHB.mort$CombD18 <- cbind(AdultSHB.mort$D18D, AdultSHB.mort$D18A)
D18Model <- mixed(CombD18 ~ Dose + (1 | CageID/Dose), family = 'binomial', method = "LRT", data = AdultSHB.mort)
nice(D18Model)

# AdultSHB.mort$CombD34 <- cbind(AdultSHB.mort$D34D, AdultSHB.mort$D34A)
# D34Model <- mixed(CombD34 ~ Dose + (1 | CageID/Dose), family = 'binomial', method = "LRT", data = AdultSHB.mort)
# nice(D34Model)

D18LarvaeModel <- glm(D18L ~ Dose, family = 'binomial' , data = AdultSHB.mort)
anova(D18LarvaeModel, test = 'Chisq')

# D34LarvaeModel <- glm(D34L ~ Dose, family = 'binomial' , data = AdultSHB.mort)
# anova(D34LarvaeModel, test = 'Chisq')

##########################

# Field Trial - Patty Experiment

PattyData <- read.csv(file = 'FieldTrialData.csv', header = T, stringsAsFactors = F)

Skrrt <- mixed(PattyMass ~ Treatment + (1 | Site), 
               data = PattyData)

nice(Skrrt)

emmeans(Skrrt, specs = 'Treatment')

# p = .113. Bees consumed more treated patty than control after 5 days, but it was not significant.

Apples <- mixed(D0Larvae ~ Treatment + (1 | Site),
                data = PattyData)
nice(Apples)

# Unsurprisingly, p < 0.001.

# Exact Test on larvae presence in patties after 14 days

ExactData <- data.frame(
  "Treated_yes" = c(0,10),
  "Treated_no" = c(10,0),
  row.names = c("Treatment", "Control"),
  stringsAsFacotrs = FALSE
)
colnames(ExactData) <- c("Larvae","No Larvae")

fisher.test(ExactData)

# Fisher test p-value p < .001. Treated pots significantly less likely to host larvae than control pots.

fisher.test(matrix(c(10,0,0,10), ncol = 2))


##########################

# # # Adult Honey Bee Mortality Analyses

# # Initial chlorantraniliprole test

HonBe.mort <- read.csv(file='HBAdultOralMortalityASSAY1.csv', header = T, stringsAsFactors = F)

# 24hr mortality analysis

#Not done as no mortality observed at all.

# 48Hr Mortality Analysis
HonBe.mort$Comb <- cbind(HonBe.mort$D48H, HonBe.mort$A48H)
HonBe.mort$Starvation <- HonBe.mort$Consumed48H
for( H in 1:NROW(HonBe.mort)){
  if(HonBe.mort$Starvation[H] == 2.00)
  {HonBe.mort$Starvation[H] = TRUE}
  else
  {HonBe.mort$Starvation[H] = FALSE}
}

MM.2 <- mixed(Comb ~ Dose + Starvation + Methanol + (1|Colony/Cage), family = 'binomial', 
            method = "LRT", 
            data = HonBe.mort)
# Dose p-value: .615. No significant effect of dose on bee mortality.

nice(MM.2)

# # Full second chlorantraniliprole test

# Adult Honeybee Mortality Analysis V.2 @ 48H

HonBe2.mort <- read.csv(file = 'HBAdultOralMortalityASSAY2.csv', header = T, stringsAsFactors = F)

HonBe2.mort$Comb48h <- cbind(HonBe2.mort$Dead48h, HonBe2.mort$Alive48h)

Honk <- mixed(Comb48h ~ Dose + Methanol + Dimethoate + (1|Colony/Cage), family ='binomial',
               method = 'LRT',
               data = HonBe2.mort)
nice(Honk)
# Dose p-value = .106 --> no significant effect of a.i. on adult honeybee mortality, though this p-value is much lower than that of the first trial
# Methanol p-value = .001 --> significant effect of methanol on mortality ?
# Dimethoate p-value < .001 --> (Obviously) dimethoate has a significant effect on adult honeybee mortality (duh)

###############################

# Flubendiamide assay

Flub.mort <- read.csv(file = 'HBAdultOralMortalityFLUBENDIAMIDE.csv', header = T, stringsAsFactors = F)

Flub.mort$A24H <- Flub.mort$TotalBees - Flub.mort$D24H
Flub.mort$A48H <- Flub.mort$TotalBees - Flub.mort$D48H

Flub.mort$Comb48h <- cbind(Flub.mort$D48H, Flub.mort$A48H)

# not having 

HonBe.F.mort.48H <-  mixed(Comb48h ~ AIConc + Methanol + Dimethoate + (1|Colony/Cage), family ='binomial',
                           method = 'LRT',
                           data = Flub.mort)

nice(HonBe.F.mort.48H)

###############################

### Graphing

# Dose mortality curve for Adult Honeybees and Small Hive Beetles

Honk2 <- glm(Comb48h ~ Dose + Methanol + Dimethoate, family = 'binomial', 
               data = HonBe2.mort)


par(mar = c(5,5,5,2))

plot(x= log10(1:300),
     y = (-10:289)*0.0035, 
     type='n', 
     xlab =  expression(paste('Log'[10],'(x+1)  Dose Chlorantraniliprole', sep='')),
     ylab = 'Proportion Dead', 
     main = 'Honeybee vs. SHB Mortality',
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)



GV1 <- seq(1,300, length.out = 500)

GV2 <- predict(Honk2, 
               newdata = data.frame(Dose = GV1, 
                                    Dimethoate = 0, 
                                    Methanol = 1),
                                    type = 'response')

points(x = log10(GV1), y = GV2, type = 'l', lwd = 4, col = 'blue3', pch = 20, cex = 1.75)

HonBe2.mort$PD48 <- HonBe2.mort$Dead48h/HonBe2.mort$InitialNumber

points(x = jitter(log10((HonBe2.mort$Dose[which(HonBe2.mort$Methanol == 1 & HonBe2.mort$Dimethoate == 0)])+1), factor = 0.5),
       y = HonBe2.mort$PD48[which(HonBe2.mort$Methanol == 1 & HonBe2.mort$Dimethoate == 0)],
       type = 'p', lwd = 2 ,
       col = 'blue3',
       pch = 20, cex=2.5)


AdultSHB.mort$D4Total <- AdultSHB.mort$D4D + AdultSHB.mort$D4A
AdultSHB.mort$D4PD <- AdultSHB.mort$D4D / AdultSHB.mort$D4Total
CombDeadSHBD4 <- cbind(AdultSHB.mort$D4D, AdultSHB.mort$D4A)

JenniferCoolidge <- glm(CombDeadSHBD4 ~ Dose, # + (1 | CageID/Dose), 
                        family = 'binomial', 
                        data = AdultSHB.mort)

Vodka <- predict(JenniferCoolidge, 
                 newdata = data.frame(Dose = GV1),
                 type = 'response')

points(x = log10(GV1), y = Vodka, type = 'l', lwd = 4, col = 'orange', pch = 20, cex = 1.75)
points(x = jitter(log10(AdultSHB.mort$Dose + 1), factor = 0.5),
       y = AdultSHB.mort$D4PD,
       type = 'p', lwd = 2 ,
       col = 'orange',
       pch = 20, cex=2.5)



### Convert to ug g-1 for both

# But dose mortality curve for Adult Honeybees and Small Hive Beetles

HonBe3.mort <- HonBe2.mort

HonBe3.mort$DoseT <- (HonBe3.mort$Dose)/1.23

Honk3 <- glm(Comb48h ~ DoseT + Methanol + Dimethoate, family = 'binomial', 
             data = HonBe3.mort)

par(mar = c(5,5,5,2))

plot(x= log10(1:300),
     y = (-10:289)*0.0035, 
     type='n', 
     xlab =  expression(paste('Log'[10],'(x+1)  DoseT Chlorantraniliprole', sep='')),
     ylab = 'Proportion Dead', 
     main = 'Honeybee vs. SHB Mortality',
     cex.main = 1.75, cex.lab = 1.635, cex.axis = 1.635)



GV1 <- seq(1,300, length.out = 500)

GV2 <- predict(Honk3, 
               newdata = data.frame(DoseT = GV1, 
                                    Dimethoate = 0, 
                                    Methanol = 1),
               type = 'response')

points(x = log10(GV1), y = GV2, type = 'l', lwd = 4, col = 'blue3', pch = 20, cex = 1.75)

HonBe3.mort$PD48 <- HonBe3.mort$Dead48h/HonBe3.mort$InitialNumber

points(x = jitter(log10((HonBe3.mort$DoseT[which(HonBe3.mort$Methanol == 1 & HonBe3.mort$Dimethoate == 0)])+1), factor = 0.5),
       y = HonBe3.mort$PD48[which(HonBe3.mort$Methanol == 1 & HonBe3.mort$Dimethoate == 0)],
       type = 'p', lwd = 2 ,
       col = 'blue3',
       pch = 20, cex=2.5)


AdultSHB.mort$D4Total <- AdultSHB.mort$D4D + AdultSHB.mort$D4A
AdultSHB.mort$D4PD <- AdultSHB.mort$D4D / AdultSHB.mort$D4Total
CombDeadSHBD4 <- cbind(AdultSHB.mort$D4D, AdultSHB.mort$D4A)

JenniferCoolidge <- glm(CombDeadSHBD4 ~ Dose, # + (1 | CageID/DoseT), 
                        family = 'binomial', 
                        data = AdultSHB.mort)

Vodka <- predict(JenniferCoolidge, 
                 newdata = data.frame(Dose = GV1),
                 type = 'response')

points(x = log10(GV1), y = Vodka, type = 'l', lwd = 4, col = 'orange', pch = 20, cex = 1.75)
points(x = jitter(log10(AdultSHB.mort$Dose + 1), factor = 0.5),
       y = AdultSHB.mort$D4PD,
       type = 'p', lwd = 2 ,
       col = 'orange',
       pch = 20, cex=2.5)








