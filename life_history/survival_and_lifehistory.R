###
### original analysis by Jessica Bainbridge, condensed here by JGR
###

rm(list=ls())
library(ggplot2)
library(dplyr)
library(lmodel2)
library(lme4)
library(lmerTest)
library(ggeffects)
library(ggplot2)
library(khroma)

crix <- read.table("crixmass.csv", header=T, sep=",", fileEncoding = 'UTF-8-BOM')
# convert mass from g to mg
crix$mass <- 1000*crix$mass
# add wing_shape variable
crix$wing_shape <- crix$wing
crix$wing_shape <- as.character(crix$wing_shape)
# add wing_veins variable
crix$wing_veins <- crix$wing
crix$wing_veins <- as.character(crix$wing_veins)
crix$wing_shape[crix$wing_shape == "Cw, Nw"] <- "Cw"
crix$wing_shape[crix$wing_shape == "Cw, Fw"] <- "Cw"
crix$wing_shape[crix$wing_shape == "wt, Nw"] <- "Wt"
crix$wing_shape[crix$wing_shape == "wt, Fw"] <- "Wt"
crix$wing_shape[crix$wing_shape == "wt"] <- "Wt"
crix$wing_shape <- as.factor(crix$wing_shape)

crix$wing_veins[crix$wing_veins == "Cw, Nw"] <- "Nw"
crix$wing_veins[crix$wing_veins == "Cw, Fw"] <- "Fw"
crix$wing_veins[crix$wing_veins == "wt, Nw"] <- "Nw"
crix$wing_veins[crix$wing_veins == "wt, Fw"] <- "Fw"
crix$wing_veins <- as.factor(crix$wing_veins)

# Subset mass 14 days post eclosion
crix1 <- subset(crix, days==14)
# Subset to exclude individuals with lost limbs
crix1 <- subset(crix1, limbs == "")
fe <- subset(crix1, sex=="F")
ma <- subset(crix1, sex=="M")

# 2. Fit a SMA line to the mass length data on a natural log-log scale
library(lmodel2)
lmodel2(log(mass)~log(mpl), data = fe)
fbSMA = 2.713439
lmodel2(log(mass)~log(mpl), data = ma)
mbSMA = 2.194112 

# 3. Estimate the mean length
fx0 = mean(fe$mpl)
fx0 # 3.415682
mx0 = mean(ma$mpl)
mx0 # 3.634235

# 4. Calculate SMI for each sex
fSMI = fe$mass*(fx0/fe$mpl)^fbSMA
fe$SMI <- fSMI
mSMI = ma$mass*(mx0/ma$mpl)^mbSMA
ma$SMI <- mSMI
crix1 <- rbind(fe, ma)

# Calculate SMI for all weeks
# exclude individuals with lost limbs
crixl <- subset(crix, limbs == "")
fel <- subset(crixl,sex == "F" )
mal <- subset(crixl, sex == "M")
# calculate SMI using week 2 as formula
fSMI = fel$mass*(fx0/fel$mpl)^fbSMA
fel$SMI <- fSMI
mSMI = mal$mass*(mx0/mal$mpl)^mbSMA
mal$SMI <- mSMI
crix1 <- rbind(fel, mal)
all<-rbind(fel,mal)
fe <- subset(crix1, sex=="F")
ma <- subset(crix1, sex=="M")

#female mass
f.mass.lm <- lmer(mass ~ scale(days) + scale(I(days^2)) + wing_shape + scale(days)*wing_shape + scale(I(days^2))*wing_shape +
                (1 | id),
              REML = T, data = fe)
Anova(f.mass.lm,type="III")
qqnorm(resid(f.mass.lm))
qqline(resid(f.mass.lm))

f.smi.lm <- lmer(SMI ~ scale(days) + scale(I(days^2)) + wing_shape + scale(days)*wing_shape + scale(I(days^2))*wing_shape +
                    (1 | stock/id),
                  REML = T, data = fe)
Anova(f.smi.lm,type="III")
qqnorm(resid(f.smi.lm))
qqline(resid(f.smi.lm))

#female pl
f.pl.lm <- lmer((mpl) ~ wing_shape +
                    (1 | stock),
                  REML = T, data = fe[fe$days==14,])
qqnorm(resid(f.pl.lm))
qqline(resid(f.pl.lm))
Anova(f.pl.lm,type="II")

#male mass
m.mass.lm <- lmer(mass ~  wing_shape + wing_veins +
                    (1 | stock),
                  REML = T, data = ma[ma$days==14,])
Anova(m.mass.lm,type="II")
qqnorm(resid(m.mass.lm))
qqline(resid(m.mass.lm))

m.smi.lm <- lmer(SMI ~  wing_shape + wing_veins +
                    (1 | stock),
                  REML = T, data = ma[ma$days==14,])
Anova(m.smi.lm,type="II")
qqnorm(resid(m.smi.lm))
qqline(resid(m.smi.lm))

m.pl.lm <- lmer(mpl ~ wing_shape + wing_veins +
                  (1 | stock),
                REML = T, data = ma[ma$days==14,])
qqnorm(resid(m.pl.lm))
qqline(resid(m.pl.lm))
Anova(m.pl.lm,type="II")

#plot female/male mass
predf <- ggpredict(f.mass.lm, terms=c("days [all]", "wing_shape"),  type = "random")
f.mass<-ggplot(predf, aes(x, predicted, colour = group)) +
  #geom_point(data=fe,aes(x=days,y=mass,colour=wing_shape),size=0.2,alpha=0.5)+
  geom_line(size=0.8) + 
  geom_ribbon(aes(ymin = (predicted - std.error), ymax = (predicted + std.error), x = x, fill = group), linetype = 0, alpha = 0.3) +
  scale_x_continuous(breaks = seq(0,100, by = 7)) +
  scale_y_continuous(breaks = seq(0,1000, by = 100)) +
  guides(fill="none") +
  labs(x = "Time (Days)", y = "Mass (mg)") +  theme(aspect.ratio = 0.8) + 
  theme_bw()+scale_colour_manual(values=c('#355C7D','#C06C84'))+
  scale_fill_manual(values=c('#355C7D','#C06C84'))+
  theme(panel.grid=element_blank(),legend.position='top',legend.title=element_blank())

m.mass<-ggplot(ma[ma$days==14,], aes(x=wing_shape, y=mass, colour = wing_shape)) +
  geom_quasirandom(alpha=0.5,size=0.5)+
  stat_summary()+
  labs(x = "Wing_shape", y = "Mass (mg)") +  theme(aspect.ratio = 0.8) + 
  theme_bw() + 
  scale_y_continuous(breaks = seq(0,1000, by = 100)) +
  scale_colour_manual(values=c('#355C7D','#C06C84'))+
  scale_fill_manual(values=c('#355C7D','#C06C84'))+
  theme(panel.grid=element_blank(),legend.position='top',legend.title=element_blank())

####SURVIVAL
library(survival)
library(ggplot2)
library(survminer)
library(gtsummary)
# read in data
crix <- read.table("ProjectData.csv", header=T, fill = T, sep=",", fileEncoding = 'UTF-8-BOM')
crix <- crix[ which( crix$status == 1 | crix$status == 0) , ]

# remove crickets wiht lost limbs
crix <- subset(crix, limb == "")
## Calculate SMI for mass 1 (at eclosion)
crix$mass1 <- 1000*crix$mass1
fe <- subset(crix, sex=="F")
ma <- subset(crix, sex=="M")

fbSMA = 2.713439
mbSMA = 2.194112 
fx0 = 3.415682
mx0 = 3.634235

fSMI = fe$mass1*(fx0/fe$mpl)^fbSMA
fSMI
fe$SMI <- fSMI

mSMI = ma$mass1*(mx0/ma$mpl)^mbSMA
mSMI
ma$SMI <- mSMI

fe$wing_shape <- fe$wing
fe$wing_shape <- as.character(fe$wing_shape)
fe$wing_shape[fe$wing_shape == "Cw, Nw"] <- "Cw"
fe$wing_shape[fe$wing_shape == "Cw, Fw"] <- "Cw"
fe$wing_shape[fe$wing_shape == "wt, Nw"] <- "Wt"
fe$wing_shape[fe$wing_shape == "wt, Fw"] <- "Wt"
fe$wing_shape <- as.factor(fe$wing_shape)
fe$wing_veins <- fe$wing

mod1b <- coxph(Surv(survt, status == 1) ~ wing + SMI  + cluster(stock), data = fe)
summary(mod1b)

zp <- cox.zph(mod1b, transform=rank)
zp
# prop hazards assumption satisfied
mod1 <- mod1b

# make table
fe_tbl <- tbl_regression(mod1, exponentiate = TRUE) %>%
  add_n() %>%
  modify_header(label = "**Variable**")
fe_tbl

## Plotting cox-adjusted survival curves
# Create the new data  
fe_df <- with(fe,
              data.frame(wing = c("wt", "Cw"), 
                         SMI = rep(mean(SMI, na.rm = TRUE), 2)
              ))
fe_df

# Survival curves

fit <- survfit(mod1, newdata = fe_df)
f.curve<-ggsurvplot(fit = fit, data=fe, 
                    ylab = "Survival Probability", xlab= "Time (Days)",
                    legend.labs=c("F.Wt", "F.Cw"),
                    conf.int = TRUE,
                    conf.int.alpha=c(0.1),
                    break.x.by = 7,xlim=c(0,88),
                    legend.title=element_blank(),
                    palette = c('#C06C84','#355C7D'),
                    ggtheme= theme_bw() + 
                      theme(panel.grid = element_blank()))


#males
ma$wing_shape <- ma$wing
ma$wing_shape <- as.character(ma$wing_shape)
ma$wing_shape[ma$wing_shape == "Cw, Nw"] <- "Cw"
ma$wing_shape[ma$wing_shape == "Cw, Fw"] <- "Cw"
ma$wing_shape[ma$wing_shape == "wt, Nw"] <- "Wt"
ma$wing_shape[ma$wing_shape == "wt, Fw"] <- "Wt"
ma$wing_shape <- as.factor(ma$wing_shape)

# add wing_veins variable
ma$wing_veins <- ma$wing
ma$wing_veins <- as.character(ma$wing_veins)
ma$wing_veins[ma$wing_veins == "Cw, Nw"] <- "Nw"
ma$wing_veins[ma$wing_veins == "Cw, Fw"] <- "Fw"
ma$wing_veins[ma$wing_veins == "wt, Nw"] <- "Nw"
ma$wing_veins[ma$wing_veins == "wt, Fw"] <- "Fw"
ma$wing_veins <- as.factor(ma$wing_veins)
ma$wing<-ma$wing_shape

mod2b <- coxph(Surv(survt, status == 1) ~ wing_shape + wing_veins + SMI +
                 wing_shape*wing_veins + cluster(stock), data = ma)
summary(mod2b)
cox.zph(mod2b, transform=rank)
# assumption is met

# make table
ma_tbl <- tbl_regression(mod2b, exponentiate = TRUE) %>%
  add_n() %>%
  modify_header(label = "**Variable**")
ma_tbl

## Plotting cox-adjusted survival curves
# Create the new data  
ma_df <- with(ma,
              data.frame(wing_shape = c("Wt", "Cw"),
                         #wing_veins = c("Nw", "Fw", "Nw", "Fw"),
                         SMI = rep(mean(SMI, na.rm = TRUE), 2)
              ))
ma_df
# Survival curves
fitb <- survfit(mod2b, newdata = ma_df)
head(fitb)
m.curve<-ggsurvplot(fit = fitb, data=ma, 
                    ylab = "Survival Probability", xlab= "Time (Days)",
                    #legend.labs=c("wt, Nw", "wt, Fw", "Cw, Nw", "Cw, Fw"),
                    legend.labs=c("M.Wt","M.Cw"),
                    conf.int=TRUE,xlim=c(0,88),
                    palette = c('#C06C84','#355C7D'),
                    break.x.by = 7,
                    legend.title=element_blank(),legend.position='none',
                    ggtheme= theme_bw() + 
                      theme(panel.grid = element_blank()))



