library(ggplot2)
library(lme4)
library(car)
library(MuMIn)
library(gridExtra)
library(dplyr)

f1<-read.csv('Crosses_offspring.csv',h=T)
f1<-subset(f1,Mother!="NA")
f1$grp<-"f0"
f1[f1$MotherCw==0 & f1$FatherCw==0,]$grp<-"F.wt x M.wt"
f1[f1$MotherCw>0 & f1$FatherCw==0,]$grp<-"F.cw x M.wt"
f1[f1$MotherCw>0 & f1$FatherCw>0,]$grp<-"F.cw x M.cw"
f1[f1$MotherCw==0 & f1$FatherCw>0,]$grp<-"F.wt x M.cw"
f1$grp<-factor(f1$grp)

f1$Mother<-factor(f1$Mother)
f1$Father<-factor(f1$Father)
f1$f1.cwyn<-0
f1[f1$F1Cw>0,]$f1.cwyn<-1
f1$Mother.cwyn<-0
f1[f1$MotherCw>0,]$Mother.cwyn<-1
f1$Father.cwyn<-0
f1[f1$FatherCw>0,]$Father.cwyn<-1

f1sum<-data.frame(
  f1 %>%
    group_by(F1Sex,MotherCw,FatherCw,Mother,Father,grp,
             Mother.cwyn,Father.cwyn,MotherFw,FatherFw,density) %>%
    summarise(meancw.yn=mean(f1.cwyn,na.rm=TRUE),meanCw=mean(F1Cw,na.rm=TRUE), n = n()))

glm.cw.bin<-lmer(meancw.yn~1+(1|Father/Mother),data=f1sum)
Anova(glm.cw.bin)
r.squaredGLMM(glm.cw.bin)
qqnorm(resid(glm.cw.bin))
qqline(resid(glm.cw.bin))
plot(glm.cw.bin)

glm.cw.lin<-lmer(meanCw~1+(1|Father/Mother),data=f1sum)
Anova(glm.cw.lin)
r.squaredGLMM(glm.cw.lin)
qqnorm(resid(glm.cw.lin))
qqline(resid(glm.cw.lin))
plot(glm.cw.lin)

glm.cw.bin.int<-lmer(meancw.yn~F1Sex+density+FatherCw*MotherCw+(1|Father/Mother),data=f1sum)
Anova(glm.cw.bin.int,type="III")
r.squaredGLMM(glm.cw.bin.int)
qqnorm(resid(glm.cw.bin.int))
qqline(resid(glm.cw.bin.int))
plot(glm.cw.bin.int)

ggplot(f1sum,aes(x=Father,y=meancw.yn,fill=F1Sex))+
  geom_point(shape=21)+facet_grid(F1Sex~grp,scales='free_x')+
  scale_fill_manual(values=c("#fa8490","#9fd4fc"))+
  theme_bw()+theme(axis.text.x=element_blank(),panel.grid.major.x=element_blank())
