rm(list=ls())
library(car)
library(ggplot2)
library(ggbeeswarm)
library(gridExtra)

t<-read.csv('trials.csv',h=T)
t$MaleID<-factor(t$MaleID)
t$FemID<-factor(t$FemID)
t$CourtYN<-factor(t$CourtYN)
levels(t$CourtYN)<-c('NoCourtship','Courtship')

#get rid of trials where crickets did not interact
t1<-subset(t,DNI!=1)
#remove repeated trials
t2<-t1[!duplicated(t1$MaleID),]
t2<-t2[!duplicated(t2$FemID),]
t2<-subset(t2,MaleMorph!="NA")

t2[t2$MaleMorph=="wtFw",]$MaleMorph<-"Fw"
t2[t2$MaleMorph=="CwFw",]$MaleMorph<-"Fw"
t2$MaleMorph<-factor(t2$MaleMorph)


library(car)
glm.mp<-glm(MountYN~CourtYN*MaleMorph,data=t2,family='binomial')
Anova(glm.mp,type="III")

levels(t2$MaleMorph)<-c("CwNw","Fw (Wt & Cw)","WtNw")

lines<-ggplot(t2,aes(x=CourtYN,y=MountYN,fill=MaleMorph))+
  geom_quasirandom(shape=21,alpha=0.75,size=2)+
  stat_summary(fun.data=mean_se,aes(group=MaleMorph))+
  stat_summary(geom='line',aes(group=factor(MaleMorph)),
               colour='black',linetype='dashed')+
  ylab('Female mounting rate')+
  xlab('')+theme_bw()+scale_fill_manual(values=c("#355C7D","#d18015","#C06C84"))+
  theme(legend.position='none',legend.title=element_blank(),
        panel.grid=element_blank())+facet_grid(.~factor(MaleMorph))
lines
ggsave('mate_pref.png',plot=lines,height=3,width=6.5,dpi=600,units='in')
