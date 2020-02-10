library(grid)
#library(gridExtra)
library(plyr)
library(ggplot2)

#set your working directory
setwd("/home/samuel/Documents/Hote-parasite/HPV/PTRS/Evolution_oncoviruses/Analysis")


###load files
#sumstats<-read.csv(file="sumstats_mean_retrovir1.csv",header = TRUE)
#sumstats_median<-read.csv(file="sumstats_median_retrovir1.csv",header = TRUE)

sumstats1<-read.csv(file="sumstats_mean_retrovir1.csv",header = TRUE)
sumstats1$virus<-"retro-like"
sumstats1$rel_trans<-sumstats1$trans/max(sumstats1$trans)
sumstats_median1<-read.csv(file="sumstats_median_retrovir1.csv",header = TRUE)
sumstats_median1$virus<-"retro-like"
sumstats_median1$rel_trans<-sumstats_median1$trans/max(sumstats_median1$trans)

sumstats2<-read.csv(file="sumstats_mean_sDNA_div1.csv",header = TRUE)
sumstats_median2<-read.csv(file="sumstats_median_sDNA_div1.csv",header = TRUE)
sumstats2$virus<-"small-DNA"
sumstats_median2$virus<-"small-DNA"
sumstats2$rel_trans<-sumstats2$trans/max(sumstats2$trans)
sumstats_median2$rel_trans<-sumstats_median2$trans/max(sumstats_median2$trans)


sumstats3<-read.csv(file="sumstats_mean_lDNA1_1.csv",header = TRUE)
sumstats_median3<-read.csv(file="sumstats_median_lDNA1_1.csv",header = TRUE)
sumstats3$virus<-"large-DNA"
sumstats_median3$virus<-"large-DNA"
sumstats3$rel_trans<-sumstats3$trans/max(sumstats3$trans)
sumstats_median3$rel_trans<-sumstats_median3$trans/max(sumstats_median3$trans)


names(sumstats2)

sumstats<-rbind(sumstats1,sumstats2[,-6],sumstats3[,-c(6,13)])
sumstats_median<-rbind(sumstats_median1,sumstats_median2[,-6],sumstats_median3[,-c(6,13)])

sumstats<-as.data.frame(sumstats)
sumstats_median<-as.data.frame(sumstats_median)

head(sumstats)


#############
#Plot Log(infection duration)
#############

ggplot(subset(sumstats,(epsilon2>=1)&(epsilon2<=5)), 
       aes(as.numeric(epsilon1),as.numeric(epsilon2),z=as.numeric(log10(tfinal))))+
  geom_raster(aes(fill = as.numeric(log10(tfinal))), interpolate = FALSE) +
  #geom_contour(colour = "white",binwidth = .2)+
  labs(fill=expression(paste(log[10] (T[fin]) )),x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="inferno")+
  theme_light()+
  theme(legend.title = element_text(size=8), panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#  guides(fill=guide_legend(keywidth = 0.5))+
  facet_wrap(~ as.factor(virus))

ggsave(plot=last_plot(),filename="T_final_viruses_div1_N1l.pdf", width = 13, height = 5.5, units = "cm")



#############
#Plot BH fitness
#############

#for all the number of classes (for comparison)     
ggplot(subset(sumstats,(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=log10(trans+1)))+
  geom_raster(aes(fill = log10(trans+1)), interpolate = FALSE) +
# geom_contour(colour = "white")+
  labs(fill="infection\n fitness",x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="plasma")+
  theme_light()+
  facet_wrap(~ as.factor(virus))+
  theme(legend.title = element_text(size=8), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  

ggsave(plot=last_plot(),filename="BH_Fitness_viruses_div1_N1l.pdf", width = 13, height = 5.5, units = "cm")



#############
#Plot relative BH fitness
#############


ggplot(subset(sumstats,(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=100*rel_trans))+
  geom_raster(aes(fill = 100*rel_trans), interpolate = FALSE) +
  # geom_contour(colour = "white")+
  labs(fill="% max transm",x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  theme(axis.text=element_text(size=10),axis.title=element_text(size=14),plot.title=element_text(size=12,face="bold"))+
  scale_fill_viridis_c(option="plasma")+
  theme_minimal()+
  facet_wrap(~ as.factor(virus))+
  theme(legend.position="right")


ggsave(plot=last_plot(),filename="rel_BH_Fitness_viruses_div1.pdf", width = 17, height = 6.5, units = "cm")


#############
#Plot WH fitness
#############

ggplot(subset(sumstats,(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=log10(Vtot)))+
  geom_raster(aes(fill = log10(Vtot+1)),interpolate = FALSE) +
  #  geom_contour(colour = "white")+
  labs(fill=expression(paste(log[10] (V[tot]) )),x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="viridis")+
  theme_light()+
  theme(legend.title = element_text(size=8))+
  facet_wrap(~ as.factor(virus))


# N decreases WH fitness (less reinfection)
ggsave(plot=last_plot(),filename="WH_Fitness_viruses_div1_N1.pdf", width = 17, height = 6.5, units = "cm")



#############
#Plot proportion cancer
#############


ggplot(subset(sumstats,(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=cancer))+
  geom_raster(aes(fill = cancer),interpolate = FALSE) +
#  geom_contour(colour = "white")+
  labs(fill=expression(paste(cancer)),x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="cividis")+
  theme_light()+
  theme(legend.title = element_text(size=8),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  facet_wrap(~ as.factor(virus))



ggsave(plot=last_plot(),filename="Cancer_viruses_div1_N1l.pdf", width = 13, height = 5.5, units = "cm")




#############
#Plot time to cancer
#############




ggplot(subset(sumstats,(cancer>0)&(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=tfinal/365))+
  geom_raster(aes(fill = tfinal/365),interpolate = FALSE) +
  #  geom_contour(colour = "white")+
  labs(fill="years\n to\n cancer",x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="inferno")+
  theme_light()+
  theme(legend.title = element_text(size=8),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  facet_wrap(~ as.factor(virus))



ggsave(plot=last_plot(),filename="Time_to_cancer_viruses_div1_N1l.pdf", width = 13, height = 5.5, units = "cm")



#############
#Plot time to recovery
#############

#for all the number of classes (for comparison)     
ggplot(subset(sumstats,(as.numeric(mut)==0)&(epsilon2>=1)&(epsilon2<=5)), 
       aes(epsilon1,epsilon2,z=tfinal/365))+
  geom_raster(aes(fill = tfinal/365),interpolate = FALSE) +
  #  geom_contour(colour = "white")+
  labs(fill="time\n to\n recovery",x=expression(paste(epsilon [1])),y=expression(paste(epsilon [2])))+
  scale_fill_viridis_c(option="magma")+
  theme_light()+
  theme(legend.title = element_text(size=8))+
  facet_wrap(~ as.factor(virus))


ggsave(plot=last_plot(),filename="Time_to_recovery_viruses_div1_N1.pdf", width = 17, height = 6.5, units = "cm")




#############
#Select runs with the highest fitness (mean)
#############

fraction<-0.025  #the fraction to keep
cancer_threshold<-0.0

length(sumstats1$trans)*fraction

# topTrans1<-quantile(sumstats1$trans,(1-fraction),na.rm=TRUE) #find the threshold for 10% best fitness
# bestfit1<-sumstats1[(sumstats1$trans>=topTrans1),] #select the runs with 10% highest fitness
# 
# topTrans2<-quantile(sumstats2$trans,(1-fraction),na.rm=TRUE) #find the threshold for 10% best fitness
# bestfit2<-sumstats2[(sumstats2$trans>=topTrans2),] #select the runs with 10% highest fitness
# 
# topTrans3<-quantile(sumstats3$trans,(1-fraction),na.rm=TRUE) #find the threshold for 10% best fitness
# bestfit3<-sumstats3[(sumstats3$trans>=topTrans3),] #select the runs with 10% highest fitness


topTrans1c<-quantile(sumstats1$trans[sumstats1$cancer>cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs with cancer
topTrans1n<-quantile(sumstats1$trans[sumstats1$cancer<=cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs without cancer

bestfit1c<-sumstats1[(sumstats1$trans>=topTrans1c)&(sumstats1$cancer>cancer_threshold),] #select the runs with 10% highest fitness
bestfit1n<-sumstats1[(sumstats1$trans>=topTrans1n)&(sumstats1$cancer<=cancer_threshold),] #select the runs with 10% highest fitness

bestfit1<-rbind(bestfit1c,bestfit1n)

topTrans2c<-quantile(sumstats2$trans[sumstats2$cancer>cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs with cancer
topTrans2n<-quantile(sumstats2$trans[sumstats2$cancer<=cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs without cancer

bestfit2c<-sumstats2[(sumstats2$trans>=topTrans2c)&(sumstats2$cancer>cancer_threshold),] #select the runs with 10% highest fitness
bestfit2n<-sumstats2[(sumstats2$trans>=topTrans2n)&(sumstats2$cancer<=cancer_threshold),] #select the runs with 10% highest fitness

bestfit2<-rbind(bestfit2c,bestfit2n)


quantile(sumstats3$trans[sumstats3$cancer>cancer_threshold])
quantile(sumstats3$trans[sumstats3$cancer<=cancer_threshold])

topTrans3c<-quantile(sumstats3$trans[sumstats3$cancer>cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs with cancer
topTrans3n<-quantile(sumstats3$trans[sumstats3$cancer<=cancer_threshold],(1-fraction),na.rm=TRUE) #fitness for fittest runs without cancer

bestfit3c<-sumstats3[(sumstats3$trans>=topTrans3c)&(sumstats3$cancer>cancer_threshold),] #select the runs with 10% highest fitness
bestfit3n<-sumstats3[(sumstats3$trans>=topTrans3n)&(sumstats3$cancer<=cancer_threshold),] #select the runs with 10% highest fitness

bestfit3<-rbind(bestfit3c,bestfit3n)


bestfit<-rbind(bestfit1,bestfit2[,-6],bestfit3[,-c(6,13)])


bestfit$cancer2<-bestfit$cancer>cancer_threshold

bestfit_mean<-as.data.frame(bestfit)


# Mean duration

ggplot(bestfit_mean, aes(x=tfinal/365,colour=as.factor(cancer2>0)))+
  geom_density()+ 
  labs(x="Infection duration (in years)",colour="cancer")+
  theme_light()+
  scale_color_grey()+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,tfin=mean(tfinal/365)),aes(xintercept=tfin, colour=as.factor(cancer2)),linetype="dashed", size=0.5)+
  facet_wrap(~ as.factor(virus))
  
  

#ggsave(plot=last_plot(),filename="best_Infection_duration_viruses_div1.pdf", width = 17, height = 6.5, units = "cm")
ggsave(plot=last_plot(),filename="best_Infection_duration_viruses_N1_025p0.pdf", width = 9, height = 6.5, units = "cm")




ggplot(bestfit_mean, aes(x=tfinal/365,colour=as.factor(virus)))+
  geom_density()+ 
  labs(x="Infection duration (in years)",colour="virus")+
  theme_light()+
  scale_color_viridis_d()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
geom_vline(data=ddply(bestfit_mean,"virus",summarize,med=mean(tfinal/365)),aes(xintercept=med, colour=as.factor(virus)),linetype="dashed", size=0.5)
  
ddply(bestfit_mean,"virus",summarize,med=mean(tfinal/365))

ggsave(plot=last_plot(),filename="best_Infection_duration_viruses_N1_025p0.pdf", width = 9, height = 6.5, units = "cm")


# Infection fitness


ggplot(bestfit_mean, aes(x=log10(trans+1),colour=as.factor(cancer2)))+
  geom_density()+ 
  labs(x="log10(infection fitness)",colour="cancer")+
  xlim(2.76,2.79)+
  theme_light()+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_grey()+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(log10(trans+1))),aes(xintercept=med, colour=as.factor(cancer2)),linetype="dashed", size=0.5)+
  facet_wrap(~ as.factor(virus),scale="free",ncol=1)

ggsave(plot=last_plot(),filename="best_fitness_viruses_N1_025p0t.pdf", width = 5.5, height = 11, units = "cm")



# total viral load


ggplot(bestfit_mean, aes(x=log10(Vtot),colour=as.factor(cancer2)))+
  geom_density()+ 
  labs(x="log10(Vtot)",colour="cancer")+
  theme_light()+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_grey()+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(log10(Vtot))),aes(xintercept=med, colour=as.factor(cancer2)),linetype="dashed", size=0.5)+
  facet_wrap(~ as.factor(virus),scale="free",ncol=1)

ggsave(plot=last_plot(),filename="best_Vtot_viruses_N1_025p0t.pdf", width = 5.5, height = 11, units = "cm")



# Relative nfection fitness

ggplot(bestfit_mean, aes(x=100*rel_trans,colour=as.factor(cancer2)))+
  geom_density()+ 
#  xlim(90,100)+
  labs(x="% maximum infection fitness",colour="cancer")+
  theme_light()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  scale_color_grey()+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(100*rel_trans)),aes(xintercept=med, colour=as.factor(cancer2)),linetype="dashed", size=0.5)+
  facet_wrap(~ as.factor(virus),scale="free")

ggsave(plot=last_plot(),filename="best_relative_fitness_viruses_N1_025p0.pdf", width = 17, height = 6.5, units = "cm")


# Proportion of cancers


ggplot(subset(bestfit_mean,cancer2>0), aes(x=cancer*100,color=virus))+
  geom_density()+ 
  labs(x="cancers (%)")+
  theme_light()+
  scale_color_viridis_d()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
geom_vline(data=ddply(subset(bestfit_mean,cancer2>0),"virus",summarize,med=mean(cancer*100)),aes(xintercept=med, colour=as.factor(virus)),linetype="dashed", size=0.5)


ggsave(plot=last_plot(),filename="best_proportion_cancers_viruses_N1_025p0.pdf", width = 5, height = 5, units = "cm")





# Time to cancer


ggplot(subset(bestfit_mean,cancer2>0), 
       aes(x=tfinal/365,color=virus))+
  geom_density()+ 
  labs(x="years to cancer")+
  scale_color_viridis_d()+
  theme_light()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position="none")+
  geom_vline(data=ddply(subset(bestfit_mean,cancer2>0),"virus",summarize,med=mean(tfinal/365)),aes(xintercept=med, colour=as.factor(virus)),linetype="dashed", size=0.5)


ggsave(plot=last_plot(),filename="best_Time_to_cancer_viruses_N1_025p0.pdf", width = 5, height = 5, units = "cm")




####"" epsilon values with or without cancer

ggplot(bestfit_mean, 
       aes(x=epsilon1,colour=as.factor(virus)))+
  geom_density()+ 
  labs(x=expression(paste(epsilon[1])),colour="cancer")+
  scale_color_viridis_d()+
  theme_light()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")+
  geom_vline(data=ddply(bestfit_mean,c("virus"),summarize,med=mean(epsilon1)),aes(xintercept=med, colour=as.factor(virus)),linetype="dashed", size=0.5)


ggsave(plot=last_plot(),filename="best_epsilon1_viruses_joined_N1_025p0.pdf", width = 5, height = 5, units = "cm")


ggplot(bestfit_mean, 
       aes(x=epsilon1,colour=as.factor(cancer2>0)))+
  geom_density()+ 
  labs(x=expression(paste(epsilon[1])),colour="cancer")+
  scale_color_grey()+
  theme_light()+
  facet_wrap(~ as.factor(virus),ncol=1)+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(epsilon1)),aes(xintercept=med, colour=as.factor(cancer2)),linetype="dashed", size=0.5)
  
ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(epsilon1))

ggsave(plot=last_plot(),filename="best_epsilon1_viruses_N1_025p0.pdf", width = 6.5, height = 11, units = "cm")


ggplot(bestfit_mean,  
       aes(x=epsilon2,colour=as.factor(cancer2>0)))+
  geom_density()+ 
  labs(x=expression(paste(epsilon[2])),colour="cancer")+
  scale_color_grey()+
  theme_light()+
  facet_wrap(~ as.factor(virus),ncol=1,scale="free")+
  theme(legend.position="none",panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  geom_vline(data=ddply(bestfit_mean,c("virus","cancer2"),summarize,med=mean(epsilon2)),aes(xintercept=med, colour=as.factor(cancer2)),linetype="dashed", size=0.5)
  


ggsave(plot=last_plot(),filename="best_epsilon2_viruses_N1_025p0.pdf", width = 5.5, height = 11, units = "cm")




