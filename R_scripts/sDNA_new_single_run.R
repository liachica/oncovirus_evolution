setwd("/home/samuel/Documents/Hote-parasite/HPV/PTRS/Evolution_oncoviruses/Gillespie/single_runs")


##########################
#load libraries
##########################
library(tictoc)


tic()

set.seed(578) 

##########################
#creating an ad hoc function to draw one number in a poisson distribution 
##########################

rpoisson<-function(x)
{
  rpois(1,x)
}


MaxDiv<-as.integer(30)    #maximum number of cellular divisions

powerrisk <- 1



epsilon1 <- 0.5

epsilon2 <- 3.5

Tmax <- 18250         #maximum number of days that is 50 years 

nu    <- 5*10^(-15) #nu = 10^(-12) - ^(-14)

omega <-0.001   
delta1 <-1.0/1.2
delta2 <-0.2/1.2
sigma <-0.02
mu    <-0.02
alpha <-0.001
gamma <-0.0000001 
K     <-50 
K2    <-10^(10) #G carrying capacity
phi   <-2 #5
theta <-0.5
lambda<-1

iruns<-10

results<-list(length=iruns)
  

for(i in c(1:iruns))  
{
 
  print(i)
  
  
    #initialise the variables for this run
  
  Gvec <- rep(0.0,MaxDiv)
  Gvec[1]<-10
  G<-sum(Gvec)
  Rvec <- rep(0.0,MaxDiv)
  R  <- sum(Rvec)
  P  <- 0
  Tc <- 1
  V  <- 0
  t  <- 0
  mut <- rep(0.0,MaxDiv)
  trans<-0
  Vtot<-0
  
  stop<-0
  
  #initialise results vector
  pop <- c(t,G,R,P,Tc,V,sum(mut),trans,Vtot)
  
  #time step is per day
  deltat <- 1.0
  
  
  #loop until max time reached or extinction or explosion
  while((t<=Tmax)&(stop==0))
  {
    #when is the next event happenning (drawn from an poisson distribution)
    
    infection   <- if(omega*V*(1-G/K2)*deltat>0){rpois(1, omega*V*(1-G/K2)*deltat)}else{0}; #G has a carrying capacity K2
    
    division_result1_tot <- if(delta1*R*deltat>0){rpois(1, delta1*R*deltat)}else{0};
    division_result2_tot <- if(delta2*R*deltat>0){rpois(1, delta2*R*deltat)}else{0};
    division_start_tot <- if(sigma*epsilon2*G*(1-G/K2)*deltat>0){rpois(1, sigma*epsilon2*G*(1-G/K2)*deltat)}else{0};#G has a carrying capacity K2
    nat_deathG_tot   <- if(mu*(1-epsilon1)*G*deltat>0){rpois(1, mu*(1-epsilon1)*G*deltat)}else{0};
    killingG_tot    <- if(alpha*G*Tc*deltat>0){rpois(1, alpha*G*Tc*deltat)}else{0};
    
    #splitting the events accross replication classes
    division_result1<-rep(0,MaxDiv)
    if(division_result1_tot>10^4)
    {
      # if there are too many cells we can't assign them individually to replication classes and just take proportions
      division_result1<-round(division_result1_tot*Rvec/R)
    }else
    {
      if(division_result1_tot>0)
      {
        temp_division_result1 <- sample(c(1:MaxDiv),size=division_result1_tot,prob=Rvec/R,replace=T);
        division_result1[as.numeric(names(table(temp_division_result1)))]<-table(temp_division_result1)
      }
      
    }
    
    
    
    division_result2<-rep(0,MaxDiv)
    if(division_result2_tot>0)
    {
      if(division_result2_tot<10^4)
      {
        temp_division_result2 <- sample(c(1:MaxDiv),size=division_result2_tot,prob=Rvec/R,replace=T);
        division_result2[as.numeric(names(table(temp_division_result2)))]<-table(temp_division_result2)
      }
      else{
        # if there are too many cells we can't assign them individually to replication classes and just take proportions
        division_result2<-round(division_result2_tot*Rvec/R)
      }
    }
    
    
    
    division_start<-rep(0,MaxDiv)
    if(division_start_tot>0)
    {
      if(division_start_tot<10^4)
      {
        temp_division_start <- sample(c(1:MaxDiv),size=division_start_tot,prob=Gvec/G,replace=T);
        division_start[as.numeric(names(table(temp_division_start)))]<-table(temp_division_start)
      }else{
        division_start<-round(division_start_tot*Gvec/G)
      }
    }      
    
    
    nat_deathG<-rep(0,MaxDiv)
    if(nat_deathG_tot>0)
    {
      if(nat_deathG_tot<10^4)
      {
        temp_nat_deathG <- sample(c(1:MaxDiv),size=nat_deathG_tot,prob=Gvec/G,replace=T);
        nat_deathG[as.numeric(names(table(temp_nat_deathG)))]<-table(temp_nat_deathG)
      }else{
        nat_deathG<-round(nat_deathG_tot*Gvec/G)
      }
    }
    
    
    
    killingG<-rep(0,MaxDiv)
    if(killingG_tot>0)
    {
      if(killingG_tot<10^4)
      {
        temp_killingG <- sample(c(1:MaxDiv),size=killingG_tot,prob=Gvec/G,replace=T);
        killingG[as.numeric(names(table(temp_killingG)))]<-table(temp_killingG)
      }else{
        killingG<-round(killingG_tot*Gvec/G)
      }
    }    
    
    
    nat_deathP  <- if(mu*P*deltat>0){rpois(1, mu*P*deltat)}else{0};
    killingP    <- if(alpha*P*Tc*deltat>0){rpois(1, alpha*P*Tc*deltat)}else{0};
    imm_growthG <- if(gamma*(G+R)*Tc*(1-Tc/K)*deltat>0){rpois(1, gamma*(G+R)*Tc*(1-Tc/K)*deltat)}else{0};
    imm_growthP <- if(gamma*P*Tc*(1-Tc/K)*deltat>0){rpois(1, gamma*P*Tc*(1-Tc/K)*deltat)}else{0};
    v_production<- if(phi*P*deltat>0){rpois(1, phi*P*deltat)}else{0};
    v_clearance <- if(theta*V*deltat>0){rpois(1, theta*V*deltat)}else{0};
    
    
    mut_process <- if(sum(Gvec)>0){unlist(lapply(nu*Gvec*(c(1:MaxDiv)^powerrisk)*deltat,rpoisson))}else{rep(0,MaxDiv)};
    
    
    Gvec[1]<-Gvec[1] + infection + 2*sum(division_result1) + sum(division_result2);      # new cells born after division are age 0
    Gvec<-Gvec+c(0,(division_result1[-MaxDiv]+division_result2[-MaxDiv]))                  # mother cells move to the next class 
    #    Gvec[MaxDiv]<-Gvec[MaxDiv]+ division_result1[MaxDiv] + division_result2[MaxDiv]               # if we remove this line, cells with MaxDiv do not produce offspring anymore
    Gvec  <- Gvec - nat_deathG - division_start - killingG;     # death terms
    
    Rvec  <- Rvec + division_start - division_result1 - division_result2;
    
    P  <- P + sum(division_result2) - nat_deathP - killingP;
    Tc <-Tc + imm_growthG + imm_growthP;
    V  <- V + v_production - v_clearance;
    mut<-mut+ mut_process;
    Vtot<-Vtot + sum(v_production);
    
    #check that numbers are positive for Rvec and Gvec (poisson can do weird stuff)
    Rvec<-Rvec*(Rvec>=0)
    Gvec<-Gvec*(Gvec>=0)
    
    #update total number of G cells and R cells    
    G<-sum(Gvec)
    R  <- sum(Rvec)
    
    #update time and store results
    t <- t+deltat
    pop <- rbind(pop,c(t,G,R,P,Tc,V,sum(mut),trans,Vtot))
    
    #check whether stopping conditions are met
    stop=(sum(mut)>0)+(G+R<=0);
    
    if(is.na(stop)){stop=100}
    
      

    }
    
    
colnames(pop) <- c("time","G","R","P","Tc","V","mut","trans","Vtot")

rownames(pop) <-c()    

results[[i]]<-pop

}


plot(rep(-10,length(c(0:round(Tmax/365))))~c(0:round(Tmax/365)),type="l",xlab="time (years)",ylab="log[10](population size +1)",xlim=c(0,10),ylim=c(0,10))
for(j in c(1:iruns))
{
  data<-as.data.frame(results[[j]])
  timeY<-as.numeric((data$time/365))
  if(sum(data$mut>0)){print(max(timeY))}
  lines(log10(1+data$G)~timeY,col=1,lwd=0.2)
  lines(log10(1+data$R)~timeY,col=2,lwd=0.2)
  lines(log10(1+data$V)~timeY,col=3,lwd=0.2)
  lines(log10(1+data$Tc)~timeY,col=4,lwd=0.2)
  lines(log10(1+data$P)~timeY,col=5,lwd=0.2)
 # abline(v=max(timeY),lty=2)
}



gdata<-c()
vdata<-c()
for(ig in c(1:iruns))
{
  data<-as.data.frame(results[[ig]])
  timeY<-as.numeric(as.character(data$time))
  Gs<-cbind(timeY,log10(1+data$G),"G",numGclasses,ig)
  Rs<-cbind(timeY,log10(1+data$R),"R",numGclasses,ig)
  Ps<-cbind(timeY,log10(1+data$P),"P",numGclasses,ig)
  Vs<-cbind(timeY,log10(1+data$V),"V",numGclasses,ig)
  Tc<-cbind(timeY,log10(1+data$Tc),"Tc",numGclasses,ig)
#  Cs<-cbind(max(timeY),0,"end",numGclasses,ig)
  gdata<-rbind(gdata,Gs,Rs,Ps,Vs,Tc)
  if(max(timeY)<Tmax)
  {
  vdata<-rbind(vdata,c(max(timeY),numGclasses,ig))
  }
}

colnames(gdata)<-c("time","popsize","type","N","run")
gdata<-as.data.frame(gdata)



write.csv(gdata, file = "single_run_sDNA_persist_(05-35).csv",row.names = F,quote=F)

