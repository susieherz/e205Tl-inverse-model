# Thallium cycle inverse model
# Herz et al., 2025
# Modified from Kipp and Tissot, 2022: https://github.com/m-kipp/d238U-inverse-model.git

# Note: This code will take longer than a day to run

# loading necessary packages 

library(msir)
library(FME)
library(doParallel)
library(LaplacesDemon)

### PART I: READING IN DATA ###

# set working directory, read in dataset
#setwd("/Users/herz/Desktop")
current_dataset=read.csv("Herz_et_al_2026_best_sections.csv",header=T,stringsAsFactors=F)


# quick plot of dataset

par(mfrow=c(1,1))
par(mar=c(4,4,1,1))
par(oma=c(0,0,0,0))
plot(current_dataset$e205Tl~current_dataset$time,axes=F,xlab="",ylab="")
loess_mod=loess.sd(current_dataset$e205Tl~current_dataset$time,span=0.3,nsigma=2) 
polygon(c(loess_mod$x,rev(loess_mod$x)),c(loess_mod$upper,rev(loess_mod$lower)),col=rgb(0,0,0.8,0.1),border=F)
lines(loess_mod$y~loess_mod$x,col="blue")
arrows(current_dataset$time,current_dataset$e205Tl+(2*current_dataset$err),current_dataset$time,current_dataset$e205Tl-(2*current_dataset$err),length=0)
points(current_dataset$e205Tl~current_dataset$time,pch=21,bg="white")
axis(side=1,mgp=c(3,0.5,0),las=1,tck=-0.02)
axis(side=2,mgp=c(3,0.5,0),las=1,tck=-0.02)
mtext("time",side=1,line=2.25)
mtext(expression(epsilon^{205}*Tl),side=2,line=2.5)
box(lwd=1)


## model parameters

# things to definitely tune

time_step=100 # size of time-steps
prop_uncert=T # propagate uncertainty on terms in mass balance? yes=T, no=F [start no for testing, use yes at end]
m=10 # number of terms in Fourier series [start low, demonstrate convergence, increase if needed, test convergence again, repeat...]
niterMCMC=1e4 # number of steps per walker
updatecovMCMC=376 # number of steps after which covariance matrix is updated (only during burn-in phase) [see table from Kipp and Tissot, 2022 "ReadMe" file]
n_walkers=100 # number of parallel walkers deployed [aim for 3-5 in initial tests, 100+ if propagating uncertainty]

# things to maybe never tune (fingers crossed)

amp_limit=1e-4 # range for amplitude terms [decrease if low acceptance persists]
start=1 # start of frequency terms; 1 for odd harmonics or continuous, 2 for even harmonics [default=1]
count=1 # spacing of frequency terms; 2 for odd or even harmonics, 1 for continuous [default=2]
jump_step_1=amp_limit*0.05 # SD of proposal distribution for amplitude terms [default=amp_limit*0.05]
jump_step_2=0.2 # SD of proposal distribution for log_fanox_0 [default=0.2]
burninlengthMCMC=niterMCMC-1e3 # number of discarded "burn in" steps [1e3 less than niterMCMC typically works well]
thin=1 # factor by which final output is "thinned"; 1=no thinning, 10=keep every 10th sample [default=1]

# things to never tune

start_time=min(current_dataset$time) - 200000 # first time in current_dataset should be 0. sets a 200 kyr spin-up period
max_time=max(current_dataset$time)+(time_step-start_time)
a_vals=seq(0,0,length.out=m)
b_vals=seq(0,0,length.out=m)
as=paste("a",1:m,sep="")
bs=paste("b",1:m,sep="")
names(a_vals)=as
names(b_vals)=bs
pars1=c(log_fanox_0=log10(0.00211),a_vals,b_vals) # starting values; gives modern mass balance
upper=c(0,seq(amp_limit, amp_limit,length=2*m)) # upper limits
lower=c(-4,seq(-amp_limit,-amp_limit,length=2*m)) # lower limits
jump=c(jump_step_2,seq(jump_step_1,jump_step_1,length=2*m)) 

# steady-state matrix with e205Tl and f_anox values

mols=8.76e10 # mols of Tl in ocean (equal to 65 pmol/L; Rehkämper and Nielsen, 2004)
J.in=4e6 # mols/yr of Tl input flux (Owens, 2019)
K.basalt=(0.65*4e6)/(mols) # 65% of Tl inputs go towards low-T basalt alteration in modern mass balance (Owens, 2019)
K.anox=(J.in-K.basalt*mols)/(mols*0.00211 + 0.0125*mols*(1-0.00211)) # 5% of Tl inputs (Owens, 2019)
K.mnoxide=0.0125*K.anox # 30% of Tl inputs (Owens, 2019)
e205Tlin=-2 # Nielsen et al., 2017
D.anox=0 # Nielsen et al., 2017
D.mnoxide=16 # Rehkämper et al., 2002 
D.basalt=-1.2 # Nielsen et al., 2017
fanox_seq=10^seq(-5,0,length=1000)
log_fanox_seq=seq(-5,0,length=1000)

e205Tl_SS = e205Tlin - (((J.in/(K.anox*fanox_seq + K.mnoxide*(1-fanox_seq) + K.basalt)) * (K.anox*fanox_seq*D.anox + K.mnoxide*(1-fanox_seq)*D.mnoxide + K.basalt*D.basalt)) / J.in)
plot(e205Tl_SS)
ss_out=cbind(e205Tl_SS, log_fanox_seq)

### PART II: FORWARD MODEL ###
e205_model=function(pars){
  derivs=function(time,y,pars){
    with(as.list(c(pars,y)),{
      de205Tlsw=(J.in*(e205Tlin - e205Tlsw) - K.anox*Nsw*(10^log_fanox)*D.anox - K.mnoxide*Nsw*(1-(10^log_fanox))*D.mnoxide - K.basalt*Nsw*D.basalt)/Nsw
      dNsw= J.in - K.anox*Nsw*(10^log_fanox) - K.mnoxide*Nsw*(1-(10^log_fanox))	- K.basalt*Nsw
      dlog_fanox=if(log_fanox>-1e-5){-1e-8} else if(abs(log_fanox)>5){1e-8} else {sum((pars[seq(2,by=1,length.out=m)]*sin(seq(start,length.out=m,by=count)*pi*time/max_time))+(pars[seq(2+m,by=1,length.out=m)]*cos(seq(start,length.out=m,by=count)*pi*time/max_time)))}       
      return(list(c(de205Tlsw,dNsw,dlog_fanox)))
    })
  }
  
  log_fanox_0=with(as.list(pars),log_fanox_0)
  Nsw_0=with(as.list(pars),J.in/(K.anox*(10^log_fanox_0)+K.mnoxide*(1-(10^log_fanox_0))+K.basalt))
  e205Tlsw_0=with(as.list(pars),as.numeric((ss_out[ss_out[,2]>=log_fanox_0,])[1,1]))
  
  y=c(e205Tlsw=e205Tlsw_0,Nsw=Nsw_0,log_fanox=log_fanox_0)
  
  times=seq(start_time,max_time,by=time_step)
  
  out=ode(y=y,parms=pars,times=times,func=derivs,method="euler")
  
  as.data.frame(out)
  
}

# plotting baseline (=starting) model run (should be flatline at modern values)

out1=e205_model(pars=pars1)

par(mfrow=c(3,1))

plot(out1$log_fanox~out1$time,cex=0)
lines(out1$log_fanox~out1$time)
plot(out1$Nsw~out1$time,cex=0)
lines(out1$Nsw~out1$time)
plot(out1$e205Tlsw~out1$time,cex=0)

lines(out1$e205Tlsw~out1$time)
points(current_dataset$e205Tl~current_dataset$time,pch=21,bg="white")

### PART III: INVERSE MODEL ###

k_basalt=K.basalt # saves the value of the rate constants as defined for the steady-state scenario

# defining cost function (negative log-likelihood, NLL)

e205Tlcost <- function(pars) {
  tryCatch({
    out1 <- e205_model(pars = pars)
    merged <- merge(out1, current_dataset, by = "time")
    NLL <- log(2 * pi * merged$err) + ((merged$e205Tl - merged$e205Tlsw)^2) / (merged$err^2)
    return(sum(NLL, na.rm = TRUE))}, error = function(e) {    # returns a large cost if the model fails
      return(1e10)
    })
}

# setting up parallel computation 

cores <- max(1, future::availableCores()-1)
cl=makeCluster(cores)
registerDoParallel(cl)

# running MCMC routine

start.time=Sys.time()
finalMCMC=foreach(i=1:n_walkers,.combine=rbind) %dopar% {
  
  library(FME)
  
  if(prop_uncert==T){
    
    K.basalt=rnorm(1,k_basalt, k_basalt*0.1)
    
    K.anox=(J.in-K.basalt*mols)/(mols*0.00211 + 0.0125*mols*(1-0.00211))
    
    K.mnoxide=0.0125*K.anox
    
    e205Tlin=rnorm(1,-2,1)
    
    D.anox=rnorm(1,0,1) 
    
    D.mnoxide=rnorm(1,16,1)
    
    D.basalt=rnorm(1,-1.2,1) 
  }
  
  fanox_seq=10^seq(-5,0,length=1000)
  log_fanox_seq=seq(-5,0,length=1000)
  
  e205Tl_SS = e205Tlin - (((J.in/(K.anox*fanox_seq + K.mnoxide*(1-fanox_seq) + K.basalt)) * (K.anox*fanox_seq*D.anox + K.mnoxide*(1-fanox_seq)*D.mnoxide + K.basalt*D.basalt)) / J.in)
  
  ss_out=cbind(e205Tl_SS, log_fanox_seq)
  
  curMCMC=modMCMC(f=e205Tlcost,p=pars1,niter=niterMCMC,lower=lower,upper=upper,jump=jump,prior=NULL,var0=NULL,wvar0=0.1,updatecov=updatecovMCMC,burninlength=burninlengthMCMC)
  
  tempMatrix=cbind(curMCMC$par)
  tempMatrix
  
}

end.time=Sys.time()
run.time1=end.time-start.time
run.time1
stopCluster(cl)

# thinning sample set

thin_count=dim(finalMCMC)[1]/thin
i=1
new_out=NULL
for (i in 1:thin_count){
  cur_vals=finalMCMC[i*thin,]
  new_out=rbind(new_out,cur_vals)
}

write.csv(new_out,"MCMC_out.csv")
MCMC_out=read.csv("MCMC_out.csv")
MCMC_out=MCMC_out[,2:dim(MCMC_out)[2]]

# using compiled MCMC output to generate confidence interval for time series

start.time=Sys.time()

i=1
n_reps=1e3 # 1e3 is sufficient for final run; can decrease to 1e2 for testing 
time_seq=seq(start_time,max_time,by=time_step)
e205Tlsw_out=NULL
Nsw_out=NULL
fanox_out=NULL
for (i in 1:n_reps){
  
  if(prop_uncert==T){
    
    K.basalt=rnorm(1,k_basalt, k_basalt*0.1)
    
    K.anox=(J.in-K.basalt*mols)/(mols*0.00211 + 0.0125*mols*(1-0.00211))
    
    K.mnoxide=0.0125*K.anox
    
    e205Tlin=rnorm(1,-2,1)
    
    D.anox=rnorm(1,0,1) 
    
    D.mnoxide=rnorm(1,16,1)
    
    D.basalt=rnorm(1,-1.2,1) 
  }
  
  fanox_seq=10^seq(-5,0,length=1000)
  log_fanox_seq=seq(-5,0,length=1000)
  
  e205Tl_SS = e205Tlin - (((J.in/(K.anox*fanox_seq + K.mnoxide*(1-fanox_seq) + K.basalt)) * (K.anox*fanox_seq*D.anox + K.mnoxide*(1-fanox_seq)*D.mnoxide + K.basalt*D.basalt)) / J.in)
  ss_out=cbind(e205Tl_SS, log_fanox_seq)
  
  sr1=sensRange(fun=e205_model,parms=NULL,parInput= MCMC_out,num=1) 
  e205Tlsw_sens=sr1[,(length(pars1)+1):(length(pars1)+length(time_seq))]
  Nsw_sens=sr1[,(length(pars1)+1+length(time_seq)):(length(pars1)+2*length(time_seq))]
  fanox_sens=sr1[,(length(pars1)+1+2*length(time_seq)):(length(pars1)+3*length(time_seq))]
  
  e205Tlsw_out=rbind(e205Tlsw_out, e205Tlsw_sens)
  Nsw_out=rbind(Nsw_out, Nsw_sens)
  fanox_out=rbind(fanox_out, fanox_sens)
}

i=1
e205Tlsw_quant_out=NULL
Nsw_quant_out=NULL
fanox_quant_out=NULL

for (i in 1:length(time_seq)){
  cur_time=time_seq[i]
  cur_e205Tlsw_quant=quantile(na.omit(e205Tlsw_out[,i]),probs=seq(0.01,1,by=0.01))
  cur_Nsw_quant=quantile(na.omit(Nsw_out[,i]),probs=seq(0.01,1,by=0.01))
  cur_fanox_quant=quantile(na.omit(fanox_out[,i]),probs=seq(0.01,1,by=0.01))
  e205Tlsw_quant_out=rbind(e205Tlsw_quant_out, cur_e205Tlsw_quant)
  Nsw_quant_out=rbind(Nsw_quant_out, cur_Nsw_quant)
  fanox_quant_out=rbind(fanox_quant_out, cur_fanox_quant)
}

end.time=Sys.time()
run.time2=end.time-start.time
run.time2

### PART IV: PLOTTING e205Tlsw and f_anox WITH MCMC-DERIVED POSTERIOR PROBABILITY DISTRIBUTIONS ####

# calculating NLL, effective sample size, acceptance fraction, psrf

mcmcout=cbind(time_seq,e205Tlsw_quant_out[,50])
colnames(mcmcout)=c("time","e205Tlsw")
mergemcmc=merge(mcmcout,current_dataset,by="time") 
total_NLL=round((sum(0.5*log(2*pi*mergemcmc$err))+sum(na.omit(((mergemcmc$e205Tl-mergemcmc$e205Tlsw)^2)/(mergemcmc$err^2)))),digits=1)  # final NLL of median output
accept_frac=100*round(length(unique(finalMCMC[,1]))/length(finalMCMC[,1]),digits=2) # fraction of steps accepted during random walks
z=1
y=1
auto_out=NULL
big_auto_out=NULL
mcmc_list_out=NULL
for (z in 1:n_walkers){
  cur_walker_data=finalMCMC[(1+(((niterMCMC-burninlengthMCMC)/thin)*(z-1))):(((niterMCMC-burninlengthMCMC)/thin)*(z)),]
  
  y=1
  for (y in 1:(dim(cur_walker_data)[2]-2)){
    cur_autocorr=IAT(as.numeric(cur_walker_data[,y]))
    auto_out=rbind(auto_out,cur_autocorr)
  }
  
  cur_avg_autocorr=mean(auto_out)
  cur_eff=((niterMCMC-burninlengthMCMC)/thin)/cur_avg_autocorr
  cur_eff_size=mean(effectiveSize(cur_walker_data))
  results=c(cur_avg_autocorr, cur_eff,cur_eff_size)
  big_auto_out=rbind(big_auto_out, results)
  
  mcmc_results=cbind(cur_walker_data,z)
  mcmc_list_out=rbind(mcmc_list_out,mcmc_results)
  
  
}
avg_autocorr_MCMC=mean(big_auto_out[,1]) # average integrated autocorrelation time (IAT)
eff=round(sum(big_auto_out[,2]),digits=0) # effective sample size (total samples/IAT)

mcmc_out_list <- list()
for (i in 1:n_walkers) {
  mcmc_out_list[[i]] <- mcmc(mcmc_list_out[mcmc_list_out[, dim(mcmc_list_out)[2]] == i, ], thin = thin)
  mcmc_out_list[[i]] <- mcmc_out_list[[i]][, 1:length(pars1)]
}
mcmc_list <- as.mcmc.list(mcmc_out_list)
psrf=gelman.diag(mcmc_list,autoburnin=FALSE) # calculating potential scale reduction factor (psrf), i.e. Gelman-Rubin statistic, or Rhat
avg_psrf=round(mean(psrf$psrf[,1]),digits=2) # average psrf for all variables (gives quick check of convergence; aiming for <<1.2)


# plotting recovered time series

par(mfrow=c(4,1))
par(mar=c(0.5,1,0.5,1))
par(oma=c(2.5,3.5,5,0))
plot(current_dataset$e205Tl~current_dataset$time,cex=0,axes=F,xlab="",ylab="", ylim = c(-8, 5), xlim=c(0, max(current_dataset$time)))
polygon(c(time_seq,rev(time_seq)),c(e205Tlsw_quant_out[,16],rev(e205Tlsw_quant_out[,84])),border=F,col="grey80")
lines(e205Tlsw_quant_out[,50]~ time_seq,col="black")
arrows(current_dataset$time,current_dataset$e205Tl+current_dataset$err,current_dataset$time,current_dataset$e205Tl-current_dataset$err,length=0)
points(current_dataset$e205Tl~current_dataset$time,pch=21,bg="white",cex=1.0)
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1,)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(epsilon^{205}*Tl),side=2,line=2.5,cex=0.9)
box(lwd=1)
mtext(side=3,line=3.5,"NLL:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=3.5,paste(total_NLL),at=max(current_dataset$time),adj=1,cex=0.8)
mtext(side=3,line=2.5,"fraction accepted:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=2.5,paste(accept_frac,"%"),at=max(current_dataset$time),col=if((accept_frac>=15)&(accept_frac<=40)){"green4"} else if((accept_frac<15)&(accept_frac>=10)){"orange"} else if((accept_frac>40)&(accept_frac<=50)){"orange"} else{"red"},adj=1,cex=0.8)
mtext(side=3,line=1.5,"avg psrf:",at=min(current_dataset$time),adj=0,cex=0.8)
mtext(side=3,line=1.5,paste(avg_psrf),at=max(current_dataset$time),col=if(avg_psrf<=1.1){"green4"} else if((avg_psrf>1.1)&(avg_psrf<=1.2)){"orange"} else{"red"},adj=1,cex=0.8)

plot(Nsw_quant_out[,50]~time_seq,xlim=c(0, max(current_dataset$time)), cex=0,axes=F,xlab="",ylab="")
polygon(c(time_seq,rev(time_seq)),c((Nsw_quant_out[,16]),rev(Nsw_quant_out[,84])),border=F,col="grey80")
lines(Nsw_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(italic("N"[sw])*" [mol]"),side=2,line=2.5,cex=0.9)
box(lwd=1)

Fmnoxide_quant_out=(K.mnoxide*(1-10^fanox_quant_out))/(10^fanox_quant_out*K.anox + K.mnoxide*(1-10^fanox_quant_out) + K.basalt)
plot(Fmnoxide_quant_out[,50]~time_seq, xlim=c(0, max(current_dataset$time)), cex=0,axes=F,xlab="",ylab="")
polygon(c(time_seq,rev(time_seq)),c((Fmnoxide_quant_out[,16]),rev(Fmnoxide_quant_out[,84])),border=F,col="grey80")
lines(Fmnoxide_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1,labels=NA)
mtext(expression(italic("Fmnox")),side=2,line=2.5,cex=0.9)
box(lwd=1)

plot(10^fanox_quant_out[,4]~time_seq, cex=0, ylim=c(0,1), xlim=c(0, max(current_dataset$time)), axes=F,xlab="", ylab="")
polygon(c(time_seq,rev(time_seq)),c((10^fanox_quant_out[,16]),rev(10^fanox_quant_out[,84])),border=F,col="grey80")
lines(10^fanox_quant_out[,50]~time_seq,col="black")
axis(side=2,tck=-0.02,mgp=c(3,0.5,0),las=1)
axis(side=1,tck=-0.02,mgp=c(3,0.5,0),las=1, labels = c("253.0", "252.8", "252.6", "252.4", "252.2", "252.0", "251.8", "251.6", "251.4", "251.2", "251.0", "250.8", "250.6", "250.4"), at = c(10000, 210000, 410000, 610000, 810000, 1010000, 1210000, 1410000, 1610000, 1810000, 2010000, 2210000, 2410000, 2610000))
mtext(expression(italic("f"[anox])),side=2,line=2.5,cex=0.9)
mtext("time [Ma]",side=1,line=2,cex=0.9)
box(lwd=1)

Fred_quant_out=(K.anox*(10^fanox_quant_out))/(10^fanox_quant_out*K.anox + K.mnoxide*(1-10^fanox_quant_out) + K.basalt)
Fbasalt_quant_out=(K.basalt)/(10^fanox_quant_out*K.anox + K.mnoxide*(1-10^fanox_quant_out) + K.basalt)

# saving model outputs

#setwd('/Users/herz/Desktop/tl_model_results') # optional: save results to separate folder

write.csv(fanox_quant_out, "fanox.csv", row.names = FALSE) # fraction of seafloor *not* oxygenated enough to support Mn oxide burial
write.csv(Nsw_quant_out, "Nsw.csv", row.names = FALSE) # mols of Tl in seawater
write.csv(e205Tlsw_quant_out, "e205Tlsw.csv", row.names = FALSE) # e205Tl of seawater
write.csv(Fmnoxide_quant_out, "Fmnoxide.csv", row.names = FALSE) # fraction of Tl inputs going towards Mn oxide burial
write.csv(Fbasalt_quant_out, "Fbasalt.csv", row.names = FALSE) # fraction of Tl inputs going towards low-T basalt alteration
write.csv(Fred_quant_out, "Fred.csv", row.names = FALSE) # fraction of Tl inputs going towards reducing sediments
write.csv(time_seq, "time_seq.csv", row.names = FALSE) # time sequence
