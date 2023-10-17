#######################################################################################################
#   SFU: Surface-Free Utility-Based Design for Dose Optimization in Cancer Drug-Combination Trials    #
#                        Jingyi Zhang, Nolan A. Wages∗ and Ruitao Lin*                                #
#                       Nolan.Wages@vcuhealth.org and rlin@mdanderson.org                             #
#######################################################################################################
# ----------------------------------------------------------------------------------------------------#
#        Main input variables                                                                         #
#        ntrial               --> number of replications                                              #
#        p.true.tox           --> True toxicity rate of drug combination                              #
#        p.true.eff           --> True efficacy rate of drug combination                              #
#        target.tox           --> Upper limit of toxicity rate                                        #   
#        target.eff           --> Lower limit of efficacy rate                                        #
#        samplesize           --> Maximum sample size                                                 #
#        cohortsize           --> Number of patients in one cohort                                    #
#        u11                  --> The specified utility value for (toxicity=1, efficacy=1)            #
#        u00                  --> The specified utility value for (toxicity=0, efficacy=0)            #
#        assessment.window    --> Assessment window for toxicity and efficacy                         #
#        accrual.rate         --> Accrual rate of patients (patients per month)                       #
#        maxpen               --> Maximum proportion of patients who have pending DLT or efficacy     #
#                                 outcomes at the current dose                                        #
# ----------------------------------------------------------------------------------------------------#

##### 5*3 scenarios #####

library(rjags)
library(coda)

ntrial <- 5000
samplesize <- 51
cohortsize <- 3
target.tox=0.35
target.eff=0.2
assessment.window <- c(1,2)
accrual.rate=1.5
u11=60
u00=40
maxpen=0.5

sfu_tite_3by5 <- function(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
                              assessment.window, accrual.rate, u11, u00, maxpen){
  
  library(rjags)
  library(coda)
  
  rr <- 0 # correlation between toxicity and efficacy
  cutoff.eli.T=0.95 # probability cutoff for toxicity
  offset=0.3 # probability cutoff for toxicity 
  cutoff.eli.E=0.9 # probability cutoff for toxicity
  
  # toxicity prior specification
  tox.prior.a=c(3,rep(3.6,6))
  tox.prior.b=c(1,rep(0.4,6))
  # efficacy prior specification
  eff.prior.mu = c(-2,rep(1,6))
  eff.prior.si=rep(2,7)
  
  model.fit <- function(tox.event, eff.event, dv, t.decision, t.enter, n,
                        tox.prior.a, tox.prior.b, eff.prior.mu, eff.prior.si, niters, target.tox, u.bench){
    modelstring <- "
        model {
          
          for(i in 1:npts){         
            tox.event[i] ~ dbern(p.patient[i])
            eff.event[i] ~ dbern(q.patient[i])
            p.patient[i]  <- (p[dv[i]])*w.t[i]
            q.patient[i]  <- (q[dv[i]])*w.e[i]
          }
          
          p[1]  <- 1-theta.t[1]
          p[2]  <- 1-theta.t[1]*theta.t[2]
          p[3]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]
          p[4]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]
          p[5]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]
          p[6]  <- 1-theta.t[1]*theta.t[6]
          p[7]  <- 1-theta.t[1]*theta.t[2]*theta.t[6]
          p[8]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[6]
          p[9]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[6]
          p[10]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]*theta.t[6]
          p[11]  <- 1-theta.t[1]*theta.t[6]*theta.t[7]
          p[12]  <- 1-theta.t[1]*theta.t[2]*theta.t[6]*theta.t[7]
          p[13]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[6]*theta.t[7]
          p[14]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[6]*theta.t[7]
          p[15]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]*theta.t[6]*theta.t[7]
          
          logit(q[1])  <- theta.e[1]
          logit(q[2])  <- theta.e[1]+theta.e[2]
          logit(q[3])  <- theta.e[1]+theta.e[2]+theta.e[3]
          logit(q[4])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]
          logit(q[5])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]
          logit(q[6])  <- theta.e[1]+theta.e[6]
          logit(q[7])  <- theta.e[1]+theta.e[2]+theta.e[6]
          logit(q[8])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[6]
          logit(q[9])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[6]
          logit(q[10])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]+theta.e[6]
          logit(q[11])  <- theta.e[1]+theta.e[6]+theta.e[7]
          logit(q[12])  <- theta.e[1]+theta.e[2]+theta.e[6]+theta.e[7]
          logit(q[13])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[6]+theta.e[7]
          logit(q[14])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[6]+theta.e[7]
          logit(q[15])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]+theta.e[6]+theta.e[7]
          
          for(i in 1:7){
            theta.t[i] ~ dbeta(a.t[i],b.t[i])T(0,0.999999)
            theta.e[i] ~ dnorm(a.e[i],1/b.e[i])T(-100,100)
          }
        }
        "
    jags.data <- list(tox.event=tox.event, eff.event=eff.event, dv=dv,npts=sum(n),
                      w.t=pmin((t.decision-t.enter)/assessment.window[1],1),w.e=pmin((t.decision-t.enter)/assessment.window[2],1),
                      a.t=tox.prior.a,b.t=tox.prior.b, a.e=eff.prior.mu,b.e=eff.prior.si)
    jags <- jags.model(textConnection(modelstring),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
    
    post.sam <- coda.samples(jags,c('q'),n.iter=niters,progress.bar="none") 
    e.samp <- as.matrix(post.sam)
    post.mean.eff <- matrix(colMeans(as.matrix(post.sam)),nrow=dim(p.true.tox)[1], byrow=T)
    post.prob.undereff <- matrix(colMeans(e.samp<target.eff),nrow=dim(p.true.tox)[1], byrow=T)
    
    post.sam <- coda.samples(jags,c('p'),n.iter=niters,progress.bar="none") 
    t.samp <- as.matrix(post.sam)
    post.mean.tox <- matrix(colMeans(as.matrix(post.sam)),nrow=dim(p.true.tox)[1], byrow=T)
    post.prob.overtox <- matrix(colMeans(t.samp>(target.tox+0.05)),nrow=dim(p.true.tox)[1], byrow=T)
    post.elimi.tox <- matrix(colMeans(t.samp>target.tox),nrow=dim(p.true.tox)[1], byrow=T)
    post.prob.u <- matrix(colMeans((100*(1-t.samp)*e.samp+u11*t.samp*e.samp+u00*(1-t.samp)*(1-e.samp))>u.bench),
                          nrow=dim(p.true.tox)[1], byrow=T)
    
    return(list(post.mean.tox=post.mean.tox, post.prob.overtox=post.prob.overtox, post.elimi.tox=post.elimi.tox,
                post.mean.eff=post.mean.eff, post.prob.undereff=post.prob.undereff,
                post.prob.u=post.prob.u
    )
    )
  }
  
  dosefinding <- function(elevel, d, elimi, u.mean){
    pr_H0 = rep(0,length(elevel) / 2)
    for (i in seq(1,length(elevel) / 2,by = 1)){
      if (d[1] + elevel[1,i] <= dim(p.true.tox)[1] &&
          d[2] + elevel[2,i] <= dim(p.true.tox)[2] &&
          d[1] + elevel[1,i] > 0 && d[2] + elevel[2,i] > 0){
        if (elimi[d[1] + elevel[1,i],d[2] + elevel[2,i]] == 0){
          pr_H0[i]<-u.mean[d[1] + elevel[1,i],d[2] + elevel[2,i]]+10^(-6)*(elevel[1,i]+elevel[2,i])}}}
    
    if (max(pr_H0) == 0) {d = d} else{
      k = which(pr_H0 == max(pr_H0))[as.integer(runif(1) * length(which(pr_H0 == max(pr_H0))) + 1)];
      d = d + c(elevel[1,k],elevel[2,k]);
    }
    return(d)
  }
  
  ncohort=samplesize/cohortsize
  runin="TRUE"
  dose.find.rule=2
  u.lp <- u11*target.eff + u00*(1-target.tox)
  u.bench <- 0.65*u.lp + 0.35*100
  u.true <- u11*p.true.eff + u00*(1-p.true.tox)
  
  ndose <- length(p.true.tox)
  TOX <- ELIMI.TOX <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # toxicity outcome
  EFF <- ELIMI.EFF <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # toxicity outcome
  N <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # number of patients
  dselect = matrix(rep(0, 2 * ntrial), ncol = 2); # the selected dose level
  locationM <- matrix(1:length(p.true.tox),nrow=dim(p.true.tox)[1],byrow=T)
  durationV <- NULL
  niters=2000
  
  for(trial in 1:ntrial){
    
    set.seed(trial + rseed*10000)
    
    tox <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]); ## the number of DLT at each dose level
    eff <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);    ## the number of Eff at each dose level
    n <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);      ## the number of patients treated at each dose level
    event <- NULL
    earlystop = 0;         ## indiate whether the trial terminates early
    elimi.tox = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);  ## indicate whether doses are eliminated due to toxicity
    elimi.eff = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);  ## indicate whether doses are eliminated due to efficacy
    elimi = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);      ## indicate whether doses are eliminated due to toxicity or efficacy
    u.mean <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2])
    
    t.enter <- NULL # time enter the study
    t.tox.event <- t.eff.event <- NULL # time to event
    tox.event <- eff.event <- NULL # patient-level outcome
    t.decision = 0; # decision making time
    dv1 <- dv2 <- dv <- NULL #dose of drug A / B for each subject
    npend <- 0
    model.fit.final <- NULL
    
    startpp <- ifelse(runin, 0,1)
    ncohort2 <- ncohort
    if(runin){
      
      traceRunin <- c(1,1)
      d <- c(1,1)
      for(i in 1:ncohort){
        if(d[1]<dim(p.true.tox)[1] & d[2]<dim(p.true.tox)[2]){
          rand2 <- runif(1)
          if(rand2 < 0.5) {d <- c(d[1],d[2]+1)}else{d <- c(d[1]+1,d[2])}
          traceRunin <- rbind(traceRunin,d)
        }
        if(d[1]==dim(p.true.tox)[1] & d[2]<dim(p.true.tox)[2]){
          d <- c(d[1],d[2]+1)
          traceRunin <- rbind(traceRunin,d)
        }
        if(d[1]<dim(p.true.tox)[1] & d[2]==dim(p.true.tox)[2]){
          d <- c(d[1]+1,d[2])
          traceRunin <- rbind(traceRunin,d)
        }
        if(d[1]==dim(p.true.tox)[1] & d[2]==dim(p.true.tox)[2]){break}
      }
      
      
      # until a toxicity OR efficacy event
      
      for(runt in 1:dim(traceRunin)[1]){
        d[1] <- traceRunin[runt,1]
        d[2] <- traceRunin[runt,2]
        d1 <- d
        for(j in 1:cohortsize){
          if(j==1) { t.enter = c(t.enter, t.decision) }else {
            t.enter = c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2/accrual.rate))} # uniform
        }
        t.decision = t.enter[length(t.enter)]
        
        # generate data
        t.outcome <- as.numeric(runif(cohortsize) < p.true.tox[d[1],d[2]]) # true toxicity outcome for each subject
        e.outcome <- as.numeric(runif(cohortsize) < p.true.eff[d[1],d[2]]) # true efficacy outcome for each subject
        t.tox=ifelse(t.outcome==1, runif(1, 0, assessment.window[1]),assessment.window[1]+0.001) # generate time to toxicity for each subject 
        t.eff=ifelse(e.outcome==1, runif(1, 0, assessment.window[2]),assessment.window[2]+0.001) # generate time to efficacy for each subject 
        
        # patient-level outcome
        t.tox.event = c(t.tox.event, t.tox)
        t.eff.event = c(t.eff.event, t.eff)
        tox.event = c(tox.event, t.outcome)
        eff.event = c(eff.event, e.outcome)
        dv1 = c(dv1, rep(d[1], cohortsize))
        dv2 = c(dv2, rep(d[2], cohortsize))
        dv = c(dv, rep(locationM[d[1],d[2]], cohortsize))
        
        tox[d[1],d[2]] = tox[d[1],d[2]] + sum(t.outcome)
        eff[d[1],d[2]] = eff[d[1],d[2]] + sum(e.outcome)
        n[d[1],d[2]] = n[d[1],d[2]] + cohortsize
        
        pending=1
        npend = npend-1
        while(pending==1)
        {
          npend = npend+1
          pending = 0
          
          t.decision = t.decision + runif(1, 0, 2/accrual.rate)#update decision time
          cset = (dv1==d[1] & dv2==d[2]) #current dose indicator
          
          # indicator for whether or not toxicity and efficacy is observed. O means pending.
          delta.tox = as.numeric(((t.enter+t.tox.event)<t.decision))[cset]
          delta.eff = as.numeric(((t.enter+t.eff.event)<t.decision))[cset]
          
          #number of patients who have pending outcomes
          n.pend = length(unique(c(which(delta.tox==0),which(delta.eff==0))))
          
          #check whether the trial should be suspended, where the number of patients on current dose is sum(cset)
          if(n.pend>sum(cset)*maxpen) {pending=1}else{
            
            ##### fitting model before transform to the main phase
            if(tox[d[1],d[2]]>0 | eff[d[1],d[2]]>0 | runt == dim(traceRunin)[1]){
              
              ##update utility values for doses that have been used to treat patients
              model.fit.runin <- model.fit(tox.event, eff.event, dv, t.decision, t.enter, n,
                                           tox.prior.a, tox.prior.b, eff.prior.mu, eff.prior.si, niters, target.tox, u.bench)
              model.fit.final <- model.fit.runin
              
              ##### continuous monitor #####
              
              # safety monitoring using toxicity model
              for(id in 1:dim(p.true.tox)[1]){
                for(jd in 1:dim(p.true.tox)[2]){
                  elimi.tox[id,jd] <- max(elimi.tox[id,jd], as.numeric(model.fit.final$post.elimi.tox[id,jd]>cutoff.eli.T))
                  if (elimi.tox[id,jd]==1){
                    for (i in min(id,dim(p.true.tox)[1]):dim(p.true.tox)[1]){
                      for (j in min(jd,dim(p.true.tox)[2]):dim(p.true.tox)[2]) {
                        elimi.tox[i,j] = 1}}
                  }
                }
              }
              elimi <- elimi.tox
              if(sum(elimi.tox)==ndose){earlystop==1;break}
              if(model.fit.final$post.elimi.tox[1,1]>(cutoff.eli.T-offset)){
                elimi.tox = matrix(rep(1, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2])
                elimi <- elimi.tox
                earlystop=1
                break
              }
              
              # efficacy monitoring using efficacy model
              for(id in 1:dim(p.true.eff)[1]){
                for(jd in 1:dim(p.true.tox)[2]){
                  if(model.fit.final$post.prob.undereff[id,jd]>cutoff.eli.E){elimi.eff[id,jd] = 1}
                }
              }
              elimi[which(elimi.eff==1)] <- 1
              if(sum(elimi)==ndose){earlystop=1; break}
              u.mean <- model.fit.final$post.prob.u*(1-elimi)
              
              ##### dose finding #####
              
              ### 1. exclude overly toxic doses
              
              tox.monitor <- 1-elimi.tox
              
              ### 2. select MTD for j-th row and k-th column
              mtdc.r <- mtdc.c <- 0
              if(sum(which(tox.monitor[d[1],]==1))>0){mtdc.r <- max(which(tox.monitor[d[1],]==1))}
              if(sum(which(tox.monitor[,d[2]]==1))>0){mtdc.c <- max(which(tox.monitor[,d[2]]==1))}
              mtdc.r <- min(which(abs(model.fit.final$post.mean.tox[d[1],]-target.tox)==min(abs(model.fit.final$post.mean.tox[d[1],]-target.tox))), mtdc.r)
              mtdc.c <- min(which(abs(model.fit.final$post.mean.tox[,d[2]]-target.tox)==min(abs(model.fit.final$post.mean.tox[,d[2]]-target.tox))), mtdc.c)
              
              ### 3. decision making
              dnext <- d
              if (mtdc.r > d[2] & mtdc.c > d[1]){
                elevel = matrix(c(1,0,0,1, -1,0,0,-1, 0,0),2)
                dnext <- dosefinding(elevel, d, elimi, u.mean)
              }
              if (mtdc.r == d[2] & mtdc.c > d[1]){
                elevel = matrix(c(-1,0,0,-1, 0,0, 1,0),2)
                dnext <- dosefinding(elevel, d, elimi, u.mean)
              } 
              if (mtdc.r > d[2] & mtdc.c == d[1]){
                elevel = matrix(c(-1,0,0,-1, 0,0, 0,1),2)
                dnext <- dosefinding(elevel, d, elimi, u.mean)
              }
              if ((mtdc.r == d[2] & mtdc.c == d[1]) | (mtdc.r == d[2] & mtdc.c < d[1]) | (mtdc.r < d[2] & mtdc.c == d[1])){
                elevel = matrix(c(-1,0,0,-1, 0,0),2)
                dnext <- dosefinding(elevel, d, elimi, u.mean)
              }
              if (mtdc.r < d[2] & mtdc.c < d[1]){
                elevel = matrix(c(-1,0,0,-1),2)
                dnext <- dosefinding(elevel, d, elimi, u.mean)
              }
              d <- dnext
              
            }
            
          }
          if(earlystop == 1){break}
        }
        if(tox[d1[1],d1[2]]>0 | eff[d1[1],d1[2]]>0 ){break}
      }
      ncohort2 <- ncohort-sum(n)/cohortsize
      
    }
    if(sum(elimi)==ndose){
      
      ##### data collection #####
      t.decision = max(t.enter)+max(assessment.window)
      TOX[,,trial] <- tox 
      EFF[,,trial] <- eff
      N[,,trial] <- n
      ELIMI.TOX[,,trial] <- elimi.tox
      ELIMI.EFF[,,trial] <- elimi.tox
      durationV <- c(durationV, t.decision)
      
      ##### final selection #####
      dselect[trial,] <- c(99,99)
      next
    }
    for(pp in 1:ncohort2){
      
      ## pp=0: final analysis of run-in stage 
      ##### data generation #####
      
      # d is current dose
      for(j in 1:cohortsize){
        if(j==1) { t.enter = c(t.enter, t.decision); }
        else {t.enter = c(t.enter, t.enter[length(t.enter)] + runif(1, 0, 2/accrual.rate))} # uniform
      }
      
      t.decision = t.enter[length(t.enter)]
      
      # generate data
      t.outcome <- as.numeric(runif(cohortsize) < p.true.tox[d[1],d[2]]) # true toxicity outcome for each subject
      e.outcome <- as.numeric(runif(cohortsize) < p.true.eff[d[1],d[2]]) # true efficacy outcome for each subject
      t.tox=ifelse(t.outcome==1, runif(1, 0, assessment.window[1]),assessment.window[1]+0.001) # generate time to toxicity for each subject 
      t.eff=ifelse(e.outcome==1, runif(1, 0, assessment.window[2]),assessment.window[2]+0.001) # generate time to efficacy for each subject 
      
      # patient-level outcome
      t.tox.event = c(t.tox.event, t.tox)
      t.eff.event = c(t.eff.event, t.eff)
      tox.event = c(tox.event, t.outcome)
      eff.event = c(eff.event, e.outcome)
      dv1 = c(dv1, rep(d[1], cohortsize))
      dv2 = c(dv2, rep(d[2], cohortsize))
      dv = c(dv, rep(locationM[d[1],d[2]], cohortsize))
      
      tox[d[1],d[2]] = tox[d[1],d[2]] + sum(t.outcome)
      eff[d[1],d[2]] = eff[d[1],d[2]] + sum(e.outcome)
      n[d[1],d[2]] = n[d[1],d[2]] + cohortsize;
      
      pending=1;
      npend = npend-1;
      while(pending==1)
      {
        npend = npend+1;
        pending = 0;
        
        #update decision time
        if(pp==ncohort2) { t.decision = t.decision + max(assessment.window)}else {
          t.decision = t.decision + runif(1, 0, 2/accrual.rate)
        }
        
        cset = (dv1==d[1] & dv2==d[2]); #current dose indicator
        
        # indicator for whether or not toxicity and efficacy is observed. O means pending.
        delta.tox = as.numeric(((t.enter+t.tox.event)<t.decision))[cset]
        delta.eff = as.numeric(((t.enter+t.eff.event)<t.decision))[cset]
        
        #number of patients who have pending outcomes
        n.pend = length(unique(c(which(delta.tox==0),which(delta.eff==0))));n.pend
        
        #check whether the trial should be suspended, where the number of patients on current dose is sum(cset)
        if(n.pend>sum(cset)*maxpen) {pending=1;}else{
          
          ##update utility values for doses that have been used to treat patients
          model.fit.main <- model.fit(tox.event, eff.event, dv, t.decision, t.enter, n,
                                      tox.prior.a, tox.prior.b, eff.prior.mu, eff.prior.si, niters, target.tox, u.bench)
          model.fit.final <- model.fit.main
          
          ##### continuous monitor #####
          
          # safety monitoring using toxicity model
          for(id in 1:dim(p.true.tox)[1]){
            for(jd in 1:dim(p.true.tox)[2]){
              
              elimi.tox[id,jd] <- max(elimi.tox[id,jd], as.numeric(model.fit.final$post.elimi.tox[id,jd]>cutoff.eli.T))
              if (elimi.tox[id,jd]==1){
                for (i in min(id,dim(p.true.tox)[1]):dim(p.true.tox)[1]){
                  for (j in min(jd,dim(p.true.tox)[2]):dim(p.true.tox)[2]) {
                    elimi.tox[i,j] = 1}}
              }
            }
          }
          elimi <- elimi.tox
          if(sum(elimi.tox)==ndose){earlystop==1;break}
          if(model.fit.final$post.elimi.tox[1,1]>(cutoff.eli.T-offset)){
            elimi.tox = matrix(rep(1, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2])
            elimi <- elimi.tox
            earlystop=1
            break
          }
          
          # efficacy monitoring using efficacy model
          for(id in 1:dim(p.true.eff)[1]){
            for(jd in 1:dim(p.true.tox)[2]){
              if(model.fit.final$post.prob.undereff[id,jd]>cutoff.eli.E){elimi.eff[id,jd] = 1}
            }
          }
          elimi[which(elimi.eff==1)] <- 1
          if(sum(elimi)==ndose){earlystop=1; break}
          u.mean <- model.fit.final$post.prob.u*(1-elimi)
          
          ##### dose finding #####
          
          ### 1. exclude overly toxic doses
          
          tox.monitor <- 1-elimi.tox
          
          ### 2. select MTD for j-th row and k-th column
          mtdc.r <- mtdc.c <- 0
          if(sum(which(tox.monitor[d[1],]==1))>0){mtdc.r <- max(which(tox.monitor[d[1],]==1))}
          if(sum(which(tox.monitor[,d[2]]==1))>0){mtdc.c <- max(which(tox.monitor[,d[2]]==1))}
          mtdc.r <- min(which(abs(model.fit.final$post.mean.tox[d[1],]-target.tox)==min(abs(model.fit.final$post.mean.tox[d[1],]-target.tox))), mtdc.r)
          mtdc.c <- min(which(abs(model.fit.final$post.mean.tox[,d[2]]-target.tox)==min(abs(model.fit.final$post.mean.tox[,d[2]]-target.tox))), mtdc.c)
          
          ### 3. decision making
          dnext <- d
          if (mtdc.r > d[2] & mtdc.c > d[1]){
            elevel = matrix(c(1,0,0,1, -1,0,0,-1, 0,0),2)
            dnext <- dosefinding(elevel, d, elimi, u.mean)
          }
          if (mtdc.r == d[2] & mtdc.c > d[1]){
            elevel = matrix(c(-1,0,0,-1, 0,0, 1,0),2)
            dnext <- dosefinding(elevel, d, elimi, u.mean)
          } 
          if (mtdc.r > d[2] & mtdc.c == d[1]){
            elevel = matrix(c(-1,0,0,-1, 0,0, 0,1),2)
            dnext <- dosefinding(elevel, d, elimi, u.mean)
          }
          if ((mtdc.r == d[2] & mtdc.c == d[1]) | (mtdc.r == d[2] & mtdc.c < d[1]) | (mtdc.r < d[2] & mtdc.c == d[1])){
            elevel = matrix(c(-1,0,0,-1, 0,0),2)
            dnext <- dosefinding(elevel, d, elimi, u.mean)
          }
          if (mtdc.r < d[2] & mtdc.c < d[1]){
            elevel = matrix(c(-1,0,0,-1),2)
            dnext <- dosefinding(elevel, d, elimi, u.mean)
          }
          d <- dnext
          
        }
        if(earlystop == 1){break}
      }
      
      if(earlystop==1){t.decision = max(t.enter)+max(assessment.window); break;}
    }
    
    ##### data collection #####
    if(earlystop==1){t.decision = max(t.enter)+max(assessment.window)}
    TOX[,,trial] <- tox 
    EFF[,,trial] <- eff
    N[,,trial] <- n
    ELIMI.TOX[,,trial] <- elimi.tox
    ELIMI.EFF[,,trial] <- elimi.eff
    durationV <- c(durationV, t.decision)
    
    ##### final selection #####
    u.mean <- model.fit.final$post.prob.u
    u.mean[n==0] = -100
    u.mean[elimi==1] = -100
    
    # select MTD for each row (fixing dosage of drug A)
    mtdc <- rep(0,dim(p.true.tox)[2])
    for(imtd in 1:dim(p.true.tox)[2]){
      mtdc[imtd] <- which(abs(model.fit.final$post.mean.tox[,imtd]-target.tox)==min(abs(model.fit.final$post.mean.tox[,imtd]-target.tox)))
      if(mtdc[imtd]<dim(p.true.tox)[1]){
        u.mean[(mtdc[imtd]+1):dim(p.true.tox)[1],imtd] <- -100
      }
    }
    if(sum(u.mean==-100)==length(p.true.tox)){earlystop=1}
    if(earlystop==1){dselect[trial,] <- c(99,99)}else{
      dopt = which(u.mean==max(u.mean),arr.ind = TRUE)
      dselect[trial,1] = dopt[1];
      dselect[trial,2] = dopt[2];
    }
  }
  
  selpercent = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);
  for (i in 1:dim(p.true.tox)[1]){
    for (j in 1:dim(p.true.tox)[2]){
      selpercent[i,j] = sum(dselect[,1] == i & dselect[,2] == j) / ntrial * 100;} }
  
  n.pts=round(apply(N,c(1,2),mean),2)
  n.tox=round(apply(TOX,c(1,2),mean),2)
  n.eff=round(apply(EFF,c(1,2),mean),2)
  
  # selection percentage of OBD
  uti.admi <- u.true*(p.true.tox<=target.tox)*(p.true.eff>=target.eff)
  target.admi <- u.true*(p.true.tox<=target.tox)*(p.true.eff>=0.4)
  true.OBD <- which(uti.admi==max(uti.admi),arr.ind=T)
  sel.OBD <- pts.OBD <- sel.tar <- pts.tar <- 0
  if(length(uti.admi)>0){
    sel.OBD <- sum(selpercent[true.OBD])
    pts.OBD <- sum(n.pts[true.OBD])
    sel.tar <- sum(selpercent[which(target.admi>0)])
    pts.tar <- sum(n.pts[which(target.admi>0)])
  }
  
  # number of patients treated at overdose
  pts.overdose <- sum((p.true.tox>target.tox)*n.pts)
  sel.overdose <- sum((p.true.tox>target.tox)*selpercent)
  
  relist <- list(u.true=u.true,
                 selpercent=selpercent,
                 n.pts=n.pts, n.tox=n.tox, n.eff=n.eff,
                 totaltox=round(sum(TOX) / ntrial,2),
                 totaleff=round(sum(EFF) / ntrial,2),
                 totaln= round(sum(N) / ntrial, 2),
                 duration = round(mean(durationV),2) )
  
  return(relist)
  
}

# p.true.tox <- matrix(c(0.15,0.21,0.30,0.42,0.44,  0.24,0.30,0.42,0.44,0.51,  0.30,0.42,0.44,0.51,0.55),nrow = 3, byrow = TRUE)
# p.true.eff <- matrix(c(0.18,0.35,0.50,0.52,0.54,  0.35,0.50,0.52,0.53,0.57,  0.50,0.52,0.54,0.56,0.60),nrow = 3, byrow = TRUE)
# 
# sfu_tite_3by5(rseed=1, ntrial=10, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
#          assessment.window, accrual.rate, u11, u00, maxpen)