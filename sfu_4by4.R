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
#        assessment.window    --> Assessment window for toxicity and efficacy                         #
#        accrual.rate         --> Accrual rate of patients (patients per month)                       #
#        u11                  --> The specified utility value for (toxicity=1, efficacy=1)            #
#        u00                  --> The specified utility value for (toxicity=0, efficacy=0)            #
# ----------------------------------------------------------------------------------------------------#

##### 4*4 scenarios #####

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

sfu_4by4 <- function(rseed, ntrial, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
                     assessment.window, accrual.rate, u11, u00){
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
  
  get.joint.p<-function(yeff,ytox,rr){
    p.ab<-function(yeff,ytox,a,b,rr){
      yeff^a*(1-yeff)^(1-a)*ytox^b*(1-ytox)^(1-b)+(-1)^(a+b)*yeff*(1-yeff)*ytox*(1-ytox)*((exp(rr)-1)/(exp(rr)+1))
    }
    p=c(p.ab(yeff,ytox,a=1,b=0,rr),p.ab(yeff,ytox,a=0,b=0,rr),p.ab(yeff,ytox,a=1,b=1,rr),p.ab(yeff,ytox,a=0,b=1,rr))
    return(p)
  }
  
  ncohort=samplesize/cohortsize
  
  runin="TRUE";
  u.lp <- u11*target.eff + u00*(1-target.tox)
  u.bench <- 0.65*u.lp + 0.35*100
  u.true <- u11*p.true.eff + u00*(1-p.true.tox)
  
  ndose <- length(p.true.tox)
  TOX <- EST.TOX <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # toxicity outcome
  EFF <- EST.EFF <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # toxicity outcome
  N <- array(matrix(rep(0,length(p.true.tox) * ntrial),dim(p.true.tox)[1]),dim = c(dim(p.true.tox),ntrial)) # number of patients
  dselect = matrix(rep(0, 2 * ntrial), ncol = 2); # the selected dose level
  DUR <- rep(0,ntrial)
  niters <- 2000
  
  for(trial in 1:ntrial){
    
    set.seed(trial + rseed*10000)
    # calculate trial duration
    t.decision <- 0
    t.enroll <- t.wait <- NULL
    
    tox <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]); ## the number of DLT at each dose level
    eff <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);    ## the number of Eff at each dose level
    n <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);      ## the number of patients treated at each dose level
    event <- NULL
    earlystop = 0;         ## indiate whether the trial terminates early
    elimi.tox = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);  ## indicate whether doses are eliminated due to toxicity
    elimi.eff = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);  ## indicate whether doses are eliminated due to efficacy
    elimi = matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]);      ## indicate whether doses are eliminated due to toxicity or efficacy
    final.est.tox <- final.est.eff <- NULL
    est.tradeoff <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2]); 
    u.mean <- matrix(rep(0, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2])
    E.out.trial <- T.out.trial <- NULL
    
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
        joint.p=get.joint.p(p.true.eff[d[1],d[2]],p.true.tox[d[1],d[2]],rr=rr)
        ys <- rmultinom(cohortsize, size = 1, prob = joint.p)
        e.outcome=unlist(lapply(1:cohortsize,function(x) ifelse((ys[1,x]==1 | ys[3,x]==1),1,0)))
        t.outcome=unlist(lapply(1:cohortsize,function(x) ifelse((ys[3,x]==1 | ys[4,x]==1),1,0)))
        T.out.trial <- c(T.out.trial,t.outcome)
        E.out.trial <- c(E.out.trial,e.outcome)
        tox[d[1],d[2]] = tox[d[1],d[2]] + sum(t.outcome)
        eff[d[1],d[2]] = eff[d[1],d[2]] + sum(e.outcome)
        n[d[1],d[2]] = n[d[1],d[2]] + cohortsize
        if(tox[d[1],d[2]]>0 | eff[d[1],d[2]]>0 ){break}
      }
      ncohort2 <- ncohort-sum(n)/cohortsize
      
    }
    
    for(pp in startpp:ncohort2){
      if(pp>0){ 
        
        ## pp=0: final analysis of run-in stage 
        ##### data generation #####
        
        joint.p=get.joint.p(p.true.eff[d[1],d[2]],p.true.tox[d[1],d[2]],rr=rr)
        ys <- rmultinom(cohortsize, size = 1, prob = joint.p)
        e.outcome=unlist(lapply(1:cohortsize,function(x) ifelse((ys[1,x]==1 | ys[3,x]==1),1,0)))
        t.outcome=unlist(lapply(1:cohortsize,function(x) ifelse((ys[3,x]==1 | ys[4,x]==1),1,0)))
        
        T.out.trial <- c(T.out.trial,t.outcome)
        E.out.trial <- c(E.out.trial,e.outcome)
        tox[d[1],d[2]] = tox[d[1],d[2]] + sum(t.outcome)
        eff[d[1],d[2]] = eff[d[1],d[2]] + sum(e.outcome)
        n[d[1],d[2]] = n[d[1],d[2]] + cohortsize;
      }
      
      ##### modeling toxicity and efficacy #####
      modelstring <- "
        model {
          
          for(i in 1:16){         
            tox[i] ~ dbin(p[i],n[i])
            eff[i] ~ dbin(q[i],n[i])
          }
          
          p[1]  <- 1-theta.t[1]
          p[2]  <- 1-theta.t[1]*theta.t[2]
          p[3]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]
          p[4]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]
          p[5]  <- 1-theta.t[1]*theta.t[5]
          p[6]  <- 1-theta.t[1]*theta.t[2]*theta.t[5]
          p[7]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[5]
          p[8]  <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]
          p[9]  <- 1-theta.t[1]*theta.t[5]*theta.t[6]
          p[10] <- 1-theta.t[1]*theta.t[2]*theta.t[5]*theta.t[6]
          p[11] <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[5]*theta.t[6]
          p[12] <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]*theta.t[6]
          p[13] <- 1-theta.t[1]*theta.t[5]*theta.t[6]*theta.t[7]
          p[14] <- 1-theta.t[1]*theta.t[2]*theta.t[5]*theta.t[6]*theta.t[7]
          p[15] <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[5]*theta.t[6]*theta.t[7]
          p[16] <- 1-theta.t[1]*theta.t[2]*theta.t[3]*theta.t[4]*theta.t[5]*theta.t[6]*theta.t[7]
          
          logit(q[1])  <- theta.e[1]
          logit(q[2])  <- theta.e[1]+theta.e[2]
          logit(q[3])  <- theta.e[1]+theta.e[2]+theta.e[3]
          logit(q[4])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]
          logit(q[5])  <- theta.e[1]+theta.e[5]
          logit(q[6])  <- theta.e[1]+theta.e[2]+theta.e[5]
          logit(q[7])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[5]
          logit(q[8])  <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]
          logit(q[9])  <- theta.e[1]+theta.e[5]+theta.e[6]
          logit(q[10]) <- theta.e[1]+theta.e[2]+theta.e[5]+theta.e[6]
          logit(q[11]) <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[5]+theta.e[6]
          logit(q[12]) <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]+theta.e[6]
          logit(q[13]) <- theta.e[1]+theta.e[5]+theta.e[6]+theta.e[7]
          logit(q[14]) <- theta.e[1]+theta.e[2]+theta.e[5]+theta.e[6]+theta.e[7]
          logit(q[15]) <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[5]+theta.e[6]+theta.e[7]
          logit(q[16]) <- theta.e[1]+theta.e[2]+theta.e[3]+theta.e[4]+theta.e[5]+theta.e[6]+theta.e[7]
          
          for(i in 1:7){
            theta.t[i] ~ dbeta(a.t[i],b.t[i])T(0,0.999999)
            theta.e[i] ~ dnorm(a.e[i],1/b.e[i])T(-100,100)
          }
        }
        "
      
      jags.data <- list(tox=c(tox[1,],tox[2,],tox[3,],tox[4,]), eff=c(eff[1,],eff[2,],eff[3,],eff[4,]), n=c(n[1,],n[2,],n[3,],n[4,]),
                        a.t=tox.prior.a,b.t=tox.prior.b, a.e=eff.prior.mu,b.e=eff.prior.si)
      jags <- jags.model(textConnection(modelstring),data =jags.data,n.chains=1,n.adapt=5000,quiet=TRUE)
      # sampling toxicity
      t.sam <- coda.samples(jags,c('p'),niters,progress.bar="none") 
      t.samp <- as.matrix(t.sam);head(t.samp)
      post.mean.tox <- matrix(colMeans(t.samp),nrow=dim(p.true.tox)[1], byrow=T)
      post.prob.overtox <- matrix(colMeans(t.samp>(target.tox+0.05)),nrow=dim(p.true.tox)[1], byrow=T)
      post.elimi.tox <- matrix(colMeans(t.samp>target.tox),nrow=dim(p.true.tox)[1], byrow=T)
      # sampling efficacy
      e.sam <- coda.samples(jags,c('q'),niters,progress.bar="none") 
      e.samp <- as.matrix(e.sam)
      post.mean.eff <- matrix(colMeans(e.samp),nrow=dim(p.true.tox)[1], byrow=T)
      post.prob.undereff <- matrix(colMeans(e.samp<target.eff),nrow=dim(p.true.tox)[1], byrow=T)
      final.est.tox <- post.mean.tox
      final.est.eff <- post.mean.eff
      post00 <- (1-post.mean.tox)*(1-post.mean.eff)
      post01 <- (1-post.mean.tox)*post.mean.eff
      post10 <- post.mean.tox*(1-post.mean.eff)
      post11 <- post.mean.tox*post.mean.eff
      post.mean.u <- 100*post01+u11*post11+u00*post00
      post.prob.u <- matrix(colMeans((100*(1-t.samp)*e.samp+u11*t.samp*e.samp+u00*(1-t.samp)*(1-e.samp))>u.bench),nrow=dim(p.true.tox)[1], byrow=T)
      
      ##### continuous monitor #####
      # safety monitoring using toxicity model
      for(id in 1:dim(p.true.tox)[1]){
        for(jd in 1:dim(p.true.tox)[2]){
          elimi.tox[id,jd] <- max(elimi.tox[id,jd], as.numeric(post.elimi.tox[id,jd]>cutoff.eli.T))
          if (elimi.tox[id,jd]==1){
            for (i in min(id,dim(p.true.tox)[1]):dim(p.true.tox)[1]){
              for (j in min(jd,dim(p.true.tox)[2]):dim(p.true.tox)[2]) {
                elimi.tox[i,j] = 1}}
          }
        }
      }
      elimi <- elimi.tox
      if(sum(elimi.tox)==ndose){earlystop==1;break}
      if(post.elimi.tox[1,1]>(cutoff.eli.T-offset)){
        elimi.tox = matrix(rep(1, ndose),dim(p.true.tox)[1],dim(p.true.tox)[2])
        elimi <- elimi.tox
        earlystop=1
        break
      }
      # efficacy monitoring using efficacy model
      for(id in 1:dim(p.true.eff)[1]){
        for(jd in 1:dim(p.true.tox)[2]){
          if(post.prob.undereff[id,jd]>cutoff.eli.E){elimi.eff[id,jd] = 1}
        }
      }
      elimi[which(elimi.eff==1)] <- 1
      if(sum(elimi)==ndose){earlystop=1; break}
      u.mean <- post.prob.u*(1-elimi)
      
      ########################################
      #####         dose finding         #####
      ########################################
      
      ### 1. exclude overly toxic doses
      
      tox.monitor <- 1-elimi.tox
      
      ### 2. select MTD for j-th row and k-th column
      mtdc.r <- mtdc.c <- 0
      if(sum(which(tox.monitor[d[1],]==1))>0){mtdc.r <- max(which(tox.monitor[d[1],]==1))}
      if(sum(which(tox.monitor[,d[2]]==1))>0){mtdc.c <- max(which(tox.monitor[,d[2]]==1))}
      mtdc.r <- min(which(abs(post.mean.tox[d[1],]-target.tox)==min(abs(post.mean.tox[d[1],]-target.tox))), mtdc.r)
      mtdc.c <- min(which(abs(post.mean.tox[,d[2]]-target.tox)==min(abs(post.mean.tox[,d[2]]-target.tox))), mtdc.c)
      
      ### 3. decision making
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
    
    for(icoh in 1:(sum(n)/cohortsize)){
      T.out <- T.out.trial[(1+(icoh-1)*cohortsize):(icoh*cohortsize)]
      E.out <- E.out.trial[(1+(icoh-1)*cohortsize):(icoh*cohortsize)]
      for(ipt in 1:cohortsize){
        T.time=ifelse(T.out[ipt]==1, runif(1, 0, assessment.window[1]),assessment.window[1]+0.001)
        E.time=ifelse(E.out[ipt]==1, runif(1, 0, assessment.window[2]),assessment.window[2]+0.001)
        if(ipt == 1){ 
          
          if(icoh==1){ t.enroll <- c(t.enroll, 0 ) }
          if(icoh>=2){ t.enroll <- c(t.enroll, t.decision) }
          
        }else{ t.enroll <- c(t.enroll,t.enroll[length(t.enroll)]+runif(1, 0, 2*accrual.rate/cohortsize) ) }
        
        t.wait <- c(t.wait,max(T.time,E.time))
        
      }
      # print(t.enroll)
      if(icoh==ncohort) { 
        t.decision = t.decision + max(assessment.window)
      }else{
        t.decision <- max(t.wait+t.enroll)
      }
    }
    
    ##### data collection #####
    TOX[,,trial] <- tox 
    EFF[,,trial] <- eff
    N[,,trial] <- n
    DUR[trial] = t.decision
    EST.TOX[,,trial] <- final.est.tox 
    EST.EFF[,,trial] <- final.est.eff
    
    ##### final selection #####
    u.mean <- post.prob.u
    u.mean[n==0] = -100
    u.mean[elimi==1] = -100
    
    # select MTD for each row (fixing dosage of drug A)
    mtdc <- rep(0,dim(p.true.tox)[2])
    for(imtd in 1:dim(p.true.tox)[2]){
      mtdc[imtd] <- which(abs(final.est.tox[,imtd]-target.tox)==min(abs(final.est.tox[,imtd]-target.tox)))
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
      selpercent[i,j] = sum(dselect[,1] == i & dselect[,2] == j) / ntrial * 100;} };selpercent;n
  
  n.pts=round(apply(N,c(1,2),mean),2)
  n.tox=round(apply(TOX,c(1,2),mean),2)
  n.eff=round(apply(EFF,c(1,2),mean),2)
  # selection percentage of OBD
  uti.admi <- u.true*(p.true.tox<=target.tox)*(p.true.eff>target.eff)
  target.admi <- u.true*(p.true.tox<=target.tox)*(p.true.eff>=0.4)
  true.OBD <- which(uti.admi==max(uti.admi),arr.ind=T)
  sel.OBD <- pts.OBD <- sel.tar <- pts.tar <- 0
  if(max(uti.admi)>0){
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
                 duration = round(mean(DUR),2)  )
  return(relist)
  
}

p.true.tox <- matrix(c(0.12,0.16,0.18,0.20,  0.14,0.18,0.20,0.40,  0.18,0.20,0.42,0.55,  0.20,0.40,0.52,0.60),nrow = 4, byrow = TRUE)
p.true.eff <- matrix(c(0.20,0.25,0.30,0.50,  0.25,0.30,0.50,0.60,  0.30,0.50,0.60,0.66,  0.40,0.60,0.65,0.70),nrow = 4, byrow = TRUE)

sfu_4by4(rseed=1, ntrial=10, p.true.tox, p.true.eff, target.tox, target.eff, samplesize, cohortsize, 
         assessment.window, accrual.rate, u11, u00)
