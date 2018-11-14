rm(list=ls())
# load libraries for data manipulations and figures 
library(plyr)
library(reshape2)
library(zoo)
library(mefa)
library(survival)
library(ggplot2)
library(grid)


###########################################################################
# Function for simulating survival data over one dose interval as a first
# step before such data are aggregated over multiple dose intervals
# as described in Section 2.4.1 of Huang et al. (2018)
###########################################################################
#' f_infect_singleDose
#' Input parameters: 
#' @param data  : simulated AMP-like infusion and non-infusion visit schedules according to steps described in Zhang et al. (2018) 
#' "Pharmacokinetics Simulations for Studying Correlates of Prevention Efficacy of Passive HIV-1 Antibody Prophylaxis in the Antibody 
#' Mediated Prevention (AMP) Study." arXiv:1801.08626, 2018. 
#' @param beta.t  : per-day change effect of the time-varying covariate on log hazard ratio
#' @param ts : time (in days) since the most recent infusion when drug concentration reaches the zero-protection threshold, s
#' @param r_incidence : annual HIV incidence rate in the placebo group
#' 
#' Output: 
#' A dataset containing variables
#' @ ID = participant ID
#' @ dose = participant randomly allocated dose group (10=10 mg/kg VRC01, 30=30 mg/kg VRC01, 0=placebo)
#' @ visit_n = protocol visit number
#' @ f_vis = flag variable indicating visit at which HIV diagnostic testing is performed (1=yes, 0=no) 
#' @ day_vis = time since enrollment (days)
#' @ f_infu = flag variable indicating visit at which infusion is administered (1=yes, 0=no)
#' @ f_event = flag variable indicating whether the survival event is observed (=1) or censored (=0)
#' @ t_event = time to event for cases or time to last follow up for non-cases 

f_infect_singleDose <- function(data
                                ,beta.t
                                ,ts
                                ,r_incidence = 0.04){
  lambda <- r_incidence/365/exp(beta.t*ts)
  
  
  ddply(data,.(ID),function(df){
              # browser()
    dat <- df
    dat <- dat[order(dat$day_vis),]
    #infection status
    dat_infu <- subset(dat,f_infu==1)
    if(dat$f_infu[nrow(dat)]==0){
      dat_lastVis <- dat[nrow(dat),]
      dat_infu <- rbind(dat_infu,dat_lastVis)
    }
    dat_infu$neg_log_u <- -log(runif(nrow(dat_infu)))
    dat_infu$cond <- lambda/beta.t*(exp(beta.t*ts)-1)
    dat_infu$t_event <- ifelse(dat_infu$neg_log_u<dat_infu$cond
                               ,log(1+beta.t*dat_infu$neg_log_u/lambda)/beta.t
                               ,dat_infu$neg_log_u/(lambda*exp(beta.t*ts))+(1-exp(beta.t*ts))/(beta.t*exp(beta.t*ts))+ts)
    
    dat_infu$t_interval <- dplyr::lead(dat_infu$day_vis)-dat_infu$day_vis
    dat_infu$f_event <- ifelse(dat_infu$t_interval>=dat_infu$t_event,1,0)
    # uninfected 
    if(all(dat_infu$f_event!=1,na.rm=T)){
      dat$t_event <- max(dat$day_vis)
      dat$f_event <- 0
    }
    # infected
    if(any(dat_infu$f_event==1,na.rm=T)){
      #             browser()
      dat_event <- subset(dat_infu,f_event==1)
      dat_event <- dplyr::sample_n(dat_event,1)
      dat$t_event <- dat_event$t_event+dat_event$day_vis
      dat$f_event <- ifelse(dat$t_event<= dat$day_vis,1,0)
    }
    return(dat)
  })
  
}

##############################################################################
# Function for simulating survival data over multiple dose intervals directly
# as described in Section 2.4.2 of Huang et al. (2018)
##############################################################################
#' gen_survival_time
#' Input parameters: 
#' @param beta.t  : per-day change effect of the time-varying covariate on log hazard ratio
#' #' Output: 
#' A vector of survival times in a length of n

# sample size 
n <- 3000
# incidence rate 
r <- 0.04/365

# perfect adherence infusion times (since enrolment) for the 10 infusions with an constant interval of 56 days
t.vec <- seq(0, 72*7, by=56)

# random uniform number  
u <- runif(n)

gen_survival_time <- function(beta.t){
  lambda <- r/exp(beta.t*56)
  t.ge.t1.lt.t2 <- log(1-beta.t*log(u)/lambda)/beta.t
  t.ge.t2.lt.t3 <- log(exp(beta.t*t.vec[2])*(-beta.t*log(u)/lambda-exp(beta.t*t.vec[2])+2))/beta.t
  t.ge.t3.lt.t4 <- log(exp(beta.t*t.vec[3])*(-beta.t*log(u)/lambda-(exp(beta.t*t.vec[2])+exp(beta.t*(t.vec[3]-t.vec[2])))+3))/beta.t
  # the following is only true when the dosing intervals are a constant 56 days. 
  t.ge.t4.lt.t5 <- log(exp(beta.t*t.vec[4])*(-beta.t*log(u)/lambda-(3* exp(beta.t*56))+4))/beta.t
  t.ge.t5.lt.t6 <- log(exp(beta.t*t.vec[5])*(-beta.t*log(u)/lambda-(4* exp(beta.t*56))+5))/beta.t
  t.ge.t6.lt.t7 <- log(exp(beta.t*t.vec[6])*(-beta.t*log(u)/lambda-(5* exp(beta.t*56))+6))/beta.t
  t.ge.t7.lt.t8 <- log(exp(beta.t*t.vec[7])*(-beta.t*log(u)/lambda-(6* exp(beta.t*56))+7))/beta.t
  t.ge.t8.lt.t9 <- log(exp(beta.t*t.vec[8])*(-beta.t*log(u)/lambda-(7* exp(beta.t*56))+8))/beta.t
  t.ge.t9.lt.t10 <-log(exp(beta.t*t.vec[9])*(-beta.t*log(u)/lambda-(8* exp(beta.t*56))+9))/beta.t
  t.ge.t10 <-      log(exp(beta.t*t.vec[10])*(-beta.t*log(u)/lambda-(9* exp(beta.t*56))+10))/beta.t
  a1 <- 0
  b1 <- lambda*(exp(beta.t*56)-1)/beta.t
  a2 <- b1
  b2 <- lambda*(2*exp(beta.t*56)-2)/beta.t
  a3 <- b2
  b3 <- lambda*(3*exp(beta.t*56)-3)/beta.t
  a4 <- b3
  b4 <- lambda*(4*exp(beta.t*56)-4)/beta.t
  a5 <- b4
  b5 <- lambda*(5*exp(beta.t*56)-5)/beta.t
  a6 <- b5
  b6 <- lambda*(6*exp(beta.t*56)-6)/beta.t
  a7 <- b6
  b7 <- lambda*(7*exp(beta.t*56)-7)/beta.t
  a8 <- b7
  b8 <- lambda*(8*exp(beta.t*56)-8)/beta.t
  a9 <- b8
  b9 <- lambda*(9*exp(beta.t*56)-9)/beta.t
  a10 <- b9
  
  t.infect <- rep(NA,n)
  
  for (i in 1:n)
    t.infect[i] <- ifelse(-log(u[i])<b1, t.ge.t1.lt.t2[i], 
                          ifelse(-log(u[i])>=a2 & -log(u[i]) < b2, t.ge.t2.lt.t3[i], 
                                 ifelse(-log(u[i])>=a3 & -log(u[i]) < b3, t.ge.t3.lt.t4[i], 
                                        ifelse(-log(u[i])>=a4 & -log(u[i]) < b4, t.ge.t4.lt.t5[i], 
                                               ifelse(-log(u[i])>=a5 & -log(u[i]) < b5, t.ge.t5.lt.t6[i], 
                                                      ifelse(-log(u[i])>=a6 & -log(u[i]) < b6, t.ge.t6.lt.t7[i],
                                                             ifelse(-log(u[i])>=a7 & -log(u[i]) < b7, t.ge.t7.lt.t8[i],
                                                                    ifelse(-log(u[i])>=a8 & -log(u[i]) < b8, t.ge.t8.lt.t9[i],
                                                                           ifelse(-log(u[i])>=a9 & -log(u[i]) < b9, t.ge.t9.lt.t10[i], t.ge.t10[i])))))))))
  
  return(t.infect)
}

###################
# Main simulation function
##################
#' sim
#' Input parameters: 
#' @param n_infu  : number of infusions
#' @param n_male  : number of male participants
#' @param n_female : number of female participants
#' @param prob_misInfu   : probability of an independently missed single infusion visit
#' @param prob_misVis : probability of an independently missed non-infusion visit
#' @param r_infuDiscont : cumulative probability of permanent infusion discontinuation
#' @param r_dropout  : annual dropout rate
#' @param beta.t.30 :  per-day change effect on log-hazard for 30 mg/Kg dose group
#' @param beta.t.10 : per-day change effect on log-hazard for 10 mg/Kg dose group
#' @param r_incidence : annual incidence rate in the placebo group
#' @param B  : number of datasets to simulate
#' 
#' #' Output: 
#' A dataset containing variables
#' @ ID = participant ID
#' @ dose = participant randomly allocated dose group (10=10 mg/kg VRC01 ,30=30 mg/kg VRC01,0=placebo)
#' @ inf = HIV infection (1=infected, 0=uninfected)
#' @ survtime = time to event for HIV-infected participants or time to last follow up for HIV-uninfected participants



sim <- function(n_infu=10
                ,n_male
                ,n_female
                ,prob_misInfu
                ,prob_misVis
                ,r_infuDiscont
                ,r_dropout
                ,beta.t.30
                ,beta.t.10
                ,r_incidence = 0.04
                ,B=1){
  #### redefine some parameters #####
  n_ppt <- n_male + n_female
  prob_misInfu <- rep(prob_misInfu,n_infu)
  prob_outWin <- rep(0,n_infu)
  prob_disContInfu_Win <- c(0/2,1-0,0/2)
 
  #directory to save the current scenario data
  folder <- paste0("m_",n_ppt
                   ,"_p1_",unique(prob_misInfu)
                   ,"_p2_",prob_misVis
                   ,"_r1_",r_infuDiscont
                   ,"_r2_",r_dropout
                   ,"_b30_",beta.t.30
                   ,"_b10_",beta.t.10
                   ,"_r3_",r_incidence
                   )                      
  ##### start the simulation #####
  res_list <- lapply(1:B
           ,function(b){
    set.seed(b)
    # browser()
    #-------- infusion visits -------------------
    dat <- data.frame(ID=1:n_ppt)
    dat$infu1 <- 1
    for(i in 2:n_infu){
      var_name <- paste0("infu",i)
      dat[,var_name] <- rbinom(n=n_ppt,size=1,prob=1-prob_misInfu[i])
    }
    #long format
    dat <- reshape2::melt(dat,id.vars="ID")
    dat <- dat[order(dat$ID),]
    dat <- rename(dat,c("value"="f_infu","variable"="infusion"))
    dat$infusion <- as.numeric(gsub("infu","",dat$infusion))
    #out of target window but within allowable window
    dat <- ddply(dat, .(infusion), function(df){
      if(unique(df$infusion)==1){
        df$out_win <- NA
      }else{
        infu <- unique(df$infusion)
        dat_infu <- subset(df,f_infu==1)
        dat_infu$out_win <- rbinom(n=nrow(dat_infu),size=1,prob=prob_outWin[infu])
        dat_misInfu <- subset(df,f_infu==0)
        if(nrow(dat_misInfu)==0){
          df <- dat_infu
        }else{
          dat_misInfu$out_win <- NA
          df <- rbind(dat_infu,dat_misInfu)
        }
      }
      return(df)
    })
    #random draw from windows
    dat_inWin <- subset(dat,out_win==0)
    dat_inWin$window <- round(runif(n=nrow(dat_inWin),min=-7,max=7),0)
    dat_outWin <- subset(dat,out_win==1)
    dat_outWin$window <- round(runif(n=nrow(dat_outWin),min=8,max=48),0)
    dat_naWin <- subset(dat,is.na(out_win))
    dat_naWin$window <- 0
    dat <- rbind(dat_inWin,dat_outWin,dat_naWin)
    rm(dat_inWin,dat_outWin,dat_naWin)
    dat <- dat[order(dat$ID,dat$infusion),]
    #infusion day
    dat <- ddply(dat,.(ID),function(df){
      #   browser()
      df$day_infu[df$infusion==1] <- 0
      for(i in 2:n_infu){
        df$day_infu[df$infusion==i] <- df$day_infu[df$infusion==(i-1)]+56+df$window[df$infusion==i]
      }
      return(df)
    })
    dat$day_infu <- ifelse(dat$f_infu==0,NA,dat$day_infu)
    dat$out_win <- NULL
    dat$window <- NULL
    #------------ post-infusion visits ------------------
#     browser()
    dat$window <- round(runif(n=nrow(dat),min=-7,max=7),0)
    dat$day_4wkPostInfu <- dat$day_infu+28+dat$window
    dat$day_8wkPostInfu10 <- dat$day_infu+56+dat$window
    dat$window_5dayPostInfu2 <- round(runif(n=nrow(dat),min=-2,max=2),0)
    dat$day_5dayPostInfu2 <- dat$day_infu+5+dat$window_5dayPostInfu2
    #missing visit
    dat <- ddply(dat,.(infusion),function(df){
      nppt_infu <- length(df$f_infu[df$f_infu==1])
      df$f_4wkPostInfuVis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)

      if(unique(df$infusion)==2){
        df$f_5dayPostInfu2Vis <- ifelse(df$f_infu==1,rbinom(n=nppt_infu,size=1,prob=1-prob_misVis),0)
      }else{
        df$f_5dayPostInfu2Vis <- 0
      }

      if(unique(df$infusion)==10){
        df$f_8wkPostInfu10Vis <- rbinom(n=nrow(df),size=1,prob=1-prob_misVis)
      }else{
        df$f_8wkPostInfu10Vis <- 0
      }
      return(df)

    })
    dat$day_4wkPostInfu <- ifelse(dat$f_4wkPostInfuVis==0,NA,dat$day_4wkPostInfu)
    dat$day_8wkPostInfu10 <- ifelse(dat$f_8wkPostInfu10Vis==0,NA,dat$day_8wkPostInfu10)
    dat$day_5dayPostInfu2 <- ifelse(dat$f_5dayPostInfu2Vis==0,NA,dat$day_5dayPostInfu2)
    dat$window <- NULL
    dat$window_5dayPostInfu2 <- NULL
    #------------ change data to long format ----------------
    #infusion visit
    dat_infu <- subset(dat,select=c(ID,infusion,f_infu,day_infu))
    dat_infu <- rename(dat_infu,c("infusion"="visit_n","f_infu"="f_vis","day_infu"="day_vis"))
    dat_infu$visit_n <- ifelse(dat_infu$visit_n%in%c(1,2),dat_infu$visit_n*2,dat_infu$visit_n*2+1)
    dat_infu$f_infu <- 1
    #4 week post infusion
    dat_4wkPostInfu <- subset(dat,select=c(ID,infusion,f_4wkPostInfuVis,day_4wkPostInfu))
    dat_4wkPostInfu <- rename(dat_4wkPostInfu,c("infusion"="visit_n"
                                                ,"f_4wkPostInfuVis"="f_vis"
                                                ,"day_4wkPostInfu"="day_vis"))
    dat_4wkPostInfu$visit_n <- ifelse(dat_4wkPostInfu$visit_n==1,3,dat_4wkPostInfu$visit_n*2+2)
    dat_4wkPostInfu$f_infu <- 0
    #8 week post 10th infusion
    dat_8wkPostInfu10 <- subset(dat
                              ,infusion==10
                              ,select=c(ID,infusion,f_8wkPostInfu10Vis,day_8wkPostInfu10))
    dat_8wkPostInfu10 <- rename(dat_8wkPostInfu10,c("infusion"="visit_n"
                                                      ,"f_8wkPostInfu10Vis"="f_vis"
                                                      ,"day_8wkPostInfu10"="day_vis"))
    dat_8wkPostInfu10$visit_n <- 23
    dat_8wkPostInfu10$f_infu <- 0
    # 5 days post 2nd infusion
    dat_5dayPostInfu2 <- subset(dat
                                ,infusion==2
                                ,select=c(ID,infusion,f_5dayPostInfu2Vis,day_5dayPostInfu2))
    dat_5dayPostInfu2 <- rename(dat_5dayPostInfu2,c("infusion"="visit_n"
                                                    ,"f_5dayPostInfu2Vis"="f_vis"
                                                    ,"day_5dayPostInfu2"="day_vis"))
    dat_5dayPostInfu2$visit_n <- 5
    dat_5dayPostInfu2$f_infu <- 0
    #combine
    rm(dat)
    dat_long <- rbind(dat_infu,dat_5dayPostInfu2,dat_4wkPostInfu,dat_8wkPostInfu10)
    rm(dat_infu,dat_5dayPostInfu2,dat_4wkPostInfu,dat_8wkPostInfu10)
    dat_long <- subset(dat_long,f_vis==1)
    dat_long <- dat_long[order(dat_long$ID,dat_long$visit_n),]
    #----------- permanent infusion discontinuation -------------------------------------------------
    #time to discontinuation of infusion
    if(r_infuDiscont==0){
      r_infuDiscont <- 0.00000001
    }
    time_discontInfu <- data.frame(ID=1:n_ppt,day_discontInfu=round(rexp(n=n_ppt,rate=r_infuDiscont/365),0))
    #don't allow discontiune at time 0
    time_discontInfu$day_discontInfu[time_discontInfu$day_discontInfu==0] <- 1
    dat_long <- merge(dat_long,time_discontInfu,by="ID")
    dat_long$f_discontInfu <- ifelse(dat_long$day_vis<dat_long$day_discontInfu,0,1)
    dat_long$f_vis <- ifelse(dat_long$f_discontInfu==1,0,dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$f_discontInfu==1, NA,dat_long$day_vis)
    # visits after the ppts discontinue infusions
    id_discontInfu <- unique(dat_long$ID[dat_long$f_discontInfu==1])
    dat_discont <- unique(dat_long[dat_long$ID%in%id_discontInfu,c("ID","day_discontInfu")])
    dat_discont <- data.frame(ID=rep(dat_discont$ID,each=6)
                              ,day_discontInfu=rep(dat_discont$day_discontInfu,each=6)
                              ,wk_post_enroll = rep(c("20wk","32wk","44wk","56wk","68wk","80wk"),nrow(dat_discont))
    )
    dat_discont$lw[dat_discont$wk_post_enroll=="20wk"] <- 99
    dat_discont$lw[dat_discont$wk_post_enroll=="32wk"] <- 183
    dat_discont$lw[dat_discont$wk_post_enroll=="44wk"] <- 267
    dat_discont$lw[dat_discont$wk_post_enroll=="56wk"] <- 351
    dat_discont$lw[dat_discont$wk_post_enroll=="68wk"] <- 435
    dat_discont$lw[dat_discont$wk_post_enroll=="80wk"] <- 519

    dat_discont$day_tar[dat_discont$wk_post_enroll=="20wk"] <- 140
    dat_discont$day_tar[dat_discont$wk_post_enroll=="32wk"] <- 224
    dat_discont$day_tar[dat_discont$wk_post_enroll=="44wk"] <- 308
    dat_discont$day_tar[dat_discont$wk_post_enroll=="56wk"] <- 392
    dat_discont$day_tar[dat_discont$wk_post_enroll=="68wk"] <- 476
    dat_discont$day_tar[dat_discont$wk_post_enroll=="80wk"] <- 560

    dat_discont$f_discont <- ifelse(dat_discont$day_discontInfu<dat_discont$lw,1,0)
    dat_discont <- subset(dat_discont,f_discont==1)

    win_type <- data.frame(t(rmultinom(nrow(dat_discont), size = 1, prob=prob_disContInfu_Win)))
    names(win_type) <- c("win_type1","win_type2","win_type3")
    dat_discont <- cbind(dat_discont,win_type)
    dat_discont$win1 <- round(runif(n=nrow(dat_discont),min=-41,max=-14),0)
    dat_discont$win2 <- round(runif(n=nrow(dat_discont),min=-14,max=14),0)
    dat_discont$win3 <- round(runif(n=nrow(dat_discont),min=14,max=42),0)
    dat_discont$win[dat_discont$win_type1==1] <- dat_discont$win1[dat_discont$win_type1==1]
    dat_discont$win[dat_discont$win_type2==1] <- dat_discont$win2[dat_discont$win_type2==1]
    dat_discont$win[dat_discont$win_type3==1] <- dat_discont$win3[dat_discont$win_type3==1]
    dat_discont$day_vis <- dat_discont$day_tar+dat_discont$win
    #missing visit
    dat_discont$f_vis <- rbinom(n=nrow(dat_discont),size=1,prob=1-prob_misVis)
    dat_discont$day_vis <- ifelse(dat_discont$f_vis==0,NA,dat_discont$day_vis)
    #visits for ppt who discontinue infusions
    dat_discont <- subset(dat_discont,select=c(ID,wk_post_enroll,f_vis,day_vis))
    dat_discont <- rename(dat_discont,c("wk_post_enroll"="visit_n"))
    dat_discont$visit_n  <- as.character(dat_discont$visit_n)
    dat_discont$visit_n[dat_discont$visit_n=="20wk"] <- 72
    dat_discont$visit_n[dat_discont$visit_n=="32wk"] <- 73
    dat_discont$visit_n[dat_discont$visit_n=="44wk"] <- 74
    dat_discont$visit_n[dat_discont$visit_n=="56wk"] <- 75
    dat_discont$visit_n[dat_discont$visit_n=="68wk"] <- 76
    dat_discont$visit_n[dat_discont$visit_n=="80wk"] <- 77
    dat_discont$f_infu[!is.na(dat_discont$ID)] <- 0 #for situation that no ppts who discontinue infusions
    #combine
    dat_long$day_discontInfu <- NULL
    dat_long$f_discontInfu <- NULL
    dat_long <- rbind(dat_long,dat_discont)
    rm(dat_discont)
    dat_long <- subset(dat_long,f_vis==1)
    #------------------- drop out ----------------------
    if(r_dropout==0){
      r_dropout <- 0.000000001
    }
    time_dropout <- data.frame(ID=1:n_ppt, time_dropout=rexp(n=n_ppt,rate = r_dropout/365))
    dat_long <- merge(dat_long,time_dropout)
    dat_long$f_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,0, dat_long$f_vis)
    dat_long$day_vis <- ifelse(dat_long$day_vis>=dat_long$time_dropout,NA, dat_long$day_vis)
    dat_long <- dat_long[order(as.numeric(dat_long$visit_n)),]
    #final data
    dat_long <- subset(dat_long,f_vis==1)
    dat_long$time_dropout <- NULL

    #------- generate dose levels --------
    dat_m <- data.frame(ID=1:n_male
                        ,dose=c(rep(10,n_male/3),rep(30,n_male/3),rep(0,n_male/3))
    )
    dat_f <- data.frame(ID=(n_male+1):(n_female+n_male)
                        ,dose=c(rep(10,n_female/3),rep(30,n_female/3),rep(0,n_female/3))
    )
    dat_cov <- rbind(dat_m,dat_f)
    # browser()
    rm(dat_m,dat_f)
    dat_long <- merge(dat_cov,dat_long,by="ID")
    rm(dat_cov)
    #---------------- infection status ------------------
    #survival time is simulated by using the single dose approach
    dat_long <- ddply(dat_long,.(dose),function(df){
      if(unique(df$dose)==0){
        t_event <- data.frame(ID=1:n_ppt,t_event=rexp(n=n_ppt,rate=r_incidence/365))
        df <- merge(df,t_event,by="ID")
        df$f_event <- ifelse(df$day_vis>=df$t_event,1,0)
        # browser()
        df <- ddply(df,.(ID),function(dat){
          # browser()
          if(all(dat$f_event==0)){
            dat$t_event <- max(dat$day_vis)
          }
          return(dat)
        })
        return(df)
      }else{
        if(unique(df$dose)==10){
          f_infect_singleDose(df,beta.t.10, ts=57)
        }else{
          f_infect_singleDose(df,beta.t.30, ts=81)
        }
      }
    })
    # browser()
    #---------- data preparation for KM --------------
    ID_cases <- unique(dat_long$ID[dat_long$f_event==1])
    dat_long$inf <- ifelse(dat_long$ID%in%ID_cases,1,0)
    dat_sur <- ddply(dat_long,.(ID,dose,inf),summarise,survtime=unique(t_event))
    dat_sur$adherence <- folder
    rm(dat_long)
   
    return(dat_sur)
    })
    res <- do.call(rbind,res_list)

    #-------------- final results -----------------------
    # browser()
    res_bind <- rbind(res,res_bind)
    return(res_bind)
   
}


###########################################################
# create different scenarios and run the simulation
###########################################################
n_ppt <- list(c(2250*100,2250*100))
adherence <- list(c(0.02,0.03,0.03,0.05)
                 ,c(0.1,0.15,0.15,0.15)
)
beta.t <- list(c(0.03,0.03))
scenario <- expand.grid(n_ppt,adherence,beta.t)
scenario <- apply(scenario,2,function(x)do.call("rbind",x))
scenario <- do.call("cbind",scenario)
scenario <- data.frame(scenario)
names(scenario) <- c("n_male","n_female","prob_misInfu","prob_misVis","r_infuDiscont","r_dropout","beta.t.30","beta.t.10")
# print(scenario)

res_bind <- NULL
for(i in 1:nrow(scenario)){
res_bind <- sim(n_male=scenario$n_male[i]
                  ,n_female=scenario$n_female[i]
                  ,prob_misInfu=scenario$prob_misInfu[i]
                  ,prob_misVis=scenario$prob_misVis[i]
                  ,r_infuDiscont=scenario$r_infuDiscont[i]
                  ,r_dropout=scenario$r_dropout[i]
                  ,beta.t.30=scenario$beta.t.30[i]
                  ,beta.t.10=scenario$beta.t.10[i])
                
}

###########################################################
# Draw KM curves as shown in Figure 2b (Huang et al., 2018)
###########################################################
dat_fit <- ddply(res_bind,.(adherence),function(df){
  # browser()
  fit <- survfit(Surv(survtime, inf) ~ dose, data = df)
  dat_fit <- data.frame(time    = fit$time
                        ,n.risk  = fit$n.risk
                        ,n.event = fit$n.event
                        ,surv = fit$surv
                        ,dose = gsub("dose=", "", summary(fit, censored = T)$strata))
  dat_fit$adherence <- unique(df$adherence)
  return(dat_fit)
  
})

dat_km_plot <- dplyr::mutate(dat_fit
                      ,adherence=ifelse(grepl("p1_0.1",adherence),"Medium adherence","High adherence")
                      ,dose = mapvalues(dose,c(0,10,30),c("placebo","10 mg/kg VRC01","30 mg/kg VRC01"))
)
ggplot(dat_km_plot,aes(x=time,y=1-surv,color=dose,linetype=dose))+
  geom_line(size=1.2)+
  xlab("Time since first infusion (weeks)")+
  ylab("Probability of HIV infection")+
  scale_color_manual("",values=c("black","blue","blue"))+
  scale_linetype_manual("",values=c(3,2,1))+
  scale_x_continuous(breaks=c((0:11)*8*7),labels=c((0:11)*8)
                     # ,limits = c(0,560)
  )+
  facet_grid(~adherence)+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom"
        ,plot.margin = unit(c(1,2,1,2),"lines")
  )

###########################################
# Draw  Probability of HIV infection within each infusion interval following
# ten 8-weekly IV infusions of VRC01 as shown in Figure 3 (Huang et al., 2018)
###########################################
T.beta.0.01 <- gen_survival_time(beta.t=0.01)
T.beta.0.03 <- gen_survival_time(beta.t=0.03)

gen_prob_inf <- function(dat_survival_time){
  dat <- data.frame(time=dat_survival_time)
  dat <- mutate(dat
                ,inf = ifelse(time<=560, 1, 0)
                ,time = ifelse(time>560,560,time)
  )
  dat_fit_list <- lapply(0:9, function(x){
    # browser()
    dat_sub <- dat
    dat_sub$time <- dat_sub$time-56*x
    dat_sub <- subset(dat_sub,time>0)
    fit <- survfit(Surv(time,inf) ~ 1,type="fleming",data=dat_sub)
    dat_fit <- data.frame(time = fit$time
                          ,n.risk  = fit$n.risk
                          ,n.event = fit$n.event
                          ,surv = fit$surv)
    dat_fit$infusion <- x+1
    dat_fit <- subset(dat_fit, time<=56)
    return(dat_fit)
    
  })
  
  dat_fit <- do.call(rbind,dat_fit_list)  
  dat_fit <- dplyr::mutate(dat_fit
                    ,prob_inf = 1 - surv
                    ,TIME = (infusion-1)*56+time)
  return(dat_fit)
}

dat.fit.beta.0.01 <- gen_prob_inf(T.beta.0.01)
dat.fit.beta.0.01$beta <- "0.01"
dat.fit.beta.0.03 <- gen_prob_inf(T.beta.0.03)
dat.fit.beta.0.03$beta <- "0.03"
dat_fit <- rbind(dat.fit.beta.0.01,dat.fit.beta.0.03)
dat_fit <- ddply(dat_fit,.(beta,infusion),transform,cum.n.event=cumsum(n.event))

#start the figure 
dat_text <- dplyr::summarise(dplyr::group_by(dat_fit,beta,infusion)
                      ,n.risk=round(min(n.risk)/1000,0)
                      ,cum.n.event=round(max(cum.n.event)/1000,0))
dat_text_time0 <- data.frame(beta=c("0.01","0.03"),infusion=0,n.risk=3000,cum.n.event=0)
dat_text <- rbind(data.frame(dat_text),dat_text_time0)
dat_text$TIME <- dat_text$infusion*56

p <- ggplot(dat_fit,aes(x=TIME,y=-log(surv),color=beta,linetype=beta))+
  geom_line(size=1.2)+
  geom_vline(xintercept = c(1:9*56),linetype="dashed")+
  # geom_blank()+
  scale_x_continuous(breaks=0:10*56,labels = 0:10*8
                     ,expand = c(0.05, 0.05)
  )+
  scale_color_manual(bquote(beta),values=c("red","blue"))+
  scale_linetype_manual(bquote(beta),values = c(1,3))+
  ylim(c(0,0.01))+
  xlab("Time since first infusion (weeks)")+
  ylab("Cumulative hazard of HIV infection till next infusion")+
  
  geom_text(data=subset(dat_text,beta=="0.01"&infusion!=10)
            ,aes(label=infusion+1, x=TIME, y=-Inf)
            ,size=4
            ,color="black"
            ,vjust=6)+
  geom_text(data=subset(dat_text,infusion==1), aes(x=-Inf,y=-Inf,label="beta[t]==0.01")
            ,parse=TRUE
            ,size=4.3
            ,color="black"
            ,vjust=7.5
            ,hjust=1
  )+
  geom_text(data=subset(dat_text,beta=="0.01")
            ,aes(label=n.risk, x=TIME, y=-Inf)
            ,size=4
            ,color="black"
            ,vjust=11)+
  geom_text(data=subset(dat_text,infusion==1), aes(x=-Inf,y=-Inf,label="beta[t]==0.03")
            ,parse=TRUE
            ,size=4.3
            ,color="black"
            ,vjust=9
            ,hjust=1
  )+
  geom_text(data=subset(dat_text,beta=="0.03")
            ,aes(label=n.risk, x=TIME, y=-Inf)
            ,size=4
            ,color="black"
            ,vjust=13)+
  
  geom_text(data=subset(dat_text,infusion==1), aes(x=-Inf,y=-Inf,label="beta[t]==0.01")
            ,parse=TRUE
            ,size=4.3
            ,color="black"
            ,vjust=12.5
            ,hjust=1
  )+
  geom_text(data=subset(dat_text,beta=="0.01")
            ,aes(label=cum.n.event, x=TIME, y=-Inf)
            ,size=4
            ,color="black"
            ,vjust=18)+
  geom_text(data=subset(dat_text,infusion==1), aes(x=-Inf,y=-Inf,label="beta[t]==0.03")
            ,parse=TRUE
            ,size=4.3
            ,color="black"
            ,vjust=14
            ,hjust=1
  )+
  geom_text(data=subset(dat_text,beta=="0.03")
            ,aes(label=cum.n.event, x=TIME, y=-Inf)
            ,size=4
            ,color="black"
            ,vjust=20)+
  theme_bw()+
  theme(plot.margin=unit(c(1,1,6,1), "cm"))

gt=ggplot_gtable(ggplot_build(p))
gt$layout$clip <- "off"
grid.draw(gt)

grid.text("Infusion:",gp = gpar(fontsize = 14,col="black",fontface="bold"),x=unit(0.02,"npc"),y=unit(0.27,"npc"),just="left")
grid.text("No. at risk:",gp = gpar(fontsize = 14,col="black",fontface="bold"),x=unit(0.02,"npc"),y=unit(0.24,"npc"),just="left")
grid.text("Cumulative No. of infections within each infusion interval:",gp = gpar(fontsize = 14,col="black",fontface="bold"),x=unit(0.02,"npc"),y=unit(0.14,"npc"),just="left")


