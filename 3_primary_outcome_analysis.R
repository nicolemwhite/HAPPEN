#3_primary_outcome_analysis.R
source('99_packages.R')
source('99_functions.R')
load('processed data/cluster_time_mapping.rda')
load('processed data/trial_data.rda')


#patient-level - take event on the stop day - date of first hai onset or discharge
dat.mod<-dat_atrisk_hap %>% slice_max(atrisk_day,by=record_id) %>% mutate_at(c('study_period','study_month_number'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))

mod.fit_hap<-NULL
mod.fit_hap[['patient:cloglog']]<-glmmTMB(hap~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_hap[['patient:logit']]<-glmmTMB(hap~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

#sensitivity
mod.fit_hap[['patient:cloglog - sens']]<-glmmTMB(hap~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_hap[['patient:logit - sens']]<-glmmTMB(hap~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

cluster_vec<-as.character(1:9)

#loo - note use of study_period over study_month_number for convergence
mod.fit_hap[['patient:cloglog - loo']]<-lapply(cluster_vec,function(h) glmmTMB(hap~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))
mod.fit_hap[['patient:logit - loo']]<-lapply(cluster_vec,function(h) glmmTMB(hap~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))

#cluster level - total admissions per month (ignores time spent on the ward)
dat.mod_c <- dat_atrisk_hap %>% reframe(n_admissions = length(unique(record_id)),n_admission_days=n(),n_hap=sum(hap),.by=c(cluster,hospital,study_month_number,study_period,intervention)) %>% mutate_at('study_month_number',~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))
mod.fit_hap[['cluster (total admissions):log']]<-glmmTMB(n_hap~(1|cluster)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_hap[['cluster (total admissions):log - sens']]<-glmmTMB(n_hap~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

mod.fit_hap[['cluster (total admission days):log']]<-glmmTMB(n_hap~(1|cluster)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_hap[['cluster (total admission days):log - sens']]<-glmmTMB(n_hap~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))


#loo - study_period used in place of study_month_number for convergence
mod.fit_hap[['cluster (total admissions):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_hap~(1|cluster)+intervention+study_period+offset(log(n_admissions/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))
mod.fit_hap[['cluster (total admission days):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_hap~(1|cluster)+intervention+study_period+offset(log(n_admission_days/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))

mod.pred_hap<-NULL
mod.pred_hap[['cluster (total admissions):log']]<-predict_response(mod.fit_hap[['cluster (total admissions):log']],'intervention',condition=list(n_admissions=100),margin = 'empirical') %>% data.frame()
mod.pred_hap[['cluster (total admission days):log']]<-predict_response(mod.fit_hap[['cluster (total admission days):log']],'intervention',condition=list(n_admission_days=100),margin = 'empirical') %>% data.frame()

save(mod.fit_hap,cluster_vec,mod.pred_hap,file='model fits/hap_model_results_r1.rda') #note R1 run on 20 Feb 2026 after final data checks

