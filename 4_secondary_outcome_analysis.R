#4_secondary_outcome_analysis.R
source('99_packages.R')
source('99_functions.R')
load('processed data/cluster_time_mapping.rda')
load('processed data/trial_data.rda')

#LURI
#patient-level - take event on the stop day - date of first hai onset or discharge
dat.mod<-dat_atrisk_luri %>% slice_max(atrisk_day,by=record_id) %>% mutate_at(c('study_period','study_month_number'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))

mod.fit_luri<-NULL
mod.fit_luri[['patient:cloglog']]<-glmmTMB(luri~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'))
mod.fit_luri[['patient:logit']]<-glmmTMB(luri~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'))

#sensitivity
mod.fit_luri[['patient:cloglog - sens']]<-glmmTMB(luri~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_luri[['patient:logit - sens']]<-glmmTMB(luri~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

cluster_vec<-as.character(1:9)

#note study_period used for convergence here
mod.fit_luri[['patient:cloglog - loo']]<-lapply(cluster_vec,function(h) glmmTMB(luri~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))
mod.fit_luri[['patient:logit - loo']]<-lapply(cluster_vec,function(h) glmmTMB(luri~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))

#cluster level - total admissions per month (ignores time spent on the ward)
dat.mod_c <- dat_atrisk_luri %>% reframe(n_admissions = length(unique(record_id)),n_admission_days=n(),n_luri=sum(luri),.by=c(cluster,hospital,study_month_number,study_period,intervention)) %>% mutate_at(c('study_month_number','study_period'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))
mod.fit_luri[['cluster (total admissions):log']]<-glmmTMB(n_luri~(1|cluster)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_luri[['cluster (total admission days):log']]<-glmmTMB(n_luri~(1|cluster)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

mod.fit_luri[['cluster (total admissions):log - sens']]<-glmmTMB(n_luri~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_luri[['cluster (total admission days):log - sens']]<-glmmTMB(n_luri~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

mod.fit_luri[['cluster (total admissions):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_luri~(1|cluster)+intervention+study_period+offset(log(n_admissions/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))
mod.fit_luri[['cluster (total admission days):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_luri~(1|cluster)+intervention+study_period+offset(log(n_admission_days/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))


mod.pred_luri<-NULL
mod.pred_luri[['cluster (total admissions):log']]<-predict_response(mod.fit_luri[['cluster (total admissions):log']],'intervention',condition=list(n_admissions=100),margin = 'empirical') %>% data.frame()
mod.pred_luri[['cluster (total admission days):log']]<-predict_response(mod.fit_luri[['cluster (total admission days):log']],'intervention',condition=list(n_admission_days=100),margin = 'empirical') %>% data.frame()

save(mod.fit_luri,cluster_vec,mod.pred_luri,file='model fits/luri_model_results_r1.rda') #after final data checks 20 feb 2026

#EENT-ORAL
#patient-level - take event on the stop day - date of first hai onset or discharge
dat.mod<-dat_atrisk_eent_oral %>% slice_max(atrisk_day,by=record_id) %>% mutate_at(c('study_period','study_month_number'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))

mod.fit_eent_oral<-NULL
mod.fit_eent_oral[['patient:cloglog']]<-glmmTMB(eent_oral~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_eent_oral[['patient:logit']]<-glmmTMB(eent_oral~(1|cluster)+intervention+study_month_number+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

#sensitivity
mod.fit_eent_oral[['patient:cloglog - sens']]<-glmmTMB(eent_oral~(1|cluster)+(0+intervention|hospital)+intervention+study_period+offset(log(atrisk_day)),data=dat.mod,family=binomial('cloglog'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_eent_oral[['patient:logit - sens']]<-glmmTMB(eent_oral~(1|cluster)+(0+intervention|hospital)+intervention+study_period+offset(log(atrisk_day)),data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

#loo - study_period in place of study_month_number; change optimiser to address convergence issues
mod.fit_eent_oral[['patient:cloglog - loo']]<-lapply(cluster_vec,function(h) glmmTMB(eent_oral~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('cloglog')))
mod.fit_eent_oral[['patient:logit - loo']]<-lapply(cluster_vec,function(h) glmmTMB(eent_oral~(1|cluster)+intervention+study_period+offset(log(atrisk_day)),data=filter(dat.mod,cluster!=h),family=binomial('logit')))

#cluster level - total admissions per month (ignores time spent on the ward)
dat.mod_c <- dat_atrisk_eent_oral %>% reframe(n_admissions = length(unique(record_id)),n_admission_days=n(),n_eent_oral=sum(eent_oral),.by=c(cluster,hospital,study_month_number,study_period,intervention)) %>% mutate_at(c('study_period','study_month_number'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))
mod.fit_eent_oral[['cluster (total admissions):log']]<-glmmTMB(n_eent_oral~(1|cluster)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_eent_oral[['cluster (total admissions):log - sens']]<-glmmTMB(n_eent_oral~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

mod.fit_eent_oral[['cluster (total admission days):log']]<-glmmTMB(n_eent_oral~(1|cluster)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_eent_oral[['cluster (total admission days):log - sens']]<-glmmTMB(n_eent_oral~(1|cluster)+(0+intervention|hospital)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

#loo - study_period in place of study_month_number
mod.fit_eent_oral[['cluster (total admissions):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_eent_oral~(1|cluster)+intervention+study_period+offset(log(n_admissions/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log')))#,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))
mod.fit_eent_oral[['cluster (total admission days):log - loo']]<-lapply(cluster_vec, function(h) glmmTMB(n_eent_oral~(1|cluster)+intervention+study_period+offset(log(n_admission_days/100)),data=filter(dat.mod_c,cluster!=h),family=poisson('log')))#,control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS"))))


mod.pred_eent_oral<-NULL
mod.pred_eent_oral[['cluster (total admissions):log']]<-predict_response(mod.fit_eent_oral[['cluster (total admissions):log']],'intervention',condition=list(n_admissions=100),margin = 'empirical') %>% data.frame()
mod.pred_eent_oral[['cluster (total admission days):log']]<-predict_response(mod.fit_eent_oral[['cluster (total admission days):log']],'intervention',condition=list(n_admission_days=100),margin = 'empirical') %>% data.frame()


save(mod.fit_eent_oral,cluster_vec,mod.pred_eent_oral,file='model fits/eent_oral_model_results_r1.rda') #after final data checks, 20.02.26
