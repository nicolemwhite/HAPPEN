#5_tertiary_outcome_analysis.R
source('99_packages.R')
source('99_functions.R')
load('processed data/hai_data_with_sap_rules.rda')
load('processed data/cluster_time_mapping.rda')
load('processed data/trial_data.rda')


#1. the proportion of patient who acquire one or more infections during their admission
#study period defined at discharge; logit model

dat.mod <- dat_atrisk_hai %>% slice_max(atrisk_day,by=record_id) %>% mutate_at(c('study_period','study_month_number'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))

mod.fit_hai<-NULL

mod.fit_hai[['patient: 1+ infections during admission']] <- glmmTMB(hai~(1|cluster)+intervention+study_month_number,data=dat.mod,family=binomial('logit'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))


dat.mod_c <- dat_atrisk_hai %>% reframe(n_admissions = length(unique(record_id)),n_admission_days=n(),n_hai=sum(hai),.by=c(cluster,hospital,study_month_number,study_period,intervention)) %>% mutate_at(c('study_month_number','study_period'),~factor(.,levels=unique(.)) %>% relevel(.,ref='1'))
mod.fit_hai[['cluster (total admissions):log']]<-glmmTMB(n_hai~(1|cluster)+intervention+study_month_number+offset(log(n_admissions/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))
mod.fit_hai[['cluster (total admission days):log']]<-glmmTMB(n_hai~(1|cluster)+intervention+study_month_number+offset(log(n_admission_days/100)),data=dat.mod_c,family=poisson('log'),control=glmmTMBControl(optimizer=optim,optArgs=list(method="BFGS")))

save(mod.fit_hai,file='model fits/hai_model_results_r1.rda')

