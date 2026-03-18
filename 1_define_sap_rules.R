#1_apply_sap_rules.R

#TODO SAP rule checks: 
##record_id: 2878 (rule ii) - intervention = 0

source('99_packages.R')

load('processed data/hai_data.rda')
load('processed data/cluster_time_mapping.rda')

#start with general inclusion/exclusion criteria
#for patients who were on the study ward for more than 48hrs, exclude records with missing ward admission and/or ward discharge dates

#exclude patients who were not on the study ward for more than 48hr
N[['On study ward for <48hr']]<-filter(dat,on_ward_gt48hrs=='No') %>% nrow()

dat_lt48hrs <- filter(dat,on_ward_gt48hrs=='No')

dat<-dat %>% filter(on_ward_gt48hrs=='Yes')

N[['No medical record available - excluded for other reasons']]<-dat %>% filter(exclude_reason=='Other') %>% nrow()
dat <- dat %>% filter(is.na(exclude_reason))

 
N[['Missing ward admission and/or discharge dates']] <- dat %>% filter(is.na(date_admission_ward)|is.na(date_discharge_ward)) %>% nrow()
exclude_idx<- dat %>% filter(is.na(date_admission_ward)|is.na(date_discharge_ward)) %>% pull(record_id)
 
 
dat <- dat %>% filter(!record_id %in% exclude_idx)
 

#check date order independent of infection outcome - should be hospital admit (1), ward admit (2), ward discharge (3), hospital_discharge (4)
dat_event_dates<-dat %>% select(record_id,date_admission_hospital,date_admission_ward,date_discharge_ward,date_discharge_hospital) %>% 
  pivot_longer(cols=starts_with('date_'),names_to = 'event',values_to='date_event') %>% drop_na() %>% 
  group_by(record_id) %>% arrange(date_event,.by_group=T) %>% ungroup() %>%
  mutate(event_code = case_match(event,'date_admission_hospital'~1,'date_admission_ward'~2,'date_discharge_ward'~3,'date_discharge_hospital'~4)) %>% reframe(event_order=paste0(event_code,collapse=''),.by=record_id)

#check outlying records based on date order
check_idx <- filter(dat_event_dates,event_order!=1234)
N[['Discharged before admitted']]<-nrow(check_idx)

#add trial dates based on date_admission_ward
dat<-dat %>% left_join(trial_dates %>% select(cluster,starts_with('date_')) %>% mutate(across(starts_with("date_"), ~ as.Date(format(.x, "%Y-%m-%d")))),by = join_by(cluster))

#discharged before control start
N[['Discharged from study ward before control start date']]<-dat %>% filter(date_discharge_ward<date_control_start) %>% nrow()

dat<-dat %>% filter(date_discharge_ward>=date_control_start)

#code hais
#study population - 
##For patients who meet the inclusion criteria, confirmed infections will be
##excluded if acquired on the day of study ward admission or within the first 48 hours of admission.
### time at risk - exclude if infection is more than 48 hrs after ward discharge

dat_hai_all<-dat %>% filter(rowSums(across(c(date_hai_onset1, date_hai_onset2, date_hai_onset3), ~ !is.na(.x)))>0) %>% select(record_id,cluster,date_admission_ward,date_discharge_ward,reasonfordischargefromward,hai1,hai2,hai3,starts_with('date_hai')) %>%
  pivot_longer(cols=c(hai1,hai2,hai3,date_hai_onset1, date_hai_onset2, date_hai_onset3),names_to=c(".value","hai_number"),names_pattern="(.)(\\d{1})") %>% rename('hai'=i,'date_hai_onset'=t) %>% drop_na() 

N[['Excluded infections: Infections acquired before study ward admission']]<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_admission_ward,'days')<0) %>% nrow()
N[['Excluded infections: Infection acquired on day or study ward admission of within the first 48 hours of admission']]<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_admission_ward,'days') %in% 0:1) %>% nrow()
N[['Excluded infections: Infection acquired more than 2 days after ward discharge']]<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_discharge_ward,'days')>=2) %>% nrow()


#apply study population exclusions
dat_hai_all<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_admission_ward,'days')>=2,as.numeric(date_hai_onset - date_discharge_ward)<2) 

#SAP rules leading to exclusion
#i If a patient is admitted to a study ward before the trial start date and the date of
#infection onset is less than 2 days after the trial start date, the infection will be excluded.
#add trial dates
dat_hai_all<-dat_hai_all %>% left_join(trial_dates %>% select(cluster,starts_with('date_')) %>% mutate(across(starts_with("date_"), ~ as.Date(format(.x, "%Y-%m-%d")))),by = join_by(cluster))

N[['Excluded infections: Infection acquired less than 2 days after control start']]<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_control_start,'days')<2) %>% nrow()
exclude_idx<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_control_start,'days')<2) %>% pull(record_id)

dat_hai_all<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_control_start,'days')>=2) 

#iiib If a patient is admitted to a study ward under the intervention condition, and acquires an infection in the same ward after the trial end date
#The infection will be excluded if the date of onset is more than 2 days after the trial 
N[['Excluded infections: Infection acquired 2+ days after intervention end']]<-dat_hai_all %>% filter(as.numeric(date_hai_onset-date_intervention_end,'days')>=2) %>% nrow()
exclude_idx<-c(exclude_idx,dat_hai_all %>% filter(as.numeric(date_hai_onset-date_intervention_end,'days')>=2) %>% pull(record_id))
#iiia as above but if within two days of intervention start date - keep in intervention condition
dat_hai_all<-dat_hai_all %>% filter(!as.numeric(date_hai_onset-date_intervention_end,'days')>=2) 


##########
#apply SAP rules ii and iv
dat_sap_ii<-dat_hai_all %>% filter(date_admission_ward<date_intervention_start & date_hai_onset>=date_intervention_start) %>% 
  select(record_id,date_admission_ward,date_intervention_start,contains('hai')) %>% mutate(intervention_to_hai = as.numeric(date_hai_onset - date_intervention_start,'days'))

#check transfers
ward_transfers <- dat %>% filter(reasonfordischargefromward=='Transferred to another ward') %>% pull(record_id)
dat_sap_iv <- dat_hai_all %>%  filter(date_hai_onset>=date_discharge_ward) %>% mutate(discharge_to_hai = as.numeric(date_hai_onset - date_discharge_ward,'days')) %>%
  select(record_id,date_admission_ward,date_discharge_ward,reasonfordischargefromward,contains('hai'))

# deprecated
# # #censoring based on at-risk timeframe
# #infection before ward admission
# dat_sap_preadmission<-dat_hai_all %>% filter(date_hai_onset<date_admission_ward) %>%
#   select(record_id,date_admission_ward,date_hai_onset,hai,date_admission_ward)
# 
# dat_sap_postdischarge<-dat_hai_all %>% filter(date_hai_onset>=(date_discharge_ward+days(2))) %>% 
#   select(record_id,date_admission_ward,date_discharge_ward,reasonfordischargefromward,date_hai_onset,hai)

#apply VAP criteria
ad_vap<-dat %>% select(record_id,cluster,vap,statedateofventilatorinsertion,statedateofventilatorremoval) %>% filter(record_id %in% dat_hai_all[['record_id']])

dat_hai_all <- dat_hai_all %>% left_join(ad_vap,by=join_by(record_id,cluster)) %>% mutate(across(starts_with('statedate'),~ymd(.))) %>% 
  mutate(exclude_vap = if_else(grepl('PN',hai) & grepl('Yes',vap) & date_hai_onset>statedateofventilatorinsertion & date_hai_onset<=(statedateofventilatorremoval+days(2)),1,0))

#exclude vaps
N[['Infections excluded: VAP']]<-dat_hai_all %>% filter(exclude_vap==1) %>% nrow()

dat_hai_all<-dat_hai_all %>% filter(exclude_vap==0)

#pivot_wider to original format following SAP exclusions
#add hai info back onto full dataset
ad_hai <-dat_hai_all %>% pivot_wider(names_from=hai_number,values_from=c(hai,date_hai_onset),names_sep = "") %>% mutate(across(starts_with('date_'),~ymd(.))) %>% mutate(date_hai = pmin(date_hai_onset1,date_hai_onset2,date_hai_onset3,na.rm=T)) %>% select(record_id,hai1:date_hai_onset3)

dat<-dat %>% select(-c(hai1,hai2,hai3,date_hai_onset1,date_hai_onset2,date_hai_onset3))
dat<- dat %>% left_join(ad_hai,by=join_by(record_id))

save(dat,N,dat_sap_ii,dat_sap_iv,dat_lt48hrs,file='processed data/hai_data_with_sap_rules.rda')
