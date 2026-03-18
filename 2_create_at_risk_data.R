#1_create_at_risk_data.R
source('99_packages.R')

#load datasets fro 0_read_data.R
load('processed data/hai_data_with_sap_rules.rda')
load('processed data/cluster_time_mapping.rda')

#create dataset for key dates per patient to define time at risk 
## based on ward admission and discharge dates
dates_atrisk <- dat %>% 
  select(record_id,cluster,date_admission_hospital,date_admission_ward,date_discharge_ward,date_discharge_hospital,reasonfordischargefromward) %>%
  rowwise() %>% mutate(date_daily=list(seq.Date(date_admission_hospital,date_discharge_ward+days(2),'1 day'))) %>% ungroup() #chenages from date_discharge_hospital to date_discharge_ward + days(2) to align with time at risk window
  
#unnest to get all dates per record_id - still rowwise
#state: whether the patient is in the hospital or a study ward each day
#discharge: 1 if date_daily==date_discharge_hospital, 0 otherwise
dates_atrisk_long <- dates_atrisk %>% unnest(cols=date_daily) %>% 
  mutate(state = case_when(
    between(date_daily,date_admission_ward,date_discharge_ward) ~ 'Study ward',TRUE~'Hospital'),
    discharged = if_else(date_daily>=date_discharge_ward,1,0)
    ) 

#add date_hai - join with date_daily 
#defined atrisk days on the ward
dat_hai_all <- dat %>% filter(rowSums(across(c(date_hai_onset1, date_hai_onset2, date_hai_onset3), ~ !is.na(.x)))>0) %>% select(record_id,cluster,hai1,hai2,hai3,starts_with('date_hai')) %>%
  pivot_longer(cols=c(hai1,hai2,hai3,date_hai_onset1, date_hai_onset2, date_hai_onset3),names_to=c(".value","hai_number"),names_pattern="(.)(\\d{1})") %>% rename('hai'=i,'date_hai_onset'=t) %>% drop_na() %>% select(-hai_number)


#collapse multiple infections for the same person with the same date of onset
dat_hai_all<-dat_hai_all %>% reframe(hai = str_c(hai,collapse = ';'),.by=c(record_id,cluster,date_hai_onset)) 

#add flags for hap, luri, eent_oral
dat_hai_all<- dat_hai_all %>% mutate(hap = if_else(grepl('PN',hai),1,0),luri = if_else(grepl('LRI|EENT-UR',hai),1,0),eent_oral=if_else(grepl('EENT-ORAL',hai),1,0)) %>% select(-c(hai)) #remove now redundant columns

dates_atrisk_long_1 <- dates_atrisk_long %>% left_join(dat_hai_all,by=join_by(record_id,cluster,date_daily==date_hai_onset)) %>% mutate(across(c(hap,luri,eent_oral),~replace_na(.x,0))) %>% mutate(hai=pmin(1,hap+luri+eent_oral))

#check for excluded infections - should be 0
anti_join(dat_hai_all,dates_atrisk_long_1,by = join_by(record_id, cluster, hap, luri, eent_oral))

dates_atrisk_long<-dates_atrisk_long_1


#add trial dates
dates_atrisk_long <- dates_atrisk_long %>% left_join(cluster_time_mapping %>% select(cluster,date_daily,study_week,intervention,date_control_start,date_intervention_end),by = join_by(cluster,date_daily))

dates_atrisk_long<-dates_atrisk_long %>% fill(c(date_control_start,date_intervention_end),.direction='up') 

#add ward day as a time scale
dates_atrisk_long <- dates_atrisk_long %>% mutate(on_ward = if_else(state=='Study ward',1,0)) %>% mutate(ward_day = cumsum(on_ward),.by=record_id)
dates_atrisk_long <- filter(dates_atrisk_long,!is.na(intervention)) %>% arrange(record_id,date_daily)

#SAP control/intervention attribution - ID 2878 (LURI infection); update intervention = 0
dates_atrisk_long <- dates_atrisk_long %>% mutate_at('intervention',~if_else(record_id==2878 & luri==1,0,.))

#add start and stop dates, at risk days
dates_atrisk_long <- dates_atrisk_long %>% 
  mutate(date_atrisk_start = pmax(date_control_start,date_admission_ward),
         date_atrisk_stop = pmin(date_discharge_ward+days(2),date_intervention_end+days(2))) %>% #+2 days post ward discharge or intervention end, whichever happens first
  mutate(atrisk=if_else(between(date_daily,date_atrisk_start,date_atrisk_stop),1,0),
         atrisk_day = pmax(0,as.numeric(date_daily - date_atrisk_start,'days'))+1) 

#filter to at risk timeframe only - excludes time off-ward
dates_atrisk_long <- filter(dates_atrisk_long,atrisk==1)

#define stopping times - infection onset or the end of the at-risk window; whichever comes first

#primary outcome - HAP
ad1<-dates_atrisk_long %>% filter(hap==1) %>% slice_min(atrisk_day,by=record_id)  #date of first hap
ad0<-dates_atrisk_long %>% filter(!record_id %in% ad1$record_id) %>% slice_max(atrisk_day,by=record_id) 
atrisk_times_hap<-bind_rows(ad0,ad1) %>% select(record_id,atrisk_day) %>% rename(stop_hap=atrisk_day) %>% arrange(record_id)

dat_atrisk_hap <- dates_atrisk_long %>% left_join(atrisk_times_hap,by=join_by(record_id)) %>% 
  filter(between(atrisk_day,1,stop_hap)) %>% select(record_id,cluster,date_daily,ward_day,atrisk_day,hap,discharged) %>% left_join(cluster_time_mapping %>% select(hospital,ward,cluster,date_daily,study_week,study_period,study_month_number,intervention),by = join_by(cluster,date_daily))

#define end status/event at end of atrisk time
dat_atrisk_hap<-dat_atrisk_hap %>% mutate(event = case_when(hap==0 & discharged==1 ~ 2,hap==1 ~ 1,TRUE~0))

#add hap flag to dat
idx_hap<-dat_atrisk_hap %>% filter(hap==1) %>% distinct(record_id) %>% pull()
dat <- dat %>% mutate(hap = if_else(record_id %in% idx_hap,1,0))

#luri
ad1<-dates_atrisk_long %>% filter(luri==1) %>% slice_min(atrisk_day,by=record_id)  #date of first luri
ad0<-dates_atrisk_long %>% filter(!record_id %in% ad1$record_id) %>% slice_max(atrisk_day,by=record_id) 
atrisk_times_luri<-bind_rows(ad0,ad1) %>% select(record_id,atrisk_day) %>% rename(stop_luri=atrisk_day) %>% arrange(record_id)

dat_atrisk_luri <- dates_atrisk_long %>% left_join(atrisk_times_luri,by=join_by(record_id)) %>% 
  filter(between(atrisk_day,1,stop_luri)) %>% select(record_id,cluster,date_daily,ward_day,atrisk_day,luri,discharged) %>% left_join(cluster_time_mapping %>% select(hospital,ward,cluster,date_daily,study_week,study_period,study_month_number,intervention),by = join_by(cluster,date_daily))

#define end status/event at end of atrisk time
dat_atrisk_luri <-dat_atrisk_luri %>% mutate(event = case_when(luri==0 & discharged==1 ~ 2,luri==1 ~ 1,TRUE~0))

#add luri flag to dat - sap rule attribution due to control/intervention overlap
idx_luri<-dat_atrisk_luri %>% filter(luri==1) %>% distinct(record_id) %>% pull()
dat <- dat %>% mutate(luri = if_else(record_id %in% idx_luri,1,0))


#eent-oral
ad1<-dates_atrisk_long %>% filter(eent_oral==1) %>% slice_min(atrisk_day,by=record_id)  #date of first eent oral
ad0<-dates_atrisk_long %>% filter(!record_id %in% ad1$record_id) %>% slice_max(atrisk_day,by=record_id) 
atrisk_times_eent_oral<-bind_rows(ad0,ad1) %>% select(record_id,atrisk_day) %>% rename(stop_eent_oral=atrisk_day) %>% arrange(record_id) 

dat_atrisk_eent_oral <- dates_atrisk_long %>% left_join(atrisk_times_eent_oral,by=join_by(record_id)) %>% 
  filter(between(atrisk_day,1,stop_eent_oral)) %>% select(record_id,cluster,date_daily,ward_day,atrisk_day,eent_oral,discharged) %>% left_join(cluster_time_mapping %>% select(hospital,ward,cluster,date_daily,study_week,study_period,study_month_number,intervention),by = join_by(cluster,date_daily))

#define end status/event at end of atrisk time
dat_atrisk_eent_oral <-dat_atrisk_eent_oral %>% mutate(event = case_when(eent_oral==0 & discharged==1 ~ 2,eent_oral==1 ~ 1,TRUE~0))

#add eent-oral flag to dat
idx_eent_oral<-dat_atrisk_eent_oral %>% filter(eent_oral==1) %>% distinct(record_id) %>% pull()
dat <- dat %>% mutate(eent_oral = if_else(record_id %in% idx_eent_oral,1,0))


#one or more hais - teritiary outcome; NB: hai = 1 corresponds to one or more hap based on first date of onset
ad1<-dates_atrisk_long %>% filter(hai==1) %>% slice_min(atrisk_day,by=record_id)  #date of first hai
ad0<-dates_atrisk_long %>% filter(!record_id %in% ad1$record_id) %>% slice_max(atrisk_day,by=record_id) #no hai during at-risk window
atrisk_times_hai<-bind_rows(ad0,ad1) %>% select(record_id,atrisk_day) %>% rename(stop_hai=atrisk_day) %>% arrange(record_id)

dat_atrisk_hai<- dates_atrisk_long %>% left_join(atrisk_times_hai,by=join_by(record_id)) %>% 
  filter(between(atrisk_day,1,stop_hai)) %>% select(record_id,cluster,date_daily,ward_day,atrisk_day,hai,discharged) %>% left_join(cluster_time_mapping %>% select(hospital,ward,cluster,date_daily,study_week,study_period,study_month_number,intervention),by = join_by(cluster,date_daily))

#define end status/event at end of atrisk time
dat_atrisk_hai<-dat_atrisk_hai %>% mutate(event = case_when(hai==0 & discharged==1 ~ 2,hai==1 ~ 1,TRUE~0))

#add hai flag to dat
idx_hai<-dat_atrisk_hai %>% filter(hai==1) %>% distinct(record_id) %>% pull()
dat <- dat %>% mutate(hai = if_else(record_id %in% idx_hai,1,0))


#define dataset for time on ward
dat_atrisk <- dates_atrisk_long %>% select(record_id,cluster,date_daily,ward_day,atrisk_day) %>% left_join(cluster_time_mapping %>% select(hospital,ward,cluster,date_daily,study_week,study_period,study_month_number,intervention),by = join_by(cluster,date_daily))

save(dat_atrisk,dat_atrisk_hai,dat_atrisk_hap,dat_atrisk_luri,dat_atrisk_eent_oral,N,dat,dat_sap_ii,dat_sap_iv,dat_lt48hrs,file='processed data/trial_data.rda')

