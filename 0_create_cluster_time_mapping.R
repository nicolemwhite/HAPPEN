#0_create_clustser_time_mapping.R
source('99_packages.R')

#read in trial info
trial_dates<-xl.read.file('../../Data/Site info/trial dates.xlsx') %>% janitor::clean_names() %>% mutate_at('hospital',~as.character(.))
trial_dates <- trial_dates %>% 
  mutate(cluster = paste(hospital,ward,sep='_') %>% factor(levels=c('1_1','1_2','1_3','2_1','2_2','2_3','3_1','3_2','3_3'),labels=(1:9)) %>% as.character) %>%
  mutate_at(c('hospital','ward'),~as.character(.))


trial_start_date <- min(trial_dates$date_control_start) %>% ymd()

cluster_time_mapping <- trial_dates %>% rowwise() %>% mutate(date_daily = list(seq(date_control_start, date_intervention_end, by = "day"))) %>% 
  ungroup() %>% unnest(cols=date_daily) %>%
  mutate(across(starts_with("date_"), ~ as.Date(format(.x, "%Y-%m-%d")))) %>% 
  mutate(study_weekday=weekdays(date_daily),study_week = floor(as.numeric(date_daily - trial_start_date,'weeks'))+1,
         study_monthyear = paste(month(date_daily),year(date_daily),sep='_') %>% factor(.),
         intervention = if_else(date_daily<date_intervention_start,0,1))


#filter based on start and end dates
cluster_time_mapping<-filter(cluster_time_mapping,between(date_daily,date_control_start,date_intervention_end))


#add time steps as per planned study design + 1 extra period to account for delayed start in hospital 3

design_periods_start <- trial_start_date+weeks(seq(0,52,13))
design_periods <- tibble(start=design_periods_start ,stop=c(design_periods_start[-1]-days(1),"2025-08-24")) %>% add_column('period'=1:5,.before=1)

cluster_time_mapping <- cluster_time_mapping %>%
  mutate(
    idx = findInterval(date_daily, design_periods$start),            # index of most recent start
    # valid if idx > 0 and date_admission <= corresponding stop
    in_interval = idx > 0 & date_daily <= design_periods$stop[pmax(idx, 1)],
    study_period = ifelse(in_interval, design_periods$period[idx], NA_character_))

#add study periods defined as the first Monday of each month
trial_end_date <- as.Date("2025-08-24") #original end date

# first day of each month spanned by the range
month_starts <- seq(floor_date(trial_start_date, "month"),floor_date(trial_end_date,   "month"),by = "month")

# compute offset to the first Monday (with Monday as day 1); wday(..., week_start = 1): Monday=1, Tuesday=2, ..., Sunday=7
offset_days <- (8 - wday(month_starts, week_start = 1)) %% 7
first_mondays <- month_starts + days(offset_days)

# filter to the original range
first_mondays <- first_mondays[first_mondays >= trial_start_date & first_mondays <= trial_end_date]

periods_mths <- tibble(start=first_mondays ,stop=c(first_mondays[-1]-days(1),"2025-08-24")) %>% add_column('study_month_number'=1:length(first_mondays),.before=1)

cluster_time_mapping <- cluster_time_mapping %>%
  mutate(
    idx = findInterval(date_daily, periods_mths$start),  
    # valid if idx > 0 and date_admission <= corresponding stop
    in_interval = idx > 0 & date_daily <= periods_mths$stop[pmax(idx, 1)],
    study_month_number = ifelse(in_interval, periods_mths$study_month_number[idx], NA_character_))

save(trial_dates,cluster_time_mapping,file='processed data/cluster_time_mapping.rda')