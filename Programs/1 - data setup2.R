rm(list=ls())
library(fanetc)
library(haven)

# County dictionary 
county_dic <- readxl::read_xlsx("../Data/Copy of ipums_usa_identified_counties.xlsx", skip=1)
county_dic <- county_dic[,1:5]
county_dic <- county_dic[,-4]
colnames(county_dic) <- c("state","county","statefip","countyfip")
county_dic <- county_dic %>% 
  mutate_at(c("statefip","countyfip"), as.numeric)

# Puma dictionary
puma_dic <- readxl::read_xls("../Data/puma_migpuma1_pwpuma00.xls", col_names = TRUE, skip = 1)
colnames(puma_dic) <- c("state","puma","migplac1","migpuma1")
puma_dic <- puma_dic %>% 
  mutate_all(as.numeric) 

raw_d <- read_dta('../Data/usa_00008.dta')

### Match county names to fips 
work_d <- raw_d %>% 
  left_join(county_dic %>%
              dplyr::select(county, countyfip, statefip), by = c("countyfip","statefip")) 

rm('raw_d')



### Data setup based on PUMAs 


reg_d <- work_d %>% filter(age!=0)

# > puma_dic %>% filter(migpuma1==1300 & state==36)
# # A tibble: 1 × 4
# state  puma migplac1 migpuma1
# <dbl> <dbl>    <dbl>    <dbl>
#   36  1300       36     1300

reg_d %>% 
  filter(statefip==36) %>% 
  filter(puma==1300 & age==27) %>% 
  filter(year==2015) %>%
  dplyr::select(year,age,puma,migpuma1, pernum, perwt, migrate1,countyfip,migcounty1) %>% 
  group_by(year) %>% 
  mutate(pop = sum(perwt*pernum)) %>% 
  mutate(samepuma = puma == migpuma1) %>% 
  mutate(migrate = NA, 
         migrate = ifelse(migrate1 == 1, 0, migrate),
         migrate = ifelse(migrate1 == 2 & !samepuma, 1, migrate), 
         migrate = ifelse(migrate1 == 2 & samepuma, 0, migrate), 
         migrate = ifelse(migrate1 %in% c(3,4), 1, migrate)) %>% 
  mutate(migrant = sum(pernum*perwt*migrate)) %>% 
  # slice(1) %>% 
  mutate(rate = migrant/pop) 




reg_d <- reg_d %>% 
  filter(year >= 2013) %>% 
  mutate(perwt = perwt*pernum)  ## experiment

reg_d %>% filter(age==27 & statefip == 36)

# Migration PUMA codes are state-dependent, so MIGPUMA1 must be combined with 
# state codes (see MIGPLAC1) to distinguish Migration PUMAs in different states.
# 
# A Migration PUMA may correspond to either one or multiple PUMAs, so the codes and 
# extents for PUMAs of residence 1 year ago (in MIGPUMA1) and the PUMA of current residence (in PUMA) occasionally differ.

## There's more than one puma in one migpuma
## one migpuma for every puma
reg_d <- reg_d %>% 
  left_join(puma_dic %>% 
              dplyr::select(state, puma, migpuma1) %>% 
              rename("newpuma"=migpuma1), by = c("statefip"="state", "puma")) 

reg_d %>% 
  filter(countyfip!=0) %>% 
  group_by(newpuma,statefip) %>% 
  summarise(ll = length(unique(countyfip))) %>% 
  filter(ll!=1)

reg_d %>% 
  filter(countyfip!=0) %>% 
  group_by(countyfip,statefip) %>% 
  summarise(ll = length(unique(newpuma))) %>% 
  filter(ll!=1)


obs_counties <- reg_d %>% filter(countyfip!=0) %>%group_by(statefip, newpuma,countyfip) %>% slice(1) %>%  dplyr::select(statefip, countyfip, newpuma)




#migrate1: migration status, 1 year [general version]
# Labels:
#   value                label
# 0                  n/a
# 1           same house
# 2   moved within state
# 3 moved between states
# 4  abroad one year ago
# 9              unknown

which(reg_d$migrate1==0)

#check (fine)
reg_d %>% 
  filter(migrate1 == 1 & (newpuma != migpuma1)) %>% 
  dplyr::select(migrate1, puma, newpuma, migpuma1)

reg_d %>% 
  filter(migrate1 == 2 & (newpuma != migpuma1)) %>% 
  dplyr::select(migrate1, puma, newpuma, migpuma1)

reg_d %>% 
  filter(migrate1 == 2 & (newpuma == migpuma1)) %>% 
  dplyr::select(migrate1, puma, newpuma, migpuma1)


reg_d <- reg_d %>% 
  mutate(samepuma = newpuma == migpuma1) %>% 
  mutate(migrant = NA, 
         migrant = ifelse(migrate1 == 1, 0, migrant),
         migrant = ifelse(migrate1 == 2 & !samepuma, 1, migrant), 
         migrant = ifelse(migrate1 == 2 & samepuma, 0, migrant), 
         migrant = ifelse(migrate1 %in% c(3,4), 1, migrant)) %>% 
  mutate(migrant = perwt*migrant) %>% 
  group_by(year, age, newpuma, statefip) %>% 
  summarise(migrant = sum(migrant), 
         totalpop = sum(perwt)) 

reg_d <- split(reg_d, reg_d$statefip)


reg_d <- lapply(reg_d, function(dd){
  expand.grid(year = sort(unique(dd$year)),
              newpuma = as.numeric(sort(unique(dd$newpuma))),
              age = as.numeric(sort(unique(dd$age)))) %>%
    left_join(dd, by = c("year","age","newpuma"))
})


reg_d$`6` %>% head(5)
reg_d$`36` %>% head(5)
reg_d$`48` %>% head(5)

lapply(reg_d, function(dd){dd %>% filter(is.na(migrant)|is.na(totalpop))})

reg_d$`36` %>% filter(age==27)

reg_d$`36` %>% mutate(dd = totalpop - migrant) %>% arrange(-desc(dd)) %>% filter(age<=75)
reg_d$`36` %>% filter(newpuma==1300 & age==27)


save(file="data_setup2.RData", list = c("work_d", "county_dic","puma_dic", "reg_d","obs_counties"))
