library(tidyverse)
library(zoo)
library(readxl)
library(xlsx)

#function for rounding flow data to the nearest corresponding flow in the Gard lookup table
round_to_nearest <- function(x, basis_vector) {
  nearest_value <- basis_vector[which.min(abs(basis_vector - x))]
  return(nearest_value)
}

#pulling in csv files
flow_data <- read_excel('InputData/alternatives.xlsx', sheet = 'Alternatives') #pull in flow scenario data
spawn <- read.csv('InputData/spawnTimingAll.csv')
lookup <- read.csv('InputData/model_dewater_flows.csv') %>% #pull in Gard lookup table
  rename('GardDewater' = 'Dewater', 'GardSpawn' = 'Spawn_Flows') %>%
  mutate(Run = str_to_lower(Run))

#Joining the flow and spawn data together
all <- flow_data %>% 
  gather(key = 'Alts', value = 'Discharge', 2:5) %>%
  left_join(select(spawn, c(4,6,7,13,3,5)), by = c('Day' = 'doy')) %>% #joining by Day
  rename('Spawn_Flows' = 'Discharge') %>%  #changing name for latter joining
  filter(!is.na(run)) %>% #getting rid of lines of dataframe with no spawning data
  group_by(Alts, year) %>% 
  mutate(Dewater = sapply(EmergDate, function(edate) { #adds minimum flow between spawn and emergence date
    if (!is.na(edate)) {
      min(Spawn_Flows[byDay <= edate & byDay >= byDay], na.rm = TRUE)
    } else {
      NA
    }
  })) %>% 
  ungroup()

all2 <- all %>%
  mutate(Boards = if_else(Day < 289, 'in', 'out')) %>% #adding ACID gates in or out;pulled out approx Oct. 15
  mutate(GardSpawn = sapply(all$Spawn_Flows, round_to_nearest, lookup$GardSpawn)) %>% #rounds spawn flows to nearest corresponding spawn flow in gard lookup
  mutate(GardDewater = sapply(all$Dewater, round_to_nearest, lookup$GardDewater)) %>% #rounds dewater flows to nearest corresponding dewater flow in gard lookup
  left_join(lookup, by = c('GardDewater', 'GardSpawn', 'Boards', c('run' = 'Run'))) %>% #join lookup table
  mutate(Prop_dewater = if_else(is.na(Prop_dewater), 0, Prop_dewater)) %>% #converts NA values to 0
  mutate(TotalDewater = prop * Prop_dewater) #calculates total dewatered by multiplyling proportion spawning with proportion expected to dewater

#summarizes total water volume by alternative
flow <- flow_data %>% gather(key = Alts, value = Flow, 2:5) %>%
  group_by(Alts) %>% summarize(FlowCFS = sum(Flow)) %>%
  mutate(FlowTAF = (FlowCFS*1.98)/1000)

#summarizes total dewatering by alternative and run
summary <- all2 %>%
  group_by(year, run, Alts) %>% summarize(TotalDewater = round(sum(TotalDewater),2)) %>%
  spread(key = run, value = TotalDewater) %>% 
  left_join(flow, by = 'Alts') #adds water volume to table

write.csv(summary, 'InputData/dewaterSummary.csv', row.names = FALSE) #save summary table of dewater consequences
