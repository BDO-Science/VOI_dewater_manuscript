library(tidyverse)
library(CDECRetrieve)
library(dplyr)
library(zoo)

######function for calculating brood year day starting on April 1
byDay <- function(x, start.month = 4L) {
  start.yr <- as.integer(format(x, "%Y"))
  start.yr <- ifelse(as.integer(format(x, "%m")) < start.month, start.yr - 1, start.yr)
  start.date <- as.Date(paste(start.yr, start.month, "01", sep = "-"))
  as.integer(x - start.date + 1L)
}

######read in Fall run carcass data and clean up
FRdata <- read_csv('InputData/fallCarcass.csv') %>%
  mutate(Date = dmy(Date)) %>%
  filter(!is.na(Date)) %>%
  select(1, count = 2) %>%
  mutate(count = as.numeric(count)) %>%
  replace(is.na(.), 0) %>%
  mutate(run = 'fall')

######read in winter run carcass data and clean up
WRdata <- read.csv('InputData/winterCarcass.csv', check.names = FALSE) %>% 
  gather(key = 'Year', value = Count, 2:22) %>% 
  mutate(Year = as.integer(Year), Date = as.Date(Day, origin = "1970-01-01")) %>%
  mutate(Month = month(Date)) %>% 
  mutate(Year2 = year(Date)) %>% 
  mutate(Diff = Year - Year2) %>%
  mutate(Date = Date %m+% years(Diff)) %>%
  mutate(run = 'winter') %>%
  select(Date, count = Count, run) %>%
  replace(is.na(.), 0)

######combine winter run and fall run data
alldata <- bind_rows(FRdata, WRdata) %>%
  mutate(Date = Date - 7) %>% #subtracting 7 days based on guidance that carcass likely spawned 7 days prior to being sampled
  mutate(byDay = byDay(Date)) %>% #adding brood year day starting April 1
  mutate(year = factor(year(Date))) %>% #adding calendar year
  mutate(doy = yday(Date)) %>% #adding calendar day
  mutate(doy = if_else(doy > 100, doy, doy+365)) %>% #making calendar day extend past Jan 1
  arrange(byDay) %>% #arranging by byDay
  group_by(year,run) %>% #grouping by year and run
  mutate(prop = prop.table(count)) %>% #calculating proportion of carcasses by year and by run
  mutate(cumul_prop = cumsum(prop)) %>% #calculating cumulative proportion
  ungroup() %>%
  filter(year != 2024) #filtering out 2024

spawnGraph <- ggplot(alldata, aes(x = byDay, y = cumul_prop, color = run)) + #graph to visualize cumulative spawning data
  geom_line() +
  facet_wrap(~year)
spawnGraph

######fitting a logistic regression line to spawning data
split_data <- split(alldata, alldata$run) #Split data by 'group' for logisitic regression model

models <- lapply(split_data, function(subset) { #model logistic regression to each group
  glm(cumul_prop ~ byDay + year, 
      data = subset, 
      family = binomial(link = 'logit'))
})

fallmodel <- summary(models$fall) #summarize fall spawning model
wintermodel <- summary(models$winter) #summarize winter spawning model

fitted <- alldata %>% #use model to fit line to dataset
  mutate(pred = case_when(
    run == 'winter' ~ predict(models$winter, newdata = ., type = 'response'),
    run == 'fall' ~ predict(models$fall, newdata = ., type = 'response')
  ))

allyears <- fitted %>% #summarizing all years for graphing
  mutate(group = paste0(year,run)) %>%
  select(byDay, pred, group)

fittedGraph <- ggplot() + #graph for visualizing fitted data
  geom_line(allyears, mapping = aes(x = byDay, y = pred, group = group), color = 'grey') +
  geom_line(fitted, mapping = aes(x = byDay, y = pred, color = run), linewidth = 1) +
  facet_wrap(~year) +
  labs(x = 'Date', y = 'Cumulative Proportion') +
  theme_bw()
fittedGraph

######add emergent days based on water temperature
waterTemp <- cdec_query('KWK',
                        '25',
                        'H',
                        min(alldata$Date),
                        max(alldata$Date)) %>%
  mutate(Date = as.Date(datetime)) %>%
  mutate(parameter_value = if_else(parameter_value <= 32, NA_real_, parameter_value)) %>%
  arrange(datetime) %>%
  mutate(parameter_value = zoo::na.approx(parameter_value)) %>%
  group_by(Date) %>%
  summarize(tempF = mean(parameter_value, na.rm = TRUE)) %>% #calculating mean temperature by day
  ungroup() %>%
  mutate(tempC = (tempF-32)/1.8) #convert to Celsius

final <- fitted %>%
  left_join(waterTemp, by = c('Date')) %>%
  mutate(incubDays = exp((10.404-2.043*log(tempC + 7.575)))) %>% #calculating incubation days
  mutate(EmergDate = round(byDay + incubDays,0))

write.csv(final, file = 'InputData/spawnTimingAll.csv', row.names = FALSE)
