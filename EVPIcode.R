library(tidyverse)
library(patchwork)
######values for calculations
WRlow <- 1000
WRhigh <- 10000
FRlow <- 10000
FRhigh <- 100000
BH.a <- 1216
DD <- .00004803 # strong density dependence
DI <- DD * .4 #weak density dependence

######dewater estimates and water usage dataframe
dewater <- data.frame(Alts = c(1,2,3,4),
                      Fall = c(0.31,0,0.03,0.06),
                      Winter = c(0.06,0.27,0.09,0.08),
                      FlowTAF = c(4840.75, 4090.197, 4470.124, 4531.821))

######consequence table
consequence <- expand.grid(Alts = c(1,2,3,4),
                           hypothesis = factor(c('H1','H2','H3','H4'))) %>%
  left_join(dewater, by = 'Alts') %>%
  mutate(WRfry = case_when(hypothesis == 'H1'~ (WRhigh * (1-Winter) * BH.a)/(1+DI*WRhigh*(1-Winter)),
                           hypothesis == 'H2'~ (WRhigh * (1-Winter) * BH.a)/(1+DD*WRhigh*(1-Winter)),
                           hypothesis == 'H3'~ (WRlow * (1-Winter) * BH.a)/(1+DI*WRlow*(1-Winter)),
                           hypothesis == 'H4'~ (WRlow * (1-Winter) * BH.a)/(1+DD*WRlow*(1-Winter)))) %>% 
  mutate(FRfry = case_when(hypothesis == 'H1'~ (FRhigh * (1-Fall) * BH.a)/(1+DI*FRhigh*(1-Fall)),
                           hypothesis == 'H2'~ (FRhigh * (1-Fall) * BH.a)/(1+DD*FRhigh*(1-Fall)),
                           hypothesis == 'H3'~ (FRlow * (1-Fall) * BH.a)/(1+DI*FRlow*(1-Fall)),
                           hypothesis == 'H4'~ (FRlow * (1-Fall) * BH.a)/(1+DD*FRlow*(1-Fall)))) %>%
  mutate(NormalWR = (WRfry-min(WRfry))/(max(WRfry)-min(WRfry)),
         NormalFR = (FRfry-min(FRfry))/(max(FRfry)-min(FRfry)),
         Normalwater = (FlowTAF-max(FlowTAF))/(min(FlowTAF)-max(FlowTAF)))

######creates a table of objective weights
objectiveWeights <- expand.grid(Alts = c(1,2,3,4), hypothesis = factor(c('H1', 'H2', 'H3', 'H4')), 
                                Obj1 = seq(0, 1,.05),
                                Obj2 = seq(0, 1,.05),
                                Obj3 = seq(0, 1,.05)) %>%
  rowwise() %>%
  mutate(n = sum(c(Obj1, Obj2, Obj3))) %>% #sum all objective weights for filtering
  filter(n == 1) %>% #filter for only rows that have sum(objective weights) = 1
  filter_if(is.numeric, all_vars((.) != 0)) %>% #filter out any rows that have objective weights = 0
  select(-6)

######a table computing composite scores for each objective weight combo
complist = list() #creating a list for storing values
for(i in 1:nrow(objectiveWeights)){ #for loop for looping through objective weights table
  filter <- objectiveWeights[i,] #setting each row of objective weights column as a filter object
  temp <- consequence %>% 
    filter(hypothesis == filter$hypothesis & Alts == filter$Alts) %>% #filtering consequence table by the filter object
    select(8,9,10) #selecting desired columns from consequence table for a calculating composite scores
  filter2 <- filter %>% select(3,4,5) #selecting desired columns from objective weights table for calculating composite scores
  comp <- sum(unlist(filter2) * unlist(temp)) #calculating composite scores for each row (similar to SUMPRODUCT)
  complist[[i]] <- comp #storing that into a list
}

######final table with a composite score for each combination of objective weights
finaltemp <- bind_cols(objectiveWeights, do.call(rbind,complist)) %>%
  rename('Composite' = 6) #binding the complist with the objective weights table

######WR EVPI
WR <- finaltemp %>% 
  filter(Obj2 == Obj3) %>% 
  mutate(Objective = 'Winter-run') %>%
  mutate(Density = factor(if_else(hypothesis %in% c('H1', 'H3'), 'Weak DD', 'Strong DD'), 
                          levels = c('Weak DD', 'Strong DD')),
         State = if_else(hypothesis %in% c('H1', 'H2'), 'High abundance', 'Low abundance'))

WRcertain <- WR %>%
  group_by(Obj1, hypothesis) %>% 
  summarize(Certainty = max(Composite)) %>%
  ungroup() %>% 
  group_by(Obj1) %>% 
  summarize(Certainty = mean(Certainty))

WRuncertain <- WR %>%
  group_by(Obj1, Alts) %>% 
  summarize(uncertainty = mean(Composite)) %>%
  ungroup() %>% 
  group_by(Obj1) %>% 
  summarize(uncertainty = max(uncertainty))

wrEVPI <- WRcertain %>% left_join(WRuncertain, by = 'Obj1') %>%
  mutate(EVPI = abs(Certainty - uncertainty)) %>%
  mutate(EVPIperc = (EVPI/Certainty)*100) %>%
  mutate(Objective = 'Winter-run')

######FR EVPI
FR <- finaltemp %>% 
  filter(Obj1 == Obj3) %>% 
  mutate(Objective = 'Fall-run') %>%
  mutate(Density = factor(if_else(hypothesis %in% c('H1', 'H3'), 'Weak DD', 'Strong DD'), 
                          levels = c('Weak DD', 'Strong DD')),
         State = if_else(hypothesis %in% c('H1', 'H2'), 'High abundance', 'Low abundance'))

FRcertain <- FR %>%
  group_by(Obj2, hypothesis) %>% 
  summarize(Certainty = max(Composite)) %>%
  ungroup() %>% 
  group_by(Obj2) %>% 
  summarize(Certainty = mean(Certainty))

FRuncertain <- FR %>%
  group_by(Obj2, Alts) %>% 
  summarize(uncertainty = mean(Composite)) %>%
  ungroup() %>% 
  group_by(Obj2) %>% 
  summarize(uncertainty = max(uncertainty))

frEVPI <- FRcertain %>% left_join(FRuncertain, by = 'Obj2') %>%
  mutate(EVPI = abs(Certainty - uncertainty)) %>%
  mutate(EVPIperc = (EVPI/Certainty)*100) %>%
  mutate(Objective = 'Fall-run')

######Water-cost EVPI
Water <- finaltemp %>% 
  filter(Obj1 == Obj2) %>% 
  mutate(Objective = 'Water-cost') %>%
  mutate(Density = factor(if_else(hypothesis %in% c('H1', 'H3'), 'Weak DD', 'Strong DD'), 
                          levels = c('Weak DD', 'Strong DD')),
         State = if_else(hypothesis %in% c('H1', 'H2'), 'High abundance', 'Low abundance'))

watercertain <- Water %>%
  group_by(Obj3, hypothesis) %>% 
  summarize(Certainty = max(Composite)) %>%
  ungroup() %>% 
  group_by(Obj3) %>% 
  summarize(Certainty = mean(Certainty))

wateruncertain <- Water %>%
  group_by(Obj3, Alts) %>% 
  summarize(uncertainty = mean(Composite)) %>%
  ungroup() %>% 
  group_by(Obj3) %>% 
  summarize(uncertainty = max(uncertainty))

waterEVPI <- watercertain %>% left_join(wateruncertain, by = 'Obj3') %>%
  mutate(EVPI = abs(Certainty - uncertainty)) %>%
  mutate(EVPIperc = (EVPI/Certainty)*100) %>%
  mutate(Objective = 'Water Cost')

#######annotation for labeling graph
annotation <- data.frame(
  Obj1 = c(.20, .8, .20, .8),
  Composite = c(.1, .1, .9, .9),
  label = c("H1", "H2", "H3", "H4"),
  Density = factor(c('Weak DD', 'Strong DD', 'Weak DD', 'Strong DD'),
                   levels = c('Weak DD', 'Strong DD')),
  State = c('High abundance', 'High abundance', 'Low abundance', 'Low abundance')
)

#######graphing composite scores with objective weights
WRobj <- ggplot() + 
  geom_line(data = WR, mapping = aes(x = Obj1, y = Composite, color = factor(Alts))) +
  #geom_text(data = annotation, mapping = aes(x = Obj1, y = Composite, label = label)) +
  labs(x = '', y = 'Composite Score', title = 'Winter-run Objective', color = 'Alternatives') +
  facet_grid(State~Density) +
  scale_color_manual(values = c('1','2','3','4'), labels = c('Winter-run focused', 'Fall-run focused',
                                                             'Both runs focused', 'Baseline average'))
WRobj

FRobj <- FR %>%
  ggplot() + 
  geom_line(aes(x = Obj2, y = Composite, color = factor(Alts))) +
  labs(x = '', y = '', title = 'Fall-run Objective', color = 'Alternatives') +
  facet_grid(State~Density) +
  scale_color_manual(values = c('1','2','3','4'), labels = c('Winter-run focused', 'Fall-run focused',
                                                             'Both runs focused', 'Baseline average'))
FRobj

Waterobj <- Water %>%
  ggplot() + 
  geom_line(aes(x = Obj3, y = Composite, color = factor(Alts))) +
  labs(x = '', y = '', title = 'Water Cost Objective', color = 'Alternatives') +
  facet_grid(State~Density) +
  scale_color_manual(values = c('1','2','3','4'), labels = c('Winter-run focused', 'Fall-run focused',
                                                             'Both runs focused', 'Baseline average'))
Waterobj

#graphing EVI with objective weights
wrEVPIgraph <- wrEVPI %>% ggplot(aes(x = Obj1, y = EVPIperc)) +
  geom_col() +
  ylim(0,4) +
  labs(x = '', y = 'Percent EVPI')
wrEVPIgraph

frEVPIgraph <- frEVPI %>% ggplot(aes(x = Obj2, y = EVPIperc)) +
  geom_col() +
  ylim(0,4) +
  labs(x = 'Objective Weight', y = '')
frEVPIgraph

waterEVPIgraph <- waterEVPI %>% ggplot(aes(x = Obj3, y = EVPIperc)) +
  geom_col() +
  ylim(0,4) +
  labs(x = '', y = '')
waterEVPIgraph

finalGraphobj <- (WRobj|FRobj|Waterobj) & theme_bw() & theme(legend.position = "bottom", 
                                                             legend.text = element_text(size = 7)) &
  guides(color = guide_legend(nrow = 2))
finalGraphEVPI <- (wrEVPIgraph|frEVPIgraph|waterEVPIgraph) & theme_bw()


final2 <- finalGraphobj + plot_layout(guides = 'collect')
final2
finalGraph <- final2/finalGraphEVPI + plot_layout(heights = c(6,2))
finalGraph
ggsave(finalGraph, file = 'CompositeScore.png', height = 6, width = 11)