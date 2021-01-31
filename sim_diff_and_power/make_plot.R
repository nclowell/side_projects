###########################################################################
# Make plot from sim_diff_and_power.py output

###########################################################################

library(tidyverse)

###########################################################################

setwd("C:/Users/Natalie Lowell/SHARED_FOLDER/hyak_files/nl_minmig_pc04")
alpha = 0.05

###########################################################################
### functions

test_sig <- function(pval, alpha){
  if(pval < alpha){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

get_thresh <- function(FST, low_thresh, high_thresh){if(FST <= low_thresh){return(0)}
  if(FST > low_thresh & FST <= high_thresh){return(1)}
  if(FST > high_thresh){return(2)}}

###########################################################################
### threshold parameters for plot

low_thresh <- 0.0068 # empirical global FST
high_thresh <- 0.068 # 10X global FST

###########################################################################


res <- read_table2("minmig_pc04_all_results.txt") %>%
  select(-PropTestsDropped) %>%
  drop_na() %>%
  rowwise() %>%
  mutate(test_sig = test_sig(pval=Chi2Pval, alpha=alpha)) %>%
  mutate_at(vars(Popsize, mig_rate), funs(factor))

res$Popsize <- factor(res$Popsize, levels = c(500, 2500, 10000))


# short-term 10 generations drift
meanFST_gen10 <- res %>%
  filter(gens_drift == 10) %>%
  group_by(Popsize, mig_rate) %>%
  summarise(meanFST = mean(Fst)) %>%
  rowwise() %>%
  mutate(thresh = get_thresh(meanFST, low_thresh, high_thresh))
numTests_gen10 <- res %>%
  filter(gens_drift == 10) %>%
  group_by(Popsize, mig_rate) %>%
  count() %>%
  rename(NTests = n)
numTestsSig_gen10 <- res %>%
  filter(gens_drift == 10) %>%
  group_by(Popsize, mig_rate) %>%
  count(test_sig == TRUE) %>%
  rename(NTestsSig = n)
prop_gen10 <- full_join(numTests_gen10, 
                       numTestsSig_gen10, 
                       by = c("Popsize", "mig_rate")) %>%
  rowwise() %>%
  mutate(propSig = NTestsSig / NTests)
gen10 <- full_join(meanFST_gen10, prop_gen10,
                   by = c("Popsize", "mig_rate")) %>%
  select(Popsize, mig_rate, meanFST, propSig, thresh)

# long-term 100 generations drift
meanFST_gen100 <- res %>%
  filter(gens_drift == 100) %>%
  group_by(Popsize, mig_rate) %>%
  summarise(meanFST = mean(Fst)) %>%
  rowwise() %>%
  mutate(thresh = get_thresh(meanFST, low_thresh, high_thresh))
numTests_gen100 <- res %>%
  filter(gens_drift == 100) %>%
  group_by(Popsize, mig_rate) %>%
  count() %>%
  rename(NTests = n)
numTestsSig_gen100 <- res %>%
  filter(gens_drift == 100) %>%
  group_by(Popsize, mig_rate) %>%
  count(test_sig == TRUE) %>%
  rename(NTestsSig = n)
prop_gen100 <- full_join(numTests_gen100, 
                         numTestsSig_gen100, 
                         by = c("Popsize", "mig_rate")) %>%
  rowwise() %>%
  mutate(propSig = NTestsSig / NTests)
gen100 <- full_join(meanFST_gen100, prop_gen100,
                    by = c("Popsize", "mig_rate")) %>%
  select(Popsize, mig_rate, meanFST, propSig, thresh)

# equilibrium 1k generations drift
meanFST_gen1k <- res %>%
  filter(gens_drift == 1000) %>%
  group_by(Popsize, mig_rate) %>%
  summarise(meanFST = mean(Fst)) %>%
  rowwise() %>%
  mutate(thresh = get_thresh(meanFST, low_thresh, high_thresh))
numTests_gen1k <- res %>%
  filter(gens_drift == 1000) %>%
  group_by(Popsize, mig_rate) %>%
  count() %>%
  rename(NTests = n)
numTestsSig_gen1k <- res %>%
  filter(gens_drift == 1000) %>%
  group_by(Popsize, mig_rate) %>%
  count(test_sig == TRUE) %>%
  rename(NTestsSig = n)
prop_gen1k <- full_join(numTests_gen1k, 
                        numTestsSig_gen1k, 
                        by = c("Popsize", "mig_rate")) %>%
  rowwise() %>%
  mutate(propSig = NTestsSig / NTests)
gen1k <- full_join(meanFST_gen1k, prop_gen1k,
                   by = c("Popsize", "mig_rate")) %>%
  select(Popsize, mig_rate, meanFST, propSig, thresh)


### ------ plotting data

y_axis_labels <- c("1.95e-04", "3.91e-04", "7.81e-04", "1.56e-03",
                   "3.12e-03", "6.25e-03", "1.25e-02", "2.5e-02",
                   "5e-02", "1e-01")
                   

gen10_p <- ggplot(gen10, aes(x=Popsize, y=mig_rate)) +
  geom_tile(aes(fill=thresh)) +
  labs(x="Population size", y="Migration rate") +
  theme(legend.position = "none") +
  theme_classic() +
  scale_y_discrete(labels=y_axis_labels)
gen10_p

gen100_p <- ggplot(gen100, aes(x=Popsize, y=mig_rate)) +
  geom_tile(aes(fill=thresh)) +
  labs(x="Population size", y="") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none") 
gen100_p

gen1k_p <- ggplot(gen1k, aes(x=Popsize, y=mig_rate)) +
  geom_tile(aes(fill=thresh)) +
  labs(x="Population size", y="Migration rate") +
  theme(legend.position = "none") +
  theme_classic() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(legend.position = "none") 
gen1k_p




