library(data.table)
library(tidyverse)
library(Matrix)
pairwiseType <- function(x){
  x <- t(x)
  mm <- as(x, "dgCMatrix")
  d <- tcrossprod(mm)
  denom <- matrix(rep(rowSums(x), ncol(d)), ncol = ncol(d), byrow = FALSE)
  denom <- denom + t(denom)
  return(as.matrix(2*d/denom))
}


# Load data
bs_data <- read_csv('Data/Bootstrapping_isolate_data.csv')
bs_data$age_simple_2 <- str_sub(bs_data$age_simple,1,3)
varmat <- as.tibble(fread('Data/For_Shai_1099_isolate_42399_DBLa_types_at_96%_seq_ID.csv'))
ncol(varmat)

# Randomly select one isolate from each age group and calculate their PTS
bootstraped_PTS <- function(season, ages, nsim=1000){
  # season=1;ages=c('Chi','Adu')
  bs_PTS <- c()
  for (n in 1:nsim){
    # print(n)
    x <- sample(subset(bs_data, Season==season & age_simple_2==ages[1])$SampleID,1)
    y <- sample(subset(bs_data, Season==season & age_simple_2==ages[2])$SampleID,1)
    if (x==y){next} # For within-group comparison. To avoid sampling the same isolate
    PTStable <- data.matrix(varmat[,c(x,y)])
    bs_PTS <- c(bs_PTS, pairwiseType(PTStable)[2,1])
  }
  bs_PTS <- as.tibble(bs_PTS)
  bs_PTS$comparison <- paste(ages,collapse='_')
  bs_PTS$Season <- season
  names(bs_PTS)[1] <- 'PTS'
  return(bs_PTS)
}
 

# comparison_list <- expand.grid(age1=c('Ado','Chi'),age2=c('Adu','Chi'),season=1:2, stringsAsFactors = F)
comparison_list <- expand.grid(age1=unique(bs_data$age_simple_2),age2=unique(bs_data$age_simple_2),season=1:2, stringsAsFactors = F)
comparison_list <- comparison_list[-c(4,7,8,13,16,17),]
# comparison_list <- comparison_list[!duplicated(comparison_list),]
# comparison_list <- subset(comparison_list, age1 != age2)

comparison_results <- NULL
for (l in 1:nrow(comparison_list)){
  age1 <- comparison_list[l,'age1']
  age2 <- comparison_list[l,'age2']
  season <- comparison_list[l,'season']
  print(paste(season,age1,age2,sep=' | '))
  do_bootstrap <- bootstraped_PTS(season = season, ages=c(age1,age2), 5000)
  comparison_results <- rbind(comparison_results, do_bootstrap)
}




comparison_results %>% 
  mutate(comparison=factor(comparison,levels=c('Chi_Chi','Ado_Ado','Adu_Adu','Ado_Chi','Adu_Chi','Adu_Ado'))) %>% 
  ggplot(aes(x=comparison,y=PTS))+
  geom_violin()+
  geom_boxplot(width=0.3)+
  scale_y_continuous(breaks=seq(0,0.3,0.05))+
  facet_wrap(~Season)
