##### Load Libraries #####

library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(plyr)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(reshape2)
library(huxtable)
library(tidyr)
source("../functions/funs.R")


options(mc.cores = parallel::detectCores())
set_cmdstan_path("~/StanModels/Torsten/cmdstan/")

##### The previous dataset, after filtering for only administration lines, become the "Treatment dataset" #####

xdata.AUC <- read.table(file = "datafiles/dataAUC.csv",header=T) |> arrange(id,time)
xdata.AUC[which(xdata.AUC$id == "A6181170   88" & xdata.AUC$evid == 0),]$time <- xdata.AUC[which(xdata.AUC$id == "A6181170   88" & xdata.AUC$evid == 0),]$time - 16
xdata.AUC[which(xdata.AUC$id == "A6181170  119" & xdata.AUC$evid == 0),]$time <- xdata.AUC[which(xdata.AUC$id == "A6181170  119" & xdata.AUC$evid == 0),]$time - 24

xdata.AUC <- xdata.AUC |> filter(id != "A6181170  414")

xdata.AUC <- xdata.AUC |>
  filter(id %in% c(unique(xdata.AUC$id))) |>
  filter(evid != 0)

##### Here I have to add the individual lesion level ######

xdata <- read.table(file = "datafiles/DataFileLesion.csv",header=T,sep=",") |> arrange(id,TMMLNM,time)

#### Map T -> T1:5

xdata <- xdata |> group_by(id) |>  mutate(G = dense_rank(interaction(id, TMMLNM))) |>  ungroup() |> mutate(TMMLNM = paste0("T",as.character(G))) |> select(!G)  |> filter(TMMLNM != "T6")

# Fix Error patient "A6181170  137"
xdata[which(xdata$id ==  "A6181170  137" & xdata$TMMLNM == "T3"),]$TMMLNM = "T2"

xdata <- xdata |> filter(id != "A6181170  414")

xdata <- xdata |>
  filter(id %in% c(unique(xdata$id))) |> drop_na()

NLES <- (xdata |> 
           group_by(id) |>
           reframe(nles = length(unique(TMMLNM))))$nles


#### build the full dataset #####

xdata <- rbind(
  xdata,
  xdata.AUC |> left_join(xdata |> select(id,TMMDIS,TMMLNM) |> distinct(id,TMMDIS,TMMLNM),by="id",relationship = "many-to-many") |> mutate(loqcens = 0)
) |> arrange(id,TMMLNM,time)


###### Survival #####

xdata.surv <- read.table(file = "datafiles/datafileSurv.csv",header=T)
xdata.surv <- xdata.surv |> filter(id != "A6181170  414")
xdata.surv[(which(xdata.surv$id == "A6181170  344")),]$time <- 1


#### GK ####

x<-c(0.991455371120813, -0.991455371120813, 0.949107912342759, -0.949107912342759, 0.864864423359769, -0.864864423359769,-0.741531185599394	, 0.741531185599394	, 
     0.586087235467691, -0.586087235467691, 0.405845151377397, -0.405845151377397, -0.207784955007898, 0.207784955007898, 0 )

p<-c(0.022935322010529,0.022935322010529, 0.063092092629979,0.063092092629979, 0.104790010322250, 0.104790010322250,
     0.140653259715525,0.140653259715525, 	0.169004726639267, 	0.169004726639267, 0.190350578064785,0.190350578064785, 
     0.204432940075298, 0.204432940075298, 0.209482141084728
)


df.GK <- data.frame()
for(i in seq(1,nrow(xdata.surv))) {
  
  for(j in seq(1,length(unique((xdata |> filter(id == xdata.surv$id[i]))$TMMLNM)))) {
    
    TMMDIS <- (xdata |> filter(id == xdata.surv$id[i]) |> distinct(TMMDIS,TMMLNM))$TMMDIS[j]
    TMMLNM <- (xdata |> filter(id == xdata.surv$id[i]) |> distinct(TMMDIS,TMMLNM))$TMMLNM[j]
    
    df.GK <- rbind(df.GK,data.frame(time = xdata.surv$time[i]*(x+1)/2, y = 0, evid = 1, amount = 0, addl = 0, ii = 0, rate = 0, cmt = 1, ss = 0, TMMDIS = TMMDIS, TMMLNM = TMMLNM,id = xdata.surv$id[i], loqcens = 0))
    
    
  }
  
}

df.GK <- df.GK |> arrange(id,time)

xdata <- rbind(xdata,df.GK) |> arrange(id,TMMLNM,time)

##### Data - Treatment start #####

xdata.tr <- read.table(file = "TreatmentStart.csv",header=T,sep=",")
xdata.tr <- xdata.tr |> filter(PID_A != "A6181170  414")

##### Data - Longitudinal  #####

nt <- nrow(xdata)                           # Save all observations length
start <- (1:nt)[!duplicated(xdata$id)]      # Save vectors from which the ID start ...
end <- c(start[-1] - 1, nt)                 # ... and where it ends
nSubjects <- length(unique(xdata$id))       # Number of Subjects

##### Compile data list that has to be parsed in Stan ##### 

nobs <- length(which(xdata$evid == 0))
obsid <- which(xdata$evid == 0)

obsperid <- xdata |>
  dplyr::mutate(n = row_number()) |>
  filter(evid == 0) |>
  group_by(id) |>
  reframe(sz = n())

#### Data List #####

TV0_reg <- xdata |> group_by(id) |> reframe(y0 = y[which(time == min(time))])

LLES <- xdata |> group_by(id,TMMLNM) |> reframe(count = n()) |> pivot_wider(names_from = TMMLNM, values_from = count)
LLES[is.na(LLES)] <- 0

LLES.matrix <- as.matrix(LLES[,2:6])

TYPELES <- xdata |> distinct(id,TMMLNM,.keep_all = T) |> pivot_wider(names_from = TMMLNM, values_from = TMMDIS) # |> select(id,T1,T2,T3,T4,T5)
TYPELES[is.na(TYPELES)] <- "NA"

TYPELES <- TYPELES |> mutate(T1 = case_when(
  T1 == "LIVER" ~ 1,
  T1 == "LUNG" ~ 2,
  T1 == "LYMPH NODE" ~ 3,
  T1 == "NA" ~ NA,
  .default = 4
),
T2 = case_when(
  T2 == "LIVER" ~ 1,
  T2 == "LUNG" ~ 2,
  T2 == "LYMPH NODE" ~ 3,
  T2 == "NA" ~ NA,
  .default = 4
),
T3 = case_when(
  T3 == "LIVER" ~ 1,
  T3 == "LUNG" ~ 2,
  T3 == "LYMPH NODE" ~ 3,
  T3 == "NA" ~ NA,
  .default = 4
),
T4 = case_when(
  T4 == "LIVER" ~ 1,
  T4 == "LUNG" ~ 2,
  T4 == "LYMPH NODE" ~ 3,
  T4 == "NA" ~ NA,
  .default = 4
),
T5 = case_when(
  T5 == "LIVER" ~ 1,
  T5 == "LUNG" ~ 2,
  T5 == "LYMPH NODE" ~ 3,
  T5 == "NA" ~ NA,
  .default = 4
))

TYPELES[is.na(TYPELES)] <- 0

TYPELES <- TYPELES |> group_by(id) |> reframe(T1 = max(T1),
                                              T2 = max(T2),
                                              T3 = max(T3),
                                              T4 = max(T4),
                                              T5 = max(T5))

TYPELES.matrix <- as.matrix(TYPELES[,2:6])


### Organ indicator ####

xdata.plot <- xdata |>
  mutate(TMMDIS = case_when(
    TMMDIS == "LIVER" ~ "LIVER",
    TMMDIS == "LUNG" ~ "LUNG",
    TMMDIS == "LYMPH NODE" ~ "LYMPH NODE",
    .default = "OTHER"
  ))

locator <- xdata.plot |> group_by(id,TMMDIS) |> reframe() |> pivot_wider(names_from = TMMDIS, values_from = TMMDIS) |> select(id,LIVER,LUNG,`LYMPH NODE`,OTHER)
loc <- matrix(as.numeric(!is.na(locator)),nrow=nSubjects,ncol=5)[,2:5] |> unname()

#### Extract GK indexes ####

GK <- xdata |> mutate(rownum = row_number()) |> filter(evid == 1) |> select(id,TMMLNM,rownum) |> pivot_wider(names_from = TMMLNM, values_from = rownum, values_fn = list)
GK$T2[sapply(GK$T2, is.null)] <- list(rep(0,15))
GK$T3[sapply(GK$T3, is.null)] <- list(rep(0,15))
GK$T4[sapply(GK$T4, is.null)] <- list(rep(0,15))
GK$T5[sapply(GK$T5, is.null)] <- list(rep(0,15))

GKmat <- array(unlist(GK[, c("T1", "T2","T3","T4","T5")], recursive = FALSE), dim = c(nSubjects, 5, 15))

#### 

data <- list(nt = nt, # No. of Total Points
             
             time = as.numeric(xdata$time),
             amt = as.numeric(xdata$amount)/14,
             evid = as.numeric(xdata$evid),
             
             nObs = nrow(xdata |> filter(evid == 0)),
             y = (xdata |> filter(evid == 0))$y,
             iObs = (xdata |> dplyr::mutate(rn = row_number()) |> filter(evid == 0))$rn,
             
             nSubjects = nSubjects,
             
             start = start,
             end = end,
             
             lenObs = obsperid$sz,

             T = xdata.surv$time,
             C = xdata.surv$status,
             Tstart = xdata.tr$start,
             
             maxles = 5,
             NLES = NLES,
             LLES = LLES.matrix,
             LESION = TYPELES.matrix,
             
             loq = (xdata |> filter(evid == 0))$loqcens,
             loc = loc,
             
             GK = array(unlist(GKmat, recursive = FALSE), dim = c(15, nSubjects, 5)),
             p = p
             
)


# Initial conditions

init <- function(){
  list(L0_pop = 0.0023, 
       K2_pop = 1.1e-4,
       TV0_pop = 40,
       R_pop = 0.00001,
       
       eta_L0 = rep(0, nSubjects),
       eta_K2 = rep(0, nSubjects),
       eta_R = rep(0, nSubjects),
       eta_TV0 = rep(0, nSubjects),
       
       xi_L0 = rep(0,3),
       xi_K2 = rep(0,3),
       xi_R = rep(0,3),
       xi_TV0 = rep(0,3),
       
       rho_L0 = matrix(rep(0,nSubjects*5),nSubjects,5),
       rho_K2 = matrix(rep(0,nSubjects*5),nSubjects,5),
       rho_R = matrix(rep(0,nSubjects*5),nSubjects,5),
       rho_TV0 = matrix(rep(0,nSubjects*5),nSubjects,5),
       
       omega = c(0.8,1.5,0.5,2.5),
       omega2 = c(0.5,0.5,0.5,0.5),
       
       lambda = -0.3,
       
       sigma1 = 0.05,
       sigma2 = 5,
       
       mu = 5.5,
       sigma = 1,
       
       betaLIVER = 0,
       betaLUNG = 0,
       betaLYN = 0,
       betaOTHER = 0
  )
}

# Go with the fit

file <- file.path("Models/ModelMain_ContinuousTime.stan")
mod <- cmdstan_model(file,compile = F)
mod$check_syntax() #Check syntax
mod <- cmdstan_model(file,compile = T,cpp_options = list(stan_threads = TRUE))

fit <- mod$sample(data = data,
                  init = init,
                  chains = 4,
                  iter_warmup = 1000,
                  iter_sampling = 1000,
                  parallel_chains = 4,
                  refresh = 1,
                  seed = 123,
                  save_warmup = T,
                  max_treedepth = 10,
                  threads_per_chain = 2,
                  step_size = 0.1,
                  adapt_delta = .95)

fit$time()

  