functions {
  
   // Function to convert from real to integer
   int bin_search(real x, int min_val, int max_val){ 
    // This assumes that min_val >= 0 is the minimum integer in range, 
    //  max_val > min_val,
    // and that x has already been rounded. 
    //  It should find the integer equivalent to x.
    int range = (max_val - min_val+1)/2; // We add 1 to make sure that truncation doesn't exclude a number
    int mid_pt = min_val + range;
    int out;
    while(range > 0) {
      if(x == mid_pt){
        out = mid_pt;
        range = 0;
      } else {
        // figure out if range == 1
        range =  (range+1)/2; 
        mid_pt = x > mid_pt ? mid_pt + range: mid_pt - range; 
        }
    }
    return out;
  }
  
   array[] int sub(array[] int x, int y) {
        int x_size = size(x);
        array[x_size] int z;
        for (i in 1:x_size){
          z[i] = x[i] - y;
        }
        return z;
  }
  
   // Longitudinal model ODE
   real SimeoniModel(
     
     real t,    // time of the observation,
     real tp,   // time of previous observation
     real x,    // state vector
     real AUC,
     real AUCRes,
     array[] real theta

     ) {
    
    real L0; 
    real K2;
    real R;
    
    // PD params
    L0 = theta[1];
    K2 = theta[2];
    R = theta[4];
    
    // PD model
    real y;
    y = x*exp((L0-K2*exp(-R*AUCRes)*AUC)*(t-tp));
    
    return y;
        
    }


  // Survival 
  
    /*basic hazard function*/
  real h0(real mu,
          real sigma,
          real time
  ){
      
      real h;
      
      h = exp(lognormal_lpdf(time|mu,sigma))/(1-exp(lognormal_lcdf(time|mu,sigma)));
  
      return(h);
    }

  real hazard(
           real time,
           real mu,
           real sigma,
           real LIV,
           real LUN,
           real LN,
           real OTH,
           real Tstart
           ){
  
  real h_base;
  real h;
  
  real mu_cov = mu + LIV + LUN + LN + OTH;
  h_base=fmax(machine_precision(),h0(mu_cov,sigma,time)); // Here I have to add all the contribution to calculate the h0 -> Just the parameters that characterize the baseline hazard function
  
  if(time < Tstart) { h = 0;}
  else {h=h_base;}
  
  return(h);
  
  }
  
  
  real LogSurvival(array[] real t,
                   real time,
                   real mu,
                   real sigma,
                   vector LIV,
                   vector LUN,
                   vector LYN,
                   vector OTH,
                   real Tstart,
                   array[] real p
                ){  
                  

    real inth = 0;
    
    for(i in 1:size(t)) {
      inth = inth + p[i]*hazard(t[i],mu,sigma,LIV[i],LUN[i],LYN[i],OTH[i],Tstart);
    }
 
    return(-inth*time/2);
    
  }
  
  // Within chain parallization function
  
  real partial_sum(
    
                   array[] int subject,
                   int start_subject,
                   int end_subject,
                   array[] int start,
                   array[] int end,
                   array[] real time,
                   array[] real amt,
                   array[] real evid,
                   array[] int iObs,
                   real sigma1,
                   real sigma2,
                   array[] real cObs,
                   array[,,] real rho,
                   array[] int NLES,
                   array[,] int LLES,
                   array[,] int LESION,
                   real mu,
                   real sigma,
                   real betaLIVER,
                   real betaLUNG,
                   real betaLYN,
                   real betaOTHER,
                   array[] int lenObsOcc,
                   array[] real T,
                   array[] int C,
                   array[] real Tstart,
                   array[] int loq,
                   array[,,] int GK,
                   array[] real p
                   ) { 
      
    int index_offset =  start[start_subject] - 1;                               // Index offset
    
    int nEvents_sub = end[end_subject] - start[start_subject] + 1;              // Number of the events in the thread
    int nSubjects_sub = size(subject);                                              // Number of subject in the thread
    
    int nObs_sub = 0;                                                           
    for(i in 1:nSubjects_sub) {
      nObs_sub = nObs_sub + lenObsOcc[start_subject + i - 1];                   // Compute the number of observation in the thread 
    }
    
    // Here goes the longitudinal model //
    
    vector[nEvents_sub] tumor;

    for (j in 1:nSubjects_sub) {
          
      for(nl in 1:NLES[start_subject + j - 1]) {
          
        int idx_start;
        int start_offset;
        if(nl == 1){
          idx_start = start[start_subject + j - 1];
          start_offset =idx_start - index_offset;
        } else {
          idx_start = start[start_subject + j - 1] + sum(LLES[start_subject + j - 1,1:(nl-1)]);
          start_offset =idx_start - index_offset;
        }
        
          array[LLES[start_subject + j - 1,nl]] real t = time[idx_start:(idx_start + LLES[start_subject + j - 1,nl] - 1)];  
          array[LLES[start_subject + j - 1,nl]] real AUCs = amt[idx_start:(idx_start + LLES[start_subject + j - 1,nl] - 1)];
          array[LLES[start_subject + j - 1,nl]] real ev = evid[idx_start:(idx_start + LLES[start_subject + j - 1,nl] - 1)];
          
          real AUC = AUCs[1];
          real AUCres = AUCs[1];
          
          for (i in 1:LLES[start_subject + j - 1,nl]) {
          
            if(i == 1) {
                tumor[start_offset+i-1] = rho[start_subject + j - 1,nl,3];
            } else {
                tumor[start_offset+i-1] = SimeoniModel(t[i],t[i-1],tumor[start_offset+i-2],AUC,AUCres,to_array_1d(rho[start_subject + j - 1,nl,])) ;
            }
                      
          if(ev[i] == 4){
            if(AUCres == 0) {AUCres = AUCres + AUCs[i]*14;} // AUCs is Cmean over 14 days -> Cmean*14 is dailyAUC
            else {AUCres = AUCres - AUC*14*0.5 + AUCs[i]*14;}
            AUC = AUCs[i];
          }
        
         }
       }
    }


  real ll = 0;  // LogLikelihood - Initialized to zero
  
  int st;      // Starting index of the first observations
  
  if(start_subject == 1) { st = 1; }
  else { st = sum(lenObsOcc[1:(start_subject-1)])+1;  }
  
  // Now: Lets compute the likelihood at each observation - LONGITUDINAL PART -> OK ONLY IF LESIONS SHARE SAME ERROR MODEL
  
  for(s in 1:nObs_sub) {
      
      if(loq[st + s - 1] == 0) {ll = ll+normal_lpdf(cObs[st + s - 1] | tumor[iObs[st + s - 1] - iObs[st] + 1], sigma2+fmax(machine_precision(), tumor[iObs[st + s - 1] - iObs[st] + 1]*sigma1));} // Shifting the tumor with offset}
      if(loq[st + s - 1] == 1) {ll = ll+normal_lcdf(cObs[st + s - 1] | tumor[iObs[st + s - 1] - iObs[st] + 1], sigma2+fmax(machine_precision(), tumor[iObs[st + s - 1] - iObs[st] + 1]*sigma1));}

  }
  
  // 
  
  vector[15] livertumor;
  vector[15] lungtumor;
  vector[15] lyntumor;
  vector[15] othertumor;
  
  array[15] int indexesGK;
  
  for (j in 1:nSubjects_sub) {
    
      livertumor = rep_vector(0,15);
      lungtumor = rep_vector(0,15);
      lyntumor = rep_vector(0,15);
      othertumor = rep_vector(0,15);
  
    for(nl in 1:NLES[start_subject + j - 1]) {
     
       for (i in 1:15) {
          indexesGK[i] = GK[i,start_subject + j - 1,nl] - index_offset;
       }
     
        if(LESION[start_subject + j - 1,nl] == 1) {livertumor = livertumor + (tumor[indexesGK]);}
        if(LESION[start_subject + j - 1,nl] == 2) {lungtumor = lungtumor + (tumor[indexesGK]);}
        if(LESION[start_subject + j - 1,nl] == 3) {lyntumor = lyntumor + (tumor[indexesGK]);}
        if(LESION[start_subject + j - 1,nl] == 4) {othertumor = othertumor + (tumor[indexesGK]);}

    }
    
    if(sum(livertumor) != 0) {livertumor = (livertumor);}
    if(sum(lungtumor) != 0) {lungtumor = (lungtumor);}
    if(sum(lyntumor) != 0) {lyntumor = (lyntumor);}
    if(sum(othertumor) != 0) {othertumor = (othertumor);}

    if(C[start_subject + j - 1] == 1) {ll = ll + (log(hazard(T[start_subject + j - 1],mu,sigma,betaLIVER*livertumor[15],betaLUNG*lungtumor[15],betaLYN*lyntumor[15],betaOTHER*othertumor[15],Tstart[start_subject + j - 1]))+LogSurvival(time[GK[,start_subject + j - 1,1]],T[start_subject + j - 1],mu,sigma,betaLIVER*livertumor,betaLUNG*lungtumor,betaLYN*lyntumor,betaOTHER*othertumor,Tstart[start_subject + j - 1],p));}
    if(C[start_subject + j - 1] == 0) {ll = ll + LogSurvival(time[GK[,start_subject + j - 1,1]],T[start_subject + j - 1],mu,sigma,betaLIVER*livertumor,betaLUNG*lungtumor,betaLYN*lyntumor,betaOTHER*othertumor,Tstart[start_subject + j - 1],p);}
    
  }
  
  return ll;
  
}

}

data{
  
  // Longitudinal - Events
  
  int<lower = 1> nt;                       // No. of Total Points
  
  array[nt] real<lower = 0> time;          // Time vector               
  array[nt] real<lower = 0> amt;           // Amount
  array[nt] int<lower = 0> evid;           // Evid

  // Longitudinal - Observations
  
  int<lower = 1> nObs;                     // No. of Observations
  vector[nObs] y;                          // Observations
  array[nObs] int<lower = 1> iObs;         // Row indexes that contains observations

  int<lower = 1> nSubjects;                // No. of PDXs
  
  array[nSubjects] int<lower = 1> start;       // Starting indexes of each occasion 
  array[nSubjects] int<lower = 1> end;         // Ending indexes of each occasion 

  array[nSubjects] int <lower=1> lenObs;    // Number of Observation per Occasion
  
  // Survival //
  
  array[nSubjects] real<lower = 1> T;
  array[nSubjects] int<lower = 0> C;
  
  array[nSubjects] real Tstart; // Start of the treatment

  int<lower = 1> maxles;
  array[nSubjects] int<lower = 0> NLES;
  array[nSubjects,maxles] int<lower = 0> LLES;
  
  // Location matrixes //
  array[nSubjects,maxles] int<lower = 0> LESION; // Let's take patient i: LESION[i,] will be [4 0 0 0 0] if the patient lesion T1 is in the 4th organ

  array[nObs] int<lower = 0> loq;
  
  array[15,nSubjects,maxles] int GK;
  array[15] real p;

}

transformed data{
  
  // Calculate the length of each ID observation
  
  array[nSubjects] int len;
  for(i in 1:nSubjects){
    len[i] = end[i] - start[i] + 1;
  }
  
  // index for partial_sum function
  
  array[nSubjects] int subject;
  for (i in 1:nSubjects) subject[i] = i;
   
  int grain_size = 1;                   // No. of terms in the sum computed on each threas
  
}

parameters {
  
  // Population-level longitudinal model parameters
  real<lower = 0> L0_pop;
  real<lower = 0> K2_pop;
  real<lower = 0> TV0_pop;
  real<lower = 0> R_pop;

  // Individual-level longitudinal model parameters
  vector[nSubjects] eta_L0;
  vector[nSubjects] eta_K2;
  vector[nSubjects] eta_R;
  vector[nSubjects] eta_TV0;
  
  // Location effect longitudinal model parameters
  array[3] real xi_L0;
  array[3] real xi_K2;
  array[3] real xi_R;
  array[3] real xi_TV0;
  
  // Lesion - level longitudinal model parameters
  matrix[nSubjects,maxles] rho_L0;
  matrix[nSubjects,maxles] rho_K2;
  matrix[nSubjects,maxles] rho_R;
  matrix[nSubjects,maxles] rho_TV0;

  vector<lower=0> [4] omega;
  vector<lower=0> [4] omega2;
  real lambda;
  
  // Error model
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  
  // Survival 
  real<lower=0> mu;
  real<lower=0> sigma;
  
  // Link functions
  real betaLIVER;
  real betaLUNG;
  real betaLYN;
  real betaOTHER;
  
}

transformed parameters {
  
  // Distibution of individual parameters - mu referencing 
  
  vector<lower=0>[nSubjects] L0 = L0_pop * exp(eta_L0* omega[1]);
  vector<lower=0>[nSubjects] K2 = K2_pop * exp(eta_K2* omega[2]);
  vector<lower=0>[nSubjects] TV0 = TV0_pop * exp((exp(eta_TV0* omega[3])^lambda - 1)/lambda);
  vector<lower=0>[nSubjects] R = R_pop * exp(eta_R* omega[4]);
  
  // Matrix of parameters
  
  matrix<lower = 0>[nSubjects, 4] theta;    // nOccasion rows while the individual parameters are on patient level
  
  for(ns in 1:nSubjects) {
    theta[ns,1] = L0[ns];
    theta[ns,2] = K2[ns];
    theta[ns,3] = TV0[ns];
    theta[ns,4] = R[ns];
  }
  
  // Location effect longitudinal model parameters -> constain it to be zero
  array[4] real xi_tilde_L0;
  array[4] real xi_tilde_K2;
  array[4] real xi_tilde_R;
  array[4] real xi_tilde_TV0;
  
  xi_tilde_L0[1] = 0;
  xi_tilde_K2[1] = 0;
  xi_tilde_R[1] = 0;
  xi_tilde_TV0[1] = 0;
  
  xi_tilde_L0[2:4] = xi_L0;
  xi_tilde_K2[2:4] = xi_K2;
  xi_tilde_R[2:4] = xi_R;
  xi_tilde_TV0[2:4] = xi_TV0;
  
  array[nSubjects,maxles,4] real<lower = 0> rho;    // nOccasion rows while the individual parameters are on patient level

  for(nss in 1:nSubjects) {

    for(nls in 1:maxles) {

      if(nls <= NLES[nss]){
        rho[nss,nls,1] = theta[nss,1]*exp(xi_tilde_L0[LESION[nss,nls]] + rho_L0[nss,nls]* omega2[1]);
        rho[nss,nls,2] = theta[nss,2]*exp(xi_tilde_K2[LESION[nss,nls]] + rho_K2[nss,nls]* omega2[2]);
        rho[nss,nls,3] = theta[nss,3]*exp(xi_tilde_TV0[LESION[nss,nls]] + rho_TV0[nss,nls]* omega2[3]);
        rho[nss,nls,4] = theta[nss,4]*exp(xi_tilde_R[LESION[nss,nls]] + rho_R[nss,nls]* omega2[4]);
        
      } else {
        rho[nss,nls,1] = 0;
        rho[nss,nls,2] = 0;
        rho[nss,nls,3] = 0;
        rho[nss,nls,4] = 0;
      }
    }
  }
  
}

model{
  
  // Priors for Population-level longitudinal model parameters

  L0_pop ~ lognormal(log(0.0023),0.05);
  K2_pop ~ lognormal(log(1.1e-4),0.3);
  TV0_pop ~ lognormal(log(40),0.5);
  R_pop ~ lognormal(log(0.001),10);
  
  // Priors for Individual-level longitudinal model parameters
  
  eta_L0 ~ normal(0, 1);
  eta_K2 ~ normal(0, 1);
  eta_R ~ normal(0, 1);
  eta_TV0 ~ normal(0, 1);
  
  to_vector(rho_L0) ~ normal(0, 1);
  to_vector(rho_K2) ~ normal(0, 1);
  to_vector(rho_R) ~ normal(0, 1);
  to_vector(rho_TV0) ~ normal(0, 1);
  
  xi_L0 ~ normal(0, 1);
  xi_K2 ~ normal(0, 1);
  xi_R ~ normal(0, 1);
  xi_TV0 ~ normal(0, 1);
  
  array[4] real scale_omega = {1,
                              1,
                              1,
                              1};
                              
  omega ~ normal(scale_omega,[1,1,1,1]);
  omega2 ~ normal(scale_omega,[1,1,1,1]);
  
  lambda ~ normal(-0.3, 0.3);
  
  // Priors for Error model - Proportional

  sigma1 ~ lognormal(log(0.05), 0.5);
  sigma2 ~ lognormal(log(2), 0.5);
  
  // Survival //
  
  mu ~ lognormal(log(5), 1);
  sigma ~ lognormal(log(0.5), 0.5);
  
  betaLIVER ~ normal(0, 0.25);
  betaLUNG ~ normal(0, 0.25);
  betaLYN ~ normal(0, 0.25);
  betaOTHER ~ normal(0, 0.25);

  // Build the likelihood with the reduce sum function
  target += reduce_sum(partial_sum,
                       subject,
                       grain_size,
                       start, end, time,
                       amt,
                       evid,
                       iObs,
                       sigma1,
                       sigma2,
                       to_array_1d(y),
                       rho,
                       NLES,
                       LLES,
                       LESION,
                       mu,
                       sigma,
                       betaLIVER,
                       betaLUNG,
                       betaLYN,
                       betaOTHER,
                       lenObs,
                       T,
                       C,
                       Tstart,
                       loq,
                       GK,
                       p);
  
}