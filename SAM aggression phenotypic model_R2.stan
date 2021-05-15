data {
  
  //indices and scalars used for model specification
  int<lower=1> N_sex; //total AG observations per sex
  int<lower=0> I; //total individuals
  int<lower=0> Im; //number of males
  int<lower=0> If; //number of females
  int<lower=1> Idyad; //number of dyads
  int<lower=1> idm[N_sex]; //index of male AG observations
  int<lower=1> idf[N_sex]; //index of female AG observations
  int<lower=1> idmw[Idyad]; //index of male FS observations
  int<lower=1> idfw[Idyad]; //index of female FS observations
  int<lower=1> dyadw[Idyad]; //index of dyads for FS
  int<lower=1> partners_m [Im,5]; //index of male partner IDs, first column is focal ID
  int<lower=1> partners_f [If,5]; //index of female partner IDs, first column is focal ID
  
  //empirical data
  real AG_m[N_sex]; //male aggression measurements
  real AG_f[N_sex]; //female aggression measurements
  real w[Idyad]; //FS dyad response
  real time[N_sex];
}

parameters {
  //population effects
  real alpha_0; //aggression global intercept
  real nu_0; //fitness global intercept
  real psi_1; //expected interaction coefficient
  real beta_B1; //average partner SRN intercept
  real beta_B2; //average partner SRN slopes
  real beta_N1; //selection gradients
  real beta_N2;
  real beta_S1;
  real beta_S2;
  real beta_D1;
  real beta_D2;
  real<lower=0, upper = 1> phi; //facilitates identification
  
  
  //random effects (standard deviations)
  vector<lower=0, upper = 1>[2] sd_P; //phenotypic SRN mu & psi SDs
  vector<lower=0, upper = 1>[2] sd_R; //male & female residual SDs
  real<lower=0, upper = 1>sd_delta; //residual of fitness

  cholesky_factor_corr[2] LP; //phenotypic SRN correlationsvector[N_sex]
  cholesky_factor_corr[2] LR; //sex-specific residual correlations
  
  matrix[I,2] std_devP; //individual-level SRN deviations
}

transformed parameters {
  matrix[I,2] SRN_P; //individual phenotypic SRN parameter values
  matrix[If, 2] partner_meanm; //average SRN parameters of males' partners
  matrix[Im, 2] partner_meanf; //average SRN parameters of females' partners
  
  //non-centered parameterization of random effects
  SRN_P =  std_devP * diag_pre_multiply(sd_P, LP)' ;

  //calculate the mean SRN parameters of each male's lifetime partners
  for(i in 1:Im) partner_meanm[i] = [mean(col(SRN_P[partners_m[i,2:5]],1)),
                                mean(col(SRN_P[partners_m[i,2:5]],2))];
  
  //calculate the mean SRN parameters of each female's lifetime partners
  for(i in 1:If) partner_meanf[i] = [mean(col(SRN_P[partners_f[i,2:5]],1)),
                                mean(col(SRN_P[partners_f[i,2:5]],2))];
}

model{
  
  //separate male and female vectors for efficiency
  matrix[Im,2] SRN_Pm = SRN_P[1:Im]; //male SRN phenotypic deviations
  matrix[If,2] SRN_Pf = SRN_P[(Im+1):I]; //female SRN phenotypic deviations
  
  //separate SRN intercepts and slopes (phenotypic deviations)
  vector[Im] mu_m = col(SRN_Pm,1); //SRN intercepts
  vector[If] mu_f = col(SRN_Pf,1); 
  vector[Im] psi_m = col(SRN_Pm,2); //SRN slopes
  vector[If] psi_f = col(SRN_Pf,2); 
  
  //separate mean partner SRN intercepts and slopes (deviations)
  vector[Im] mu_meanm = col(partner_meanm,1); //mean partner SRN intercept for males
  vector[If] mu_meanf = col(partner_meanf,1); //...for females
  vector[Im] psi_meanm = col(partner_meanm,2); //mean partner SRN slope for males
  vector[If] psi_meanf = col(partner_meanf,2); //...for females 
  
  //initialize vectors for constructing individual-centered linear predictors
  vector[N_sex] eta_Wm; //within-individual centered male SRN trait value 
  vector[N_sex] eta_Wf; //within-individual centered female SRN trait value
  
  vector[N_sex] meaneta_m; //individual male SRN trait value toward average partner 
  vector[N_sex] meaneta_f; //individual female SRN trait toward average partner
  vector[N_sex] eta_meanm; //average SRN partner values for males
  vector[N_sex] eta_meanf; //average SRN partner values for females
  
  vector[N_sex] linpred_m; //expected value for male responses
  vector[N_sex] linpred_f; //expected value for female responses
  vector[N_sex] zeta_m; //residuals for male responses
  vector[N_sex] zeta_f; //residuals for male responses
  vector[Idyad] w_pred; //linear predictor of fitness
  
  //Male and female aggression response model
  for (n in 1:N_sex) {
    
    //SRN trait values
    //assumes that n = 1 in the context of an ongoing social interaction
    //if n = 1 prior to social context, then specify eta[t=1] = mu_j instead
    if (time[n]==1)
      {
        //within-individual centered eta
        //male eta[t=1] = mu_j + psi_j*(mu_k - mu_meanK)
        eta_Wm[n] = mu_m[idm[n]] +  (psi_1 + psi_m[idm[n]])*(mu_f[idf[n]] - mu_meanm[idm[n]]) ;
        
        //female eta[t=1] = mu_k + psi_k*(mu_j - mu_meanJ)
        eta_Wf[n] = mu_f[idf[n]] +   (psi_1 + psi_f[idf[n]])*(mu_m[idm[n]] - mu_meanf[idf[n]]);
        
        //average individual eta
        //male eta[t=1] = mu_j + psi_j*mu_k
        meaneta_m[n] = mu_m[idm[n]] +    (psi_1 + psi_m[idm[n]])*mu_meanm[idm[n]];
        
        //female eta[t=1] = mu_k + psi_k*mu_j
        meaneta_f[n] = mu_f[idf[n]] +   (psi_1 + psi_f[idf[n]])*mu_meanf[idf[n]];
        
        //average partner eta[t=1]
        //average eta males' partners [t=1] = mu_meanK + psi_meanK*mu_j
        eta_meanm[n] = mu_meanm[idm[n]] +   (psi_1 + psi_meanm[idm[n]])*mu_m[idm[n]];
        
        //average eta females' partners [t=1] = mu_meanJ + psi_meanJ*mu_k
        eta_meanf[n] = mu_meanf[idf[n]] +   (psi_1 + psi_meanf[idf[n]])*mu_f[idf[n]];
      }    
    else
      {
        //within-individual centered eta
        //male eta[t=2] = mu_j + psi_j*(eta_k[t=1] - eta_meanK[t=1])
        eta_Wm[n] = mu_m[idm[n]] +   (psi_1 + psi_m[idm[n]])*(eta_Wf[n-1] - eta_meanm[n-1]);
        
        //female eta[t=2] = mu_k + psi_k*(eta_j[t=1] - eta_meanJ[t=1])
        eta_Wf[n] = mu_f[idf[n]] +   (psi_1 + psi_f[idf[n]])*(eta_Wm[n-1] - eta_meanf[n-1]);
        
        //average individual eta
        //male average eta[t=2] = mu_j + psi_j*eta_meanK[t=1]
        meaneta_m[n] = mu_m[idm[n]] +    (psi_1 + psi_m[idm[n]])*eta_meanm[n-1];
        
        //female average eta[t=2] = mu_k + psi_k*eta_meanJ[t=1]
        meaneta_f[n] = mu_f[idf[n]] +   (psi_1 + psi_f[idf[n]])*eta_meanf[n-1];
        
        //average eta males' partners [t=1] = mu_meanK + psi_meanK*mean eta_j[t-1]
        eta_meanm[n] = mu_meanm[idm[n]] +  (psi_1 + psi_meanm[idm[n]])*meaneta_m[n-1];
        
        //female average partner eta
        eta_meanf[n] = mu_meanf[idf[n]] +  (psi_1 + psi_meanf[idf[n]])*meaneta_f[n-1]; 
      }
    
    //add global intercept and between-individual parameters to linear predictor
    //other fixed effects can also be added here
    linpred_m[n] = alpha_0 + eta_Wm[n] + beta_B1*mu_meanm[idm[n]] + beta_B2*psi_meanm[idm[n]];
    linpred_f[n] = alpha_0 + eta_Wf[n] + beta_B1*mu_meanf[idf[n]] + beta_B2*psi_meanf[idf[n]];

    //residual trait values
    if(time[n]==1)
      {
        zeta_m [n] = AG_m[n] - linpred_m[n];
        zeta_f [n] = AG_f[n] - linpred_f[n];
      }
    else
      {
       linpred_m[n] = linpred_m[n] + phi * zeta_f[n-1];
       zeta_m[n] = AG_m[n] - linpred_m[n]; 
       
       linpred_f[n] = linpred_f[n] + phi * zeta_m[n-1];
       zeta_f[n] = AG_f[n] - linpred_f[n]; 
      }
  
    //correlated residuals between partners
    [zeta_m[n],zeta_f[n]]' ~ multi_normal_cholesky([0,0], diag_pre_multiply(sd_R, LR));
  }
    
  //fitness response model
  w_pred = nu_0 + beta_N1*mu_m[idmw] + beta_N2*psi_m[idmw] + beta_S1*mu_f[idfw] + beta_S2*psi_f[idfw] + 
                  beta_D1*(mu_m[idmw].*mu_f[idfw]) + beta_D2*(psi_m[idmw].*psi_f[idfw]);
  
  w ~ normal(w_pred, sd_delta);

  //model priors
  
  //fixed effects
  alpha_0 ~ std_normal();
  nu_0 ~ std_normal();
  psi_1 ~ std_normal();
  beta_B1 ~ std_normal();
  beta_B2 ~ std_normal();
  beta_N1 ~ std_normal();
  beta_N2 ~ std_normal();
  beta_S1 ~ std_normal();
  beta_S2 ~ std_normal();
  beta_D1 ~ std_normal();
  beta_D2 ~ std_normal();
  phi ~ std_normal();
  
  //random effects
  to_vector(sd_P) ~ cauchy(0,1);
  to_vector(sd_R) ~ cauchy(0,1);
  sd_delta ~ cauchy(0,1);
  LP ~ lkj_corr_cholesky(2);
  LR ~ lkj_corr_cholesky(2);
  to_vector(std_devP) ~ std_normal();
}

generated quantities{
matrix[2,2] Pcor = LP * LP'; //phenotypic SRN correlation matric
matrix[2,2] Rcor = LR * LR'; //residual correlation matrix
matrix[2,2] Pcov = diag_matrix(sd_P)*Pcor*diag_matrix(sd_P); //phenotypic SRN covariance
matrix[2,2] Rcov = diag_matrix(sd_R)*Rcor*diag_matrix(sd_R); //residual covariance 
vector<lower=0>[2] V_P = sd_P .* sd_P; //SRN phenotypic variances
vector<lower=0>[2] V_R = sd_R .* sd_R; //residual variances
}
