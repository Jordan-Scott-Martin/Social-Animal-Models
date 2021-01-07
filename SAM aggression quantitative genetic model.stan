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
  matrix[I,I] A; //relatedness matrix
  vector[2] AG[N_sex]; //combined male and female AG responses
  real w[Idyad]; //FS dyad response
  real time[N_sex];
}

transformed data{
  matrix[I,I] LA = cholesky_decompose(A); //lower-triangle A matrix
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
  
  
  //random effects (standard deviations)
  vector<lower=0, upper = 1>[2] sd_P; //phenotypic SRN mu & psi SDs
  vector<lower=0, upper = 1>[2] sd_R; //male & female residual SDs
  real<lower=0, upper = 1>sd_delta; //residual of fitness

  cholesky_factor_corr[2] LG; //genetic SRN correlations
  cholesky_factor_corr[2] LE; //permanent environmental SRN correlations
  cholesky_factor_corr[2] LR; //sex-specific residual correlations
  
  matrix[I,2] std_devG; //individual-level unscaled G SRN deviations
  matrix[I,2] std_devE; //individual-level unscaled E SRN deviations
  
  //SRN heritability parmameters, i.e. Var(G_RN) / Var(P_RN) 
  //see supplementary appendix SI for further explanation of this parameter
  vector<lower=0,upper=1>[2] SRN_h2;
  
}

transformed parameters {
  vector<lower=0>[2] sd_G; //SDs of G effects (derived from sd_P)
  vector<lower=0>[2] sd_E; //SDs of E effects (derived from sd_P)
  
  matrix[I,2] SRN_P; //scaled P SRN parameter deviations
  matrix[I,2] SRN_G; //scaled G SRN parameter deviations
  matrix[I,2] SRN_E; //scaled E SRN parameter deviations
  
  matrix[If, 2] partner_meanm; //average SRN parameters of males' partners
  matrix[Im, 2] partner_meanf; //average SRN parameters of females' partners
  
  //standard deviations of genetic effects
  //simplified from sqrt ( total RN phenotype variance * h2 )
  sd_G[1] = sd_P[1] * sqrt(SRN_h2[1]); //genetic SD for RN intercepts 
  sd_G[2] = sd_P[2] * sqrt(SRN_h2[2]);  //genetic SD for RN slopes
  
  //standard deviations of environmental effects (total phenotype SD * proportion environment SD)
  sd_E[1] = sd_P[1] * sqrt(1 - SRN_h2[1]); //environment SD for RN intercepts 
  sd_E[2] = sd_P[2] * sqrt(1 - SRN_h2[2]); //environment SD for RN slopes 
  
  //matrix normal parameterization of Kronecker product between G and A
  SRN_G =  LA * std_devG  * diag_pre_multiply(sd_G, LG)' ; 
  
  //non-centered parameterization of permanent environmental effects
  SRN_E = std_devE * diag_pre_multiply(sd_E, LE);
  
  //phenotypic RN effects (P = G + E); here G = additive genetic effects
  SRN_P = SRN_G + SRN_E;

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
  vector[Im] mu_Pm = col(SRN_Pm,1); //SRN intercepts
  vector[If] mu_Pf = col(SRN_Pf,1); 
  vector[Im] psi_Pm = col(SRN_Pm,2); //SRN slopes
  vector[If] psi_Pf = col(SRN_Pf,2); 
  
  //total SRN values (population average + deviation)
  vector[Im] mu_m = alpha_0 + mu_Pm;
  vector[If] mu_f = alpha_0 + mu_Pf; 
  vector[Im] psi_m = psi_1 + psi_Pm; 
  vector[If] psi_f = psi_1 + psi_Pf; 
  
  //separate mean partner SRN intercepts and slopes (deviations)
  vector[Im] mu_meanPm = col(partner_meanm,1); //mean partner SRN intercept for males
  vector[If] mu_meanPf = col(partner_meanf,1); //...for females
  vector[Im] psi_meanPm = col(partner_meanm,2); //mean partner SRN slope for males
  vector[If] psi_meanPf = col(partner_meanf,2); //...for females
  
  //mean partner total SRN values (population average + deviation)
  vector[Im] mu_meanm = alpha_0 + mu_meanPm;
  vector[If] mu_meanf = alpha_0 + mu_meanPf;
  vector[Im] psi_meanm = psi_1 + psi_meanPm; 
  vector[If] psi_meanf = psi_1 + psi_meanPf; 
  
  //initialize vectors for constructing individual-centered linear predictors
  vector[N_sex] eta_Wm; //within-individual centered male SRN trait value 
  vector[N_sex] eta_Wf; //within-individual centered female SRN trait value
  
  vector[N_sex] meaneta_m; //individual male SRN trait value toward average partner 
  vector[N_sex] meaneta_f; //individual female SRN trait toward average partner
  vector[N_sex] eta_meanm; //average SRN partner values for males
  vector[N_sex] eta_meanf; //average SRN partner values for females
  
  vector[N_sex] linpred_m; //expected value for male responses
  vector[N_sex] linpred_f; //expected value for female responses
  vector[2] linpred[N_sex]; //combined vector of male and female linear predictors
  vector[Idyad] w_pred; //linear predictor of fitness
  
  //predict unbiased SRN trait values at t=1 and t=2
  for (n in 1:N_sex) {
    
    //assumes that n = 1 in the context of an ongoing social interaction
    //if n = 1 prior to social context, then specify eta[t=1] = mu_j instead
    if (time[n]==1)
      {
        //within-individual centered eta
        //male eta[t=1] = mu_j + psi_j*(mu_k - mu_meanK)
        eta_Wm[n] = mu_m[idm[n]] +  psi_m[idm[n]]*(mu_f[idf[n]] - mu_meanm[idm[n]]) ;
        
        //female eta[t=1] = mu_k + psi_k*(mu_j - mu_meanJ)
        eta_Wf[n] = mu_f[idf[n]] + psi_f[idf[n]]*(mu_m[idm[n]] - mu_meanf[idf[n]]);
        
        //average individual eta
        //male eta[t=1] = mu_j + psi_j*mu_k
        meaneta_m[n] = mu_m[idm[n]] +  psi_m[idm[n]]*mu_meanm[idm[n]];
        
        //female eta[t=1] = mu_k + psi_k*mu_j
        meaneta_f[n] = mu_f[idf[n]] + psi_f[idf[n]]*mu_meanf[idf[n]];
        
        //average partner eta[t=1]
        //average eta males' partners [t=1] = mu_meanK + psi_meanK*mu_j
        eta_meanm[n] = mu_meanm[idm[n]] + psi_meanm[idm[n]]*mu_m[idm[n]];
        
        //average eta females' partners [t=1] = mu_meanJ + psi_meanJ*mu_k
        eta_meanf[n] = mu_meanf[idf[n]] + psi_meanf[idf[n]]*mu_f[idf[n]];
      }    
    else
      {
        //within-individual centered eta
        //male eta[t=2] = mu_j + psi_j*(eta_k[t=1] - eta_meanK[t=1])
        eta_Wm[n] = mu_m[idm[n]] + psi_m[idm[n]]*(eta_Wf[n-1] - eta_meanm[n-1]);
        
        //female eta[t=2] = mu_k + psi_k*(eta_j[t=1] - eta_meanJ[t=1])
        eta_Wf[n] = mu_f[idf[n]] + psi_f[idf[n]]*(eta_Wm[n-1] - eta_meanf[n-1]);
        
        //average individual eta
        //male average eta[t=2] = mu_j + psi_j*eta_meanK[t=1]
        meaneta_m[n] = mu_m[idm[n]] +  psi_m[idm[n]]*eta_meanm[n-1];
        
        //female average eta[t=2] = mu_k + psi_k*eta_meanJ[t=1]
        meaneta_f[n] = mu_f[idf[n]] + psi_f[idf[n]]*eta_meanf[n-1];
        
        //average eta males' partners [t=1] = mu_meanK + psi_meanK*mean eta_j[t-1]
        eta_meanm[n] = mu_meanm[idm[n]] + psi_meanm[idm[n]]*meaneta_m[n-1];
        
        //female average partner eta
        eta_meanf[n] = mu_meanf[idf[n]] + psi_meanf[idf[n]]*meaneta_f[n-1]; 
      }
    
    //add between-individual parameters to linear predictor
    linpred_m[n] = eta_Wm[n] + beta_B1*mu_meanPm[idm[n]] + beta_B2*psi_meanPm[idm[n]];
    linpred_f[n] = eta_Wf[n] + beta_B1*mu_meanPf[idf[n]] + beta_B2*psi_meanPf[idf[n]];
    
    //put male and female linear predictors together
    linpred[n] = [linpred_m[n], linpred_f[n]]';
    }
  
  //aggression response model with correlated residuals
  AG ~ multi_normal_cholesky(linpred, diag_pre_multiply(sd_R, LR));
    
  //fitness response model
  w_pred = nu_0 + beta_N1*mu_Pm[idmw] + beta_N2*psi_Pm[idmw] + beta_S1*mu_Pf[idfw] + beta_S2*psi_Pf[idfw] + 
                  beta_D1*(mu_Pm[idmw].*mu_Pf[idfw]) + beta_D2*(psi_Pm[idmw].*psi_Pf[idfw]);
  
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
  
  //random effects
  to_vector(sd_P) ~ cauchy(0,1);
  to_vector(sd_R) ~ cauchy(0,1);
  sd_delta ~ cauchy(0,1);
  
  LG ~ lkj_corr_cholesky(2);
  LE ~ lkj_corr_cholesky(2);
  LR ~ lkj_corr_cholesky(2);
  
  to_vector(std_devG) ~ std_normal();
  to_vector(std_devE) ~ std_normal();
  
  //reaction norm heritability
  to_vector(SRN_h2) ~ beta(1.2,1.2);
  
}


generated quantities{
//cor and cov matrices of SRN parameters and residuals
matrix[2,2] Gcor = LG * LG'; //G SRN correlation matric
matrix[2,2] Ecor = LE * LG'; //E SRN correlation matric
matrix[2,2] Rcor = LR * LR'; //residual correlation matrix

matrix[2,2] Rcov = diag_matrix(sd_R)*Rcor*diag_matrix(sd_R); //residual covariance
matrix[2,2] Gcov = diag_matrix(sd_G)*Gcor*diag_matrix(sd_G); //G SRN covariance
matrix[2,2] Ecov = diag_matrix(sd_E)*Ecor*diag_matrix(sd_E); //E SRN covariance
matrix[2,2] Pcov = Gcov + Ecov; //P SRN covariance
matrix[2,2] Pcor = inverse(diag_matrix(sd_P))*Pcov*inverse(diag_matrix(sd_P)); //P SRN correlation

//variances
vector<lower=0>[2] V_P = sd_P .* sd_P;
vector<lower=0>[2] V_G = sd_G .* sd_G;
vector<lower=0>[2] V_E = sd_E .* sd_E;
vector<lower=0>[2] V_R = sd_R .* sd_R;



}
