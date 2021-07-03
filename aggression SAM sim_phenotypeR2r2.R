#This code will simulate the desired number of datasets,
#and run the aggression SAM model on each sample
#Performance metrics are collected in a for loop
#and saved to the set directory. The estimation time of a single model
#at N=300 may take 20-30 minutes or so. The full simulation may therefore
#take many hours or days depending on the number of datasets and
#parallel processing capacity of the computer.

#####################################################################
#Prepare workspace
#####################################################################

memory.limit(2e5)

#load required packages
library(Matrix); library(MASS); library(mvtnorm); library(dplyr)
library(RANN); library(MCMCglmm)

#stan settings for compilation and parallel processing
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

#set directory for saving results
setwd("C:/Users/jormar/Dropbox/JEB revision 2/agg simulation/after acceptance edits")

#####################################################################
#Parameter values
#####################################################################

#population effects
alpha_0 = 0 #aggression intercept
nu_0 = 1 #fitness intercept
psi_1 = 0.3 #interaction coefficient
beta_n1 = 0.3 #selection gradients
beta_n2 = -0.3
beta_s1 = 0.3
beta_s2 = -0.3
beta_d1 = -0.3
beta_d2 = -0.3
phi = 0.3 #residual feedback parameter

#variance
V = 0.3 #variance random effects
res_V = 1 #residual variance

#associations
r_G =  0 #correlation between additive genetic RN components
r_E =  0 #cor between environmental RN components
r_R =   0.3 #cor between residuals
r_alpha = 0.3 #correlation between partners, equivalent to assortment coefficient
#when there is constant variance between focal and social vectors

#selection differentials and response
Pcov = matrix(c(0.6,0,0,0.6),2,2)
Gcov = matrix(c(0.3,0,0,0.3),2,2)
Beta_N = matrix(c(0.3,-0.3),2,1)
Beta_S = matrix(c(0.3,-0.3),2,1)
Beta_alpha = matrix(c(0,0,0,r_alpha),2,2)
    
true_differential = Pcov %*% Beta_N + diag(diag(Pcov),2) %*%Beta_alpha %*% Beta_S
true_response = Gcov %*% Beta_N + diag(diag(Gcov),2) %*%Beta_alpha %*% Beta_S

#####################################################################
#Custom functions
#####################################################################

#simulate genetic correlation function (A matrix)
pedfun<-function(popmin, popmax, ngenerations,
                 epm, nonb, nids, I, missing=FALSE){
  
  # get list of individuals and their generations
  gener<-1:ngenerations
  
  genern <- rep(1:ngenerations, times = nids)
  ID <- 1:sum(nids)
  
  # runs on generation-by-generation basis
  for(i in 1:ngenerations){
    
    id<-ID[which(genern==i)]
    dam<-rep(NA, nids[i])
    sire<-rep(NA, nids[i])

    # randomly allocates sex (0 = male, 1 = female)
    sex<-sample(c(0,1), length(id), replace=TRUE)
    
    # for first generation, no dams or sires are known 
    # so remain NA
    
    if(i==1){

      # combine into single data frame
      pedigree<-data.frame(id=id, dam=dam, sire=sire, 
                           generation=i, sex=sex)
      
    }
    
    else if(i>1){
   
      # for all generations after first
      # list of all possible dams and sires
      # from previous generation
      pdams<-pedigree$id[which(pedigree$generation==(i-1) &
                                 pedigree$sex==1)]
      psires<-pedigree$id[which(pedigree$generation==(i-1) &
                                  pedigree$sex==0)]
      
      # determine number of pairs
      # depending on how many males and females
      # and the proportion of the population that is non-breeding
      npairs<-min(length(pdams), length(psires)) - 
        round(min(length(pdams), length(psires))*nonb)
      
      # selects breeding males and females
      pdams<-pedigree$id[which(pedigree$generation==(i-1) & 
                                 pedigree$sex==1)]
      psires<-pedigree$id[which(pedigree$generation==(i-1) & 
                                  pedigree$sex==0)]
      
      if(length(pdams)<npairs | length(psires)<npairs){
        npairs<-min(length(pdams), length(psires))
      }
      
      # selects pairs from possible dams and sires
      pairs<-data.frame(dam=sample(pdams, npairs, replace=FALSE),
                        sire=sample(psires, npairs, replace=FALSE))
      # gives each offspring their parental pair
      pairid<-as.numeric(sample(rownames(pairs), 
                                length(id), replace=TRUE))
      
      # gives each offspring their sex
      sex<-sample(c(0,1), length(id), replace=TRUE)
      
      # put into dataframe format
      addped<-data.frame(id=id, 
                         dam=pairs$dam[pairid], 
                         sire=pairs$sire[pairid],
                         generation=i, 
                         sex=sex)
      
      
      # deals with extra-pair mating (if included)
      if(!is.null(epm)){
        
        # for each individual, sample if they are extra pair
        # if 0 not extra pair
        # if 1 sire resampled from breeding population
        # if 2 dam resampled
        ext<-sample(c(0,1,2), nrow(addped), 
                    replace=TRUE, 
                    prob = c(1-epm, epm/2, epm/2))
        for(j in 1:nrow(addped)){
          if(ext[j]>0){
            if(ext[j]==1){
              addped$sire[j]<-sample(psires,1,replace=TRUE)
            }else if (ext[j]==2){
              addped$dam[j]<-sample(pdams,1,replace=TRUE)
            }
            
          }
        }
      }
      
      
      # add new generation to the whole pedigree
      pedigree<-rbind(pedigree, addped)
    }
    
    }
  
  ped <- pedigree
  
  # make id's non-numeric
  ped$id<-paste("ID",ped$id, sep="")
  ped$dam[which(!is.na(ped$dam))]<-paste("ID",ped$dam[which(!is.na(ped$dam))], sep="")
  ped$sire[which(!is.na(ped$sire))]<-paste("ID",ped$sire[which(!is.na(ped$sire))], sep="")
  ped$id<-as.character(ped$id)
  ped$dam<-as.character(ped$dam)
  ped$sire<-as.character(ped$sire)

  IDs <- sample(ped[ped$generation==ngenerations, "id"], I, replace=FALSE)
  ped <- MCMCglmm::prunePed(ped, keep = IDs, make.base=TRUE)
  inv.phylo <- MCMCglmm::inverseA(ped[,c("id","dam","sire")])
  A <- solve(inv.phylo$Ainv)
  A <- cov2cor(A)
  A = (A + t(A))/2 # Not always symmetric after inversion
  A <- as.matrix(A)
  rownames(A) <- rownames(inv.phylo$Ainv)
  colnames(A) <- rownames(inv.phylo$Ainv)
  
  #subset to final generation
  A_sub<-A[IDs,IDs]
  A_mat <- as.matrix(nearPD(A_sub)$mat)
  A_mat <- cov2cor(A_mat)
  return(A_mat)
}

#simulate single dataset for aggression model
simSAM<-function(){
  
  #####################################################################
  #Simulate aggression and fitness across male and female partners
  #####################################################################
  
  #generate pedigree
  
  #population properties
  popmin=400
  popmax=600
  ngenerations = 10
  nids<-sample(popmin:popmax, ngenerations, replace=TRUE) #N / generation
  epm = sample(seq(0.15, 0.25,by=0.05),1) #extra-pair mating
  nonb = sample(seq(0.4,0.6,by=0.05),1) #proportion of non-breeding / generation
  
  #relatedness matrix
  A_mat <- pedfun(popmin=popmin, popmax=popmax, ngenerations=ngenerations,
                  epm=epm, nonb=nonb, nids=nids, I=I, missing=FALSE)
  
  #Random effect correlations
  G_cor <- matrix(c(1,r_G,r_G,1), nrow=2, ncol=2) #A0, A1
  G_sd  <- c(sqrt(V),sqrt(V)) #G effect sds
  G_cov <- diag(G_sd) %*% G_cor %*% diag(G_sd)
  
  E_cor <- matrix(c(1,r_E,r_E,1), nrow=2, ncol=2) #E0, E1
  E_sd  <- c(sqrt(V),sqrt(V)) #E effect sds
  E_cov <- diag(E_sd) %*% E_cor %*% diag(E_sd)
  
  #block matrices for simulation
  G_block <- G_cov %x% A_mat
  E_block <- E_cov %x% diag(1,I)
  
  #generate correlated RNs
  Gvalues <- rmvnorm(1, mean=rep(0,I*2), sigma=G_block)
  G_val = data.frame(matrix(Gvalues, nrow=I, ncol=2))
  cor(G_val)
  
  Evalues <- rmvnorm(1, mean=rep(0,I*2), sigma=E_block)
  E_val = data.frame(matrix(Evalues, nrow=I, ncol=2))
  cor(E_val)
  
  #combine
  P = cbind(G_val,E_val)
  colnames(P) = c("A0", "A1", "E0", "E1")
  P$ID = seq(1:I)
  
  #calculate phenotypic reaction norms
  P$P0 = P$A0 + P$E0
  P$P1 = P$A1 + P$E1
  
  #####################################################################
  #Partner assortment within each breeding season
  #####################################################################

  pairs = list()
  for (j in 1:I_partner){
    #male additive genetic RN slopes (x I_partner for multiple lifetime partners)
    sort.m <- data.frame(P1_m = P$P1[1:(I/2)], ID_m = (1:(I/2)) )
    sort.m<-sort.m[order(sort.m[,"P1_m"]),]
    #female phenotypic RN slopes
    sort.f <- data.frame(P1_f = P$P1[(I/2 + 1):I], ID_f = ((I/2+1):I) )
    sort.f<-sort.f[order(sort.f[,"P1_f"]),]
    #generate random dataset with desired rank-order correlation
    temp_mat <- matrix(r_alpha, ncol = 2, nrow = 2) #cor of male and female values
    diag(temp_mat) <- 1 #cor matrix
    #sim values
    temp_data1<-MASS::mvrnorm(n = I/2, mu = c(0, 0), Sigma = temp_mat, empirical=TRUE)
    #ranks of random data
    rm <- rank(temp_data1[ , 1], ties.method = "first")
    rf <- rank(temp_data1[ , 2], ties.method = "first")
    #induce cor through rank-ordering of RN vectors
    cor(sort.m$P1_m[rm], sort.f$P1_f[rf])
    #sort partner ids into dataframe
    partner.id = data.frame(ID_m = sort.m$ID_m[rm], ID_f = sort.f$ID_f[rf])
    partner.id = partner.id[order(partner.id[,"ID_m"]),]
    #add to list
    pairs[[j]] = partner.id
  }
  
  partner.id = bind_rows(pairs)
  partner.id[order(partner.id$ID_m),]
  
  #put all dyads together
  partner.id$dyadn = seq(1:nrow(partner.id))
  
  #add values back to dataframe (male and joint)
  partner.id$P0m <- P$P0[match(partner.id$ID_m,P$ID)]
  partner.id$P0f <- P$P0[match(partner.id$ID_f,P$ID)]
  partner.id$P1m <- P$P1[match(partner.id$ID_m,P$ID)]
  partner.id$P1f <- P$P1[match(partner.id$ID_f,P$ID)]
  
  partner.id$A0m <- P$A0[match(partner.id$ID_m,P$ID)]
  partner.id$A0f <- P$A0[match(partner.id$ID_f,P$ID)]
  partner.id$A1m <- P$A1[match(partner.id$ID_m,P$ID)]
  partner.id$A1f <- P$A1[match(partner.id$ID_f,P$ID)]
  
  partner.id$E0m <- P$E0[match(partner.id$ID_m,P$ID)]
  partner.id$E0f <- P$E0[match(partner.id$ID_f,P$ID)]
  partner.id$E1m <- P$E1[match(partner.id$ID_m,P$ID)]
  partner.id$E1f <- P$E1[match(partner.id$ID_f,P$ID)]
  
  #calculate mean partner phenotype for each subject
  
  #average female for male partners
  mean_0m <- aggregate(P0f ~ ID_m, mean, data = partner.id)
  names(mean_0m)[2] <- "meanP0m"
  
  mean_1m <- aggregate(P1f ~ ID_m, mean, data = partner.id)
  names(mean_1m)[2] <- "meanP1m"
  
  partner.id$meanP0m <- mean_0m$meanP0m[match(partner.id$ID_m,mean_0m$ID_m)]
  partner.id$meanP1m <- mean_1m$meanP1m[match(partner.id$ID_m,mean_1m$ID_m)]
  
  #average male for female partners
  mean_0f <- aggregate(P0m ~ ID_f, mean, data = partner.id)
  names(mean_0f)[2] <- "meanP0f"
  
  mean_1f <- aggregate(P1m ~ ID_f, mean, data = partner.id)
  names(mean_1f)[2] <- "meanP1f"
  
  partner.id$meanP0f <- mean_0f$meanP0f[match(partner.id$ID_f,mean_0f$ID_f)]
  partner.id$meanP1f <- mean_1f$meanP1f[match(partner.id$ID_f,mean_1f$ID_f)]
  
  #number of dyads
  ndyad = nrow(partner.id)
  
  #expand for repeated measures
  partner.id$rep <- I_obs
  pair_df <- partner.id[rep(row.names(partner.id), partner.id$rep),]
  
  #correlations
  cor(partner.id$P0m, partner.id$P0f)
  cor(partner.id$P1m, partner.id$P0f)
  cor(partner.id$P0m, partner.id$P1f)
  cor(partner.id$P1m, partner.id$P1f)
  
  #####################################################################
  #Additional effects
  #####################################################################
  
  #correlated  residuals between male and females
  R_cor <- matrix(c(1,r_R,r_R,1), nrow=2, ncol=2)
  res_sd <- sqrt(res_V)
  R_cov <- diag(c(res_sd,res_sd)) %*% R_cor %*% diag(c(res_sd,res_sd))
  res_ind<-data.frame(rmvnorm(nrow(pair_df), c(0,0), R_cov))
  pair_df$resAGm = res_ind$X1
  pair_df$resAGf = res_ind$X2
  
  #####################################################################
  #Simulate responses over t = {1,2} per partner
  #####################################################################
  
  #add interaction number
  pair_df$turn = rep(c(1,2),ndyad) 
  
  #average male social environment at time = 1
  pair_df[pair_df$turn==1,"meaneta_m"] =  pair_df[pair_df$turn==1,"meanP0m"] + 
                                         (psi_1 + pair_df[pair_df$turn==1,"meanP1m"])*(pair_df[pair_df$turn==1,"P0m"])
  
  #average female social environment at time = 1
  pair_df[pair_df$turn==1,"meaneta_f"] =  pair_df[pair_df$turn==1,"meanP0f"] + 
                                         (psi_1 + pair_df[pair_df$turn==1,"meanP1f"])*(pair_df[pair_df$turn==1,"P0f"])
  

  #individual prediction at t = 1
  
  #males
  #eta_j{t=1} = mu_j + psi_j*(mu_k - mu_meanK)
  pair_df[pair_df$turn==1,"eta_m"] = pair_df[pair_df$turn==1,"P0m"] + 
    (psi_1 + pair_df[pair_df$turn==1,"P1m"])*(pair_df[pair_df$turn==1,"P0f"]-pair_df[pair_df$turn==1,"meanP0m"])
  #females
  #eta_k{t=1} = mu_k + psi_k*(mu_j - mu_meanJ)                                                                                                                
  pair_df[pair_df$turn==1,"eta_f"] = pair_df[pair_df$turn==1,"P0f"] + 
    (psi_1 + pair_df[pair_df$turn==1,"P1f"])*(pair_df[pair_df$turn==1,"P0m"]-pair_df[pair_df$turn==1,"meanP0f"])
  
  #individual prediction at t = 2
  
  #eta_j{t=2} = mu_j + psi_j*(eta_k{t=1} - eta_meanK{t=1})
  pair_df[pair_df$turn==2,"eta_m"] = pair_df[pair_df$turn==2,"P0m"] +
    (psi_1 + pair_df[pair_df$turn==2,"P1m"])*(pair_df[pair_df$turn==1,"eta_f"]-pair_df[pair_df$turn==1,"meaneta_m"])
  
  #females
  pair_df[pair_df$turn==2,"eta_f"] = pair_df[pair_df$turn==2,"P0f"] + 
    (psi_1 + pair_df[pair_df$turn==2,"P1f"])*(pair_df[pair_df$turn==1,"eta_m"]-pair_df[pair_df$turn==1,"meaneta_f"])
  
  #individual prediction at t = 1
  
  #males
  #eta_j{t=1} = mu_j + psi_j*(mu_k - mu_meanK)
  pair_df[pair_df$turn==1,"eta_m"] = pair_df[pair_df$turn==1,"P0m"] + 
    (psi_1 + pair_df[pair_df$turn==1,"P1m"])*(pair_df[pair_df$turn==1,"P0f"])
  #females
  #eta_k{t=1} = mu_k + psi_k*(mu_j - mu_meanJ)                                                                                                                
  pair_df[pair_df$turn==1,"eta_f"] = pair_df[pair_df$turn==1,"P0f"] + 
    (psi_1 + pair_df[pair_df$turn==1,"P1f"])*(pair_df[pair_df$turn==1,"P0m"])
  
  #individual prediction at t = 2
  
  #eta_j{t=2} = mu_j + psi_j*(eta_k{t=1} - eta_meanK{t=1})
  pair_df[pair_df$turn==2,"eta_m"] = pair_df[pair_df$turn==2,"P0m"] +
    (psi_1 + pair_df[pair_df$turn==2,"P1m"])*(pair_df[pair_df$turn==1,"eta_f"])
  
  #females
  pair_df[pair_df$turn==2,"eta_f"] = pair_df[pair_df$turn==2,"P0f"] + 
    (psi_1 + pair_df[pair_df$turn==2,"P1f"])*(pair_df[pair_df$turn==1,"eta_m"])
  
  
  #add intercept and residual
  pair_df$AG_m = alpha_0 + pair_df$eta_m + pair_df$resAGm
  pair_df$AG_f = alpha_0 + pair_df$eta_f + pair_df$resAGf
  
  #add residual feedback
  pair_df[pair_df$turn==2,"AG_m"] = pair_df[pair_df$turn==2,"AG_m"] + phi * pair_df[pair_df$turn==1,"resAGf"]
  pair_df[pair_df$turn==2,"AG_f"] = pair_df[pair_df$turn==2,"AG_f"] + phi * pair_df[pair_df$turn==1,"resAGm"]

  
  #####################################################################
  
  #dyad fitness (nu_0 = 1 so that w is relative fitness w = W/W_mean for the unbiased population mean fitness)
  pair_df$w_mu = nu_0 + beta_n1*pair_df$P0m + beta_n2*pair_df$P1m + beta_s1*pair_df$P0f + beta_s2*pair_df$P1f +
                              beta_d1*(pair_df$P0m*pair_df$P0f) + beta_d2*(pair_df$P1m*pair_df$P1f)
  
  w_mu<-pair_df[seq(1, nrow(pair_df), by=I_obs),"w_mu"]
  
  #add stochastic effects
  w = w_mu + rnorm(length(w_mu),0, res_sd)
  
  #####################################################################
  #Prepare data for Stan
  #####################################################################
  
  #individual indices
  idm<-pair_df$ID_m #male ID
  idf<-pair_df$ID_f #female ID
  idf<-idf - (Im) #index within female vector
  dyadAG <- pair_df$dyadn
  dyadw <- seq(1:ndyad)
  obs <- seq(1:nrow(pair_df))
  
  #partner IDs for male individuals
  partners_m<-data.frame(idfocal = rep(1:(I/2)), #all partners ID
                   partner1 = NA, partner2 = NA, partner3 = NA, partner4 = NA)
                   for(i in 1:(I/2)){partners_m[i,c(2:5)] <-partner.id[partner.id$ID_m==i,"ID_f"]}
  
  #partner IDs for female individuals
  partners_f<-data.frame(idfocal = rep((I/2+1):I), #all partners ID
                   partner1 = NA, partner2 = NA, partner3 = NA, partner4 = NA)
                   for(i in (I/2+1):I){partners_f[i-(I/2),c(2:5)] <-partner.id[partner.id$ID_f==i,"ID_m"]}
  
  ######################
  
  #data prep for Stan
  stan_data <-
    list(N_sex = N_sex, I = I, Im=Im, If = If, Idyad=ndyad, idm = idm, idf = idf,
         idmw = partner.id$ID_m, idfw = c(partner.id$ID_f-Im),
         dyadw = dyadw, partners_m = partners_m, partners_f = partners_f,
         AG_m = pair_df$AG_m, AG_f = pair_df$AG_f, time = pair_df$turn, w = w, A = A_mat)
  
  return(stan_data)
  
}

#####################################################################
#Simulate datasets
#####################################################################

nD = 200 #number of datasets

#common settings
I_partner = 4 #partners/individual
I_obs = 2 #observations/individual/seasonal partner
I_sample = I_partner*I_obs #samples/individual
P = 2 #SRN random effects (intercept, slope)

#initialize lists for dataframes
data1 <- vector("list", length=nD)
data2 <- vector("list", length=nD)
data3 <- vector("list", length=nD)

#simulate datasets (will take some time)

I = 100 #total sample
Im = I/2 #total males
If = I/2 #total females
N <- I*I_sample #total observations
N_sex <- N/2 #total observations/sex

for (i in 1:(nD)) {
  data1[[i]] = simSAM() }

I = 200
Im = I/2 #total males
If = I/2 #total females
N <- I*I_sample
N_sex <- N/2
for (i in 1:(nD)) {
  data2[[i]] = simSAM() }

I = 300 
Im = I/2 #total males
If = I/2 #total females
N <- I*I_sample 
N_sex <- N/2 
for (i in 1:(nD)) {
  data3[[i]] = simSAM() }

#save simulated datasets
saveRDS(data1, "data1.RDS")
saveRDS(data2, "data2.RDS")
saveRDS(data3, "data3.RDS")

#####################################################################
#Analyze data
#####################################################################

#MCMC settings
n_iter <- 1500
n_warm <- 1000
n_chains <- 4

#compile model
SAM <- rstan::stan_model(file="SAM aggression phenotypic model_R3.stan")

#number of datasets/data file to analyze
run = 200 #nD #i.e. all datasets per data file; set to 1 for demonstration

#####################################################################
#Data 1

  #load simulated datasets
  data1<-readRDS("data1.RDS")
  
  #initialize lists for results
  time1p <- vector("list", length=run)
  median1p <- vector("list", length=run)
  bias1p <- vector("list", length=run)
  UC1p <- vector("list", length=run)
  PP1p <- vector("list", length=run)
  Rhat1p <- vector("list", length=run)

  #run SAM across datasets (this will take awhile)
  counter = 0 #track progress
  for(i in 1:run){
    
    #single dataset
    df = data1[[i]]
    
    #run model
    start_time <- Sys.time()
    SAM_m <- rstan::sampling(SAM, data=df, init=0, iter=n_iter, warmup=n_warm, seed = i,
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.80, max_treedepth=10 ))
    end_time <- Sys.time()
    
    summary(SAM_m)
    
    #computation time
    time1p[[i]] = as.numeric(end_time - start_time)
    
    ###############################################################
    #Calculate additional posterior quantities
    ###############################################################
    
    #extract posteriors
    post <- rstan::extract(SAM_m)
    
    #temporary vectors for assortment coefficients
    SRN_PV = post$V_P
    SRN_Psd = post$sd_P
    SRN_PVmean = post$V_P / I_partner #expected variance for mean of partners
    SRN_Psdmean = sqrt(SRN_PVmean) #expected SD for mean of partners
    SRN_focal1 <- post$SRN_P[,,1] #individual intercepts
    SRN_focal2 <- post$SRN_P[,,2] #individual slopes
    SRN_partner1 <- cbind(post$partner_meanm[,,1], post$partner_meanf[,,1])
    SRN_partner2 <- cbind(post$partner_meanm[,,2], post$partner_meanf[,,2])
    
    #scale mean partner variance to variance of single partner
    #see appendix SI for details
    
    SRN_partner1s = SRN_partner1
    for(j in 1:nrow(SRN_partner1))
    {SRN_partner1s[j,] = ( SRN_partner1[j,] / SRN_Psdmean[j,1] ) * SRN_Psd[j,1] }
   
    SRN_partner2s = SRN_partner2
    for(j in 1:nrow(SRN_partner2))
    {SRN_partner2s[j,] = ( SRN_partner2[j,] / SRN_Psdmean[j,2] ) * SRN_Psd[j,2] }
    
    #assortment matrix
    Beta_alpha = list()
    
    #generate matrices across each posterior sample
    for(j in 1:nrow(SRN_focal1))
    {
            Beta_mat = matrix(NA,2,2)
            
            #mu' ~ mu
            Beta_mat[1,1] =  cov(SRN_focal1[j,], SRN_partner1s[j,])/var(SRN_focal1[j,])
            #mu' ~ psi
            Beta_mat[2,1] =  cov(SRN_focal2[j,], SRN_partner1s[j,])/var(SRN_focal2[j,])
            #psi' ~ mu                                                        
            Beta_mat[1,2] =  cov(SRN_focal1[j,], SRN_partner2s[j,])/var(SRN_focal1[j,])
            #psi' ~ psi
            Beta_mat[2,2] =  cov(SRN_focal2[j,], SRN_partner2s[j,])/var(SRN_focal2[j,])
            
            Beta_alpha[[j]] = Beta_mat
    } 

    #extract beta_psi'psi (assortment on social plasticity)
    Beta_psi = unlist(lapply(Beta_alpha, function(x) x[2,2]))
    median(Beta_psi)
    
    #generate other relevant matrices
    Beta_N = matrix(c(post$beta_N1,post$beta_N2),ncol=2)
    Beta_S = matrix(c(post$beta_S1,post$beta_S2),ncol=2)
    P = post$Pcov

    #selection differential

    #initialize dataframe
    s_SRN = data.frame(s_mu = rep(NA,nrow(Beta_N)), s_psi = rep(NA,nrow(Beta_N)))
   
    #populate with selection differentials
    for(j in 1:nrow(P)){
      s_SRN[j,] = P[j,,] %*% t(t(Beta_N[j,])) + diag(diag(P[j,,]),2,) %*% Beta_alpha[[j]] %*% t(t(Beta_S[j,])) }
    
    #organize posterior for ease of use
    postl <-
      list(P0V = post$V_P[,1],
           P1V = post$V_P[,2],
           RmV = post$V_R[,1],
           RfV = post$V_R[,2],
           Rcor = post$Rcor[,2,1],
           psi = post$psi_1,
           Beta_psi = Beta_psi,
           bN1 = post$beta_N1,
           bN2 = post$beta_N2,
           bS1 = post$beta_S1,
           bS2 = post$beta_S2,
           bD1 = post$beta_D1,
           bD2 = post$beta_D2,
           s_mu = s_SRN$s_mu,
           s_psi = s_SRN$s_psi,
           phi = post$phi)
    
    #median estimates
    median1p[[i]]<-list(
      P0V_med = median(postl$P0V),
      P1V_med = median(postl$P1V),
      RmV_med = median(postl$RmV),
      RfV_med = median(postl$RfV),
      Rcor_med = median(postl$Rcor),
      psi_med = median(postl$psi),
      Beta_psi_med = median(postl$Beta_psi),
      bN1_med = median(postl$bN1),
      bN2_med = median(postl$bN2),
      bS1_med = median(postl$bS1),
      bS2_med = median(postl$bS2),
      bD1_med = median(postl$bD1),
      bD2_med = median(postl$bD2),
      s_mu_med = median(postl$s_mu),
      s_psi_med = median(postl$s_psi),
      phi_med = median(postl$phi) )
    
    
    #bias
    bias1p[[i]]<-list(
      P0V_bias = (median(postl$P0V) - V*2) / (V*2) ,
      P1V_bias = (median(postl$P1V) - V*2) / (V*2) ,
      RmV_bias = (median(postl$RmV) - res_V) / (res_V),
      RfV_bias = (median(postl$RfV) - res_V) / (res_V),
      Rcor_bias = (median(postl$Rcor) - r_R ) / (r_R),
      psi_bias = (median(postl$psi) - psi_1 ) / (psi_1),
      Beta_psi_bias = (median(postl$Beta_psi) - r_alpha) / (r_alpha),
      bN1_bias = (median(postl$bN1) - beta_n1 ) / (beta_n1 ),
      bN2_bias = (median(postl$bN2) - beta_n2 ) / (beta_n2 ),
      bS1_bias = (median(postl$bS1) - beta_s1 ) / (beta_s1 ),
      bS2_bias = (median(postl$bS2) - beta_s2 ) / (beta_s2 ),
      bD1_bias = (median(postl$bD1) - beta_d1 ) / (beta_d1 ),
      bD2_bias = (median(postl$bD2) - beta_d2 ) / (beta_d2 ),
      s_mu_bias = (median(postl$s_mu) - true_differential[1] ) / (true_differential[1] ),
      s_psi_bias = (median(postl$s_psi) - true_differential[2]) / (true_differential[2] ),
      phi_bias = (median(postl$phi) - phi) / phi )
    
    #Uncertainty
    UC1p[[i]]<-list(
      P0V_uc = mad(postl$P0V)/median(postl$P0V),
      P1V_uc = mad(postl$P1V)/median(postl$P1V),
      RmV_uc = mad(postl$RmV)/median(postl$RmV),
      RfV_uc = mad(postl$RfV)/median(postl$RfV),
      Rcor_uc = mad(postl$Rcor)/median(postl$Rcor),
      psi_uc = mad(postl$psi)/median(postl$psi),
      Beta_psi_uc = mad(postl$Beta_psi)/median(postl$Beta_psi),
      bN1_uc = mad(postl$bN1)/median(postl$bN1),
      bN2_uc = mad(postl$bN2)/median(postl$bN2),
      bS1_uc = mad(postl$bS1)/median(postl$bS1),
      bS2_uc = mad(postl$bS2)/median(postl$bS2),
      bD1_uc = mad(postl$bD1)/median(postl$bD1),
      bD2_uc = mad(postl$bD2)/median(postl$bD2),
      s_mu_uc = mad(postl$s_mu)/median(postl$s_mu),
      s_psi_uc = mad(postl$s_psi)/median(postl$s_psi),
      phi_uc = (mad(postl$phi)/median(phi)) )
    
    #PP
    PP1p[[i]]<- list(
      Rcor_pp = sum(postl$Rcor>0)/length(postl$Rcor),
      psi_pp = sum(postl$psi>0)/length(postl$psi),
      Beta_psi_pp = sum(postl$Beta_psi>0)/length(postl$Beta_psi),
      bN1_pp = sum(postl$bN1>0)/length(postl$bN1),
      bN2_pp = sum(postl$bN2<0)/length(postl$bN2),
      bS1_pp = sum(postl$bS1>0)/length(postl$bS1),
      bS2_pp = sum(postl$bS2<0)/length(postl$bS2),
      bD1_pp = sum(postl$bD1<0)/length(postl$bD1),
      bD2_pp = sum(postl$bD2<0)/length(postl$bD2),
      s_mu_pp = sum(postl$s_mu>0)/length(postl$s_mu),
      s_psi_pp = sum(postl$s_psi<0)/length(postl$s_psi),
      phi_pp = sum(postl$phi>0)/length(postl$phi) )
    
    #Rhat
    summary = summary(SAM_m)$summary
    parnames <- c("alpha_0","nu_0","psi_1", "beta_N1", "beta_N2", "beta_S1", "beta_S2", "beta_D1", "beta_D2", 
                  "V_P[1]","V_P[2]", "V_R[1]", "V_R[2]","Rcor[1,2]", "phi")
    
    Rhat1p[[i]] = summary[parnames,"Rhat"]
    
    saveRDS(time1p, "time1p.RDS")
    saveRDS(median1p, "median1p.RDS")
    saveRDS(bias1p, "bias1p.RDS")
    saveRDS(UC1p, "UC1p.RDS")
    saveRDS(PP1p, "PP1p.RDS")
    saveRDS(Rhat1p, "Rhat1p.RDS")
    
    counter <- counter + 1; 
    print(paste(counter/run*100, "% has been processed"))
  }

#####################################################################
#Data 2

  #load simulated datasets
  data2<-readRDS("data2.RDS")
  
  #initialize lists for results
  time2p <- vector("list", length=run)
  median2p <- vector("list", length=run)
  bias2p <- vector("list", length=run)
  UC2p <- vector("list", length=run)
  PP2p <- vector("list", length=run)
  Rhat2p <- vector("list", length=run)

  #run SAM across datasets (this will take awhile)
  counter = 0 #track progress
  for(i in 1:run){
    
    #single dataset
    df = data2[[i]]
    
    #run model
    start_time <- Sys.time()
    SAM_m <- rstan::sampling(SAM, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i,
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.95, max_treedepth=10 ))
    end_time <- Sys.time()
    
    #computation time
    time2p[[i]] = as.numeric(end_time - start_time)
    
    ###############################################################
    #Calculate additional posterior quantities
    ###############################################################
    
        #extract posteriors
    post <- rstan::extract(SAM_m)
    
    #temporary vectors for assortment coefficients
    SRN_PV = post$V_P
    SRN_Psd = post$sd_P
    SRN_PVmean = post$V_P / I_partner #expected variance for mean of partners
    SRN_Psdmean = sqrt(SRN_PVmean) #expected SD for mean of partners
    SRN_focal1 <- post$SRN_P[,,1] #individual intercepts
    SRN_focal2 <- post$SRN_P[,,2] #individual slopes
    SRN_partner1 <- cbind(post$partner_meanm[,,1], post$partner_meanf[,,1])
    SRN_partner2 <- cbind(post$partner_meanm[,,2], post$partner_meanf[,,2])
    
    #scale mean partner variance to variance of single partner
    #see appendix SI for details
    
    SRN_partner1s = SRN_partner1
    for(j in 1:nrow(SRN_partner1))
    {SRN_partner1s[j,] = ( SRN_partner1[j,] / SRN_Psdmean[j,1] ) * SRN_Psd[j,1] }
   
    SRN_partner2s = SRN_partner2
    for(j in 1:nrow(SRN_partner2))
    {SRN_partner2s[j,] = ( SRN_partner2[j,] / SRN_Psdmean[j,2] ) * SRN_Psd[j,2] }
    
    #assortment matrix
    Beta_alpha = list()
    
    #generate matrices across each posterior sample
    for(j in 1:nrow(SRN_focal1))
    {
            Beta_mat = matrix(NA,2,2)
            
            #mu' ~ mu
            Beta_mat[1,1] =  cov(SRN_focal1[j,], SRN_partner1s[j,])/var(SRN_focal1[j,])
            #mu' ~ psi
            Beta_mat[2,1] =  cov(SRN_focal2[j,], SRN_partner1s[j,])/var(SRN_focal2[j,])
            #psi' ~ mu                                                        
            Beta_mat[1,2] =  cov(SRN_focal1[j,], SRN_partner2s[j,])/var(SRN_focal1[j,])
            #psi' ~ psi
            Beta_mat[2,2] =  cov(SRN_focal2[j,], SRN_partner2s[j,])/var(SRN_focal2[j,])
            
            Beta_alpha[[j]] = Beta_mat
    } 
    
    #extract beta_psi'psi (assortment on social plasticity)
    Beta_psi = unlist(lapply(Beta_alpha, function(x) x[2,2]))
    
    #generate other relevant matrices
    Beta_N = matrix(c(post$beta_N1,post$beta_N2),ncol=2)
    Beta_S = matrix(c(post$beta_S1,post$beta_S2),ncol=2)
    P = post$Pcov

    #selection differential

    #initialize dataframe
    s_SRN = data.frame(s_mu = rep(NA,nrow(Beta_N)), s_psi = rep(NA,nrow(Beta_N)))

    #populate with selection differentials
    for(j in 1:nrow(P)){
      s_SRN[j,] = P[j,,] %*% t(t(Beta_N[j,])) + diag(diag(P[j,,]),2,) %*% Beta_alpha[[j]] %*% t(t(Beta_S[j,])) }
    
    #organize posterior for ease of use
    postl <-
      list(P0V = post$V_P[,1],
           P1V = post$V_P[,2],
           RmV = post$V_R[,1],
           RfV = post$V_R[,2],
           Rcor = post$Rcor[,2,1],
           psi = post$psi_1,
           Beta_psi = Beta_psi,
           bN1 = post$beta_N1,
           bN2 = post$beta_N2,
           bS1 = post$beta_S1,
           bS2 = post$beta_S2,
           bD1 = post$beta_D1,
           bD2 = post$beta_D2,
           s_mu = s_SRN$s_mu,
           s_psi = s_SRN$s_psi,
           phi = post$phi )
    
    #median estimates
    median2p[[i]]<-list(
      P0V_med = median(postl$P0V),
      P1V_med = median(postl$P1V),
      RmV_med = median(postl$RmV),
      RfV_med = median(postl$RfV),
      Rcor_med = median(postl$Rcor),
      psi_med = median(postl$psi),
      Beta_psi_med = median(postl$Beta_psi),
      bN1_med = median(postl$bN1),
      bN2_med = median(postl$bN2),
      bS1_med = median(postl$bS1),
      bS2_med = median(postl$bS2),
      bD1_med = median(postl$bD1),
      bD2_med = median(postl$bD2),
      s_mu_med = median(postl$s_mu),
      s_psi_med = median(postl$s_psi),
      phi_med = median(postl$phi) )
    
    
    #bias
    bias2p[[i]]<-list(
      P0V_bias = (median(postl$P0V) - V*2) / (V*2) ,
      P1V_bias = (median(postl$P1V) - V*2) / (V*2) ,
      RmV_bias = (median(postl$RmV) - res_V) / (res_V),
      RfV_bias = (median(postl$RfV) - res_V) / (res_V),
      Rcor_bias = (median(postl$Rcor) - r_R ) / (r_R),
      psi_bias = (median(postl$psi) - psi_1 ) / (psi_1),
      Beta_psi_bias = (median(postl$Beta_psi) - r_alpha) / (r_alpha),
      bN1_bias = (median(postl$bN1) - beta_n1 ) / (beta_n1 ),
      bN2_bias = (median(postl$bN2) - beta_n2 ) / (beta_n2 ),
      bS1_bias = (median(postl$bS1) - beta_s1 ) / (beta_s1 ),
      bS2_bias = (median(postl$bS2) - beta_s2 ) / (beta_s2 ),
      bD1_bias = (median(postl$bD1) - beta_d1 ) / (beta_d1 ),
      bD2_bias = (median(postl$bD2) - beta_d2 ) / (beta_d2 ),
      s_mu_bias = (median(postl$s_mu) - true_differential[1] ) / (true_differential[1] ),
      s_psi_bias = (median(postl$s_psi) - true_differential[2]) / (true_differential[2] ),
      phi_bias = (median(postl$phi) - phi) / phi )
    
    #Uncertainty
    UC2p[[i]]<-list(
      P0V_uc = mad(postl$P0V)/median(postl$P0V),
      P1V_uc = mad(postl$P1V)/median(postl$P1V),
      RmV_uc = mad(postl$RmV)/median(postl$RmV),
      RfV_uc = mad(postl$RfV)/median(postl$RfV),
      Rcor_uc = mad(postl$Rcor)/median(postl$Rcor),
      psi_uc = mad(postl$psi)/median(postl$psi),
      Beta_psi_uc = mad(postl$Beta_psi)/median(postl$Beta_psi),
      bN1_uc = mad(postl$bN1)/median(postl$bN1),
      bN2_uc = mad(postl$bN2)/median(postl$bN2),
      bS1_uc = mad(postl$bS1)/median(postl$bS1),
      bS2_uc = mad(postl$bS2)/median(postl$bS2),
      bD1_uc = mad(postl$bD1)/median(postl$bD1),
      bD2_uc = mad(postl$bD2)/median(postl$bD2),
      s_mu_uc = mad(postl$s_mu)/median(postl$s_mu),
      s_psi_uc = mad(postl$s_psi)/median(postl$s_psi),
      phi_uc = (mad(postl$phi)/median(phi)) )
    
    #PP
    PP2p[[i]]<- list(
      Rcor_pp = sum(postl$Rcor>0)/length(postl$Rcor),
      psi_pp = sum(postl$psi>0)/length(postl$psi),
      Beta_psi_pp = sum(postl$Beta_psi>0)/length(postl$Beta_psi),
      bN1_pp = sum(postl$bN1>0)/length(postl$bN1),
      bN2_pp = sum(postl$bN2<0)/length(postl$bN2),
      bS1_pp = sum(postl$bS1>0)/length(postl$bS1),
      bS2_pp = sum(postl$bS2<0)/length(postl$bS2),
      bD1_pp = sum(postl$bD1<0)/length(postl$bD1),
      bD2_pp = sum(postl$bD2<0)/length(postl$bD2),
      s_mu_pp = sum(postl$s_mu>0)/length(postl$s_mu),
      s_psi_pp = sum(postl$s_psi<0)/length(postl$s_psi),
      phi_pp = sum(postl$phi>0)/length(postl$phi) )
    
    #Rhat
    summary = summary(SAM_m)$summary
    parnames <- c("alpha_0","nu_0","psi_1", "beta_N1", "beta_N2", "beta_S1", "beta_S2", "beta_D1", "beta_D2", 
                  "V_P[1]","V_P[2]", "V_R[1]", "V_R[2]","Rcor[1,2]", "phi")
    
    Rhat2p[[i]] = summary[parnames,"Rhat"]
    
    saveRDS(time2p, "time2p.RDS")
    saveRDS(median2p, "median2p.RDS")
    saveRDS(bias2p, "bias2p.RDS")
    saveRDS(UC2p, "UC2p.RDS")
    saveRDS(PP2p, "PP2p.RDS")
    saveRDS(Rhat2p, "Rhat2p.RDS")
    
    counter <- counter + 1; 
    print(paste(counter/run*100, "% has been processed"))
  }

#####################################################################
#Data 3

  #load simulated datasets
  data3<-readRDS("data3.RDS")
  
  #initialize lists for results
  time3p <- vector("list", length=run)
  median3p <- vector("list", length=run)
  bias3p <- vector("list", length=run)
  UC3p <-vector("list", length=run)
  PP3p <- vector("list", length=run)
  Rhat3p <- vector("list", length=run)

  #run SAM across datasets (this will take awhile)
  counter = 0 #track progress
  for(i in 1:run){
    
    #single dataset
    df = data3[[i]]
    
    #run model
    start_time <- Sys.time()
    SAM_m <- rstan::sampling(SAM, data=df, init="0", iter=n_iter, warmup=n_warm, seed = i,
                             chains=n_chains, cores=n_chains, control=list(adapt_delta=0.95, max_treedepth=10 ))
    end_time <- Sys.time()
    
    #computation time
    time3p[[i]] = as.numeric(end_time - start_time)
    
    ###############################################################
    #Calculate additional posterior quantities
    ###############################################################
    
    #extract posteriors
    post <- rstan::extract(SAM_m)
    
    #temporary vectors for assortment coefficients
    SRN_PV = post$V_P
    SRN_Psd = post$sd_P
    SRN_PVmean = post$V_P / I_partner #expected variance for mean of partners
    SRN_Psdmean = sqrt(SRN_PVmean) #expected SD for mean of partners
    SRN_focal1 <- post$SRN_P[,,1] #individual intercepts
    SRN_focal2 <- post$SRN_P[,,2] #individual slopes
    SRN_partner1 <- cbind(post$partner_meanm[,,1], post$partner_meanf[,,1])
    SRN_partner2 <- cbind(post$partner_meanm[,,2], post$partner_meanf[,,2])
    
    #scale mean partner variance to variance of single partner
    #see appendix SI for details
    
    SRN_partner1s = SRN_partner1
    for(j in 1:nrow(SRN_partner1))
    {SRN_partner1s[j,] = ( SRN_partner1[j,] / SRN_Psdmean[j,1] ) * SRN_Psd[j,1] }
   
    SRN_partner2s = SRN_partner2
    for(j in 1:nrow(SRN_partner2))
    {SRN_partner2s[j,] = ( SRN_partner2[j,] / SRN_Psdmean[j,2] ) * SRN_Psd[j,2] }
    
    #assortment matrix
    Beta_alpha = list()
    
    #generate matrices across each posterior sample
    for(j in 1:nrow(SRN_focal1))
    {
            Beta_mat = matrix(NA,2,2)
            
            #mu' ~ mu
            Beta_mat[1,1] =  cov(SRN_focal1[j,], SRN_partner1s[j,])/var(SRN_focal1[j,])
            #mu' ~ psi
            Beta_mat[2,1] =  cov(SRN_focal2[j,], SRN_partner1s[j,])/var(SRN_focal2[j,])
            #psi' ~ mu                                                        
            Beta_mat[1,2] =  cov(SRN_focal1[j,], SRN_partner2s[j,])/var(SRN_focal1[j,])
            #psi' ~ psi
            Beta_mat[2,2] =  cov(SRN_focal2[j,], SRN_partner2s[j,])/var(SRN_focal2[j,])
            
            Beta_alpha[[j]] = Beta_mat
    } 
    
    #extract beta_psi'psi (assortment on social plasticity)
    Beta_psi = unlist(lapply(Beta_alpha, function(x) x[2,2]))
    
    #generate other relevant matrices
    Beta_N = matrix(c(post$beta_N1,post$beta_N2),ncol=2)
    Beta_S = matrix(c(post$beta_S1,post$beta_S2),ncol=2)
    P = post$Pcov

    #selection differential

    #initialize dataframe
    s_SRN = data.frame(s_mu = rep(NA,nrow(Beta_N)), s_psi = rep(NA,nrow(Beta_N)))
    
    #populate with selection differentials
    for(j in 1:nrow(P)){
      s_SRN[j,] = P[j,,] %*% t(t(Beta_N[j,])) + diag(diag(P[j,,]),2,) %*% Beta_alpha[[j]] %*% t(t(Beta_S[j,])) }
    
    #organize posterior for ease of use
    postl <-
      list(P0V = post$V_P[,1],
           P1V = post$V_P[,2],
           RmV = post$V_R[,1],
           RfV = post$V_R[,2],
           Rcor = post$Rcor[,2,1],
           psi = post$psi_1,
           Beta_psi = Beta_psi,
           bN1 = post$beta_N1,
           bN2 = post$beta_N2,
           bS1 = post$beta_S1,
           bS2 = post$beta_S2,
           bD1 = post$beta_D1,
           bD2 = post$beta_D2,
           s_mu = s_SRN$s_mu,
           s_psi = s_SRN$s_psi,
           phi = post$phi )
    
    #median estimates
    median3p[[i]]<-list(
      P0V_med = median(postl$P0V),
      P1V_med = median(postl$P1V),
      RmV_med = median(postl$RmV),
      RfV_med = median(postl$RfV),
      Rcor_med = median(postl$Rcor),
      psi_med = median(postl$psi),
      Beta_psi_med = median(postl$Beta_psi),
      bN1_med = median(postl$bN1),
      bN2_med = median(postl$bN2),
      bS1_med = median(postl$bS1),
      bS2_med = median(postl$bS2),
      bD1_med = median(postl$bD1),
      bD2_med = median(postl$bD2),
      s_mu_med = median(postl$s_mu),
      s_psi_med = median(postl$s_psi), 
      phi_med = median(postl$phi))
    
    
    #bias
    bias3p[[i]]<-list(
      P0V_bias = (median(postl$P0V) - V*2) / (V*2) ,
      P1V_bias = (median(postl$P1V) - V*2) / (V*2) ,
      RmV_bias = (median(postl$RmV) - res_V) / (res_V),
      RfV_bias = (median(postl$RfV) - res_V) / (res_V),
      Rcor_bias = (median(postl$Rcor) - r_R ) / (r_R),
      psi_bias = (median(postl$psi) - psi_1 ) / (psi_1),
      Beta_psi_bias = (median(postl$Beta_psi) - r_alpha) / (r_alpha),
      bN1_bias = (median(postl$bN1) - beta_n1 ) / (beta_n1 ),
      bN2_bias = (median(postl$bN2) - beta_n2 ) / (beta_n2 ),
      bS1_bias = (median(postl$bS1) - beta_s1 ) / (beta_s1 ),
      bS2_bias = (median(postl$bS2) - beta_s2 ) / (beta_s2 ),
      bD1_bias = (median(postl$bD1) - beta_d1 ) / (beta_d1 ),
      bD2_bias = (median(postl$bD2) - beta_d2 ) / (beta_d2 ),
      s_mu_bias = (median(postl$s_mu) - true_differential[1] ) / (true_differential[1] ),
      s_psi_bias = (median(postl$s_psi) - true_differential[2]) / (true_differential[2] ),
      phi_bias = (median(postl$phi) - phi) / phi )
    
    #Uncertainty
    UC3p[[i]]<-list(
      P0V_uc = mad(postl$P0V)/median(postl$P0V),
      P1V_uc = mad(postl$P1V)/median(postl$P1V),
      RmV_uc = mad(postl$RmV)/median(postl$RmV),
      RfV_uc = mad(postl$RfV)/median(postl$RfV),
      Rcor_uc = mad(postl$Rcor)/median(postl$Rcor),
      psi_uc = mad(postl$psi)/median(postl$psi),
      Beta_psi_uc = mad(postl$Beta_psi)/median(postl$Beta_psi),
      bN1_uc = mad(postl$bN1)/median(postl$bN1),
      bN2_uc = mad(postl$bN2)/median(postl$bN2),
      bS1_uc = mad(postl$bS1)/median(postl$bS1),
      bS2_uc = mad(postl$bS2)/median(postl$bS2),
      bD1_uc = mad(postl$bD1)/median(postl$bD1),
      bD2_uc = mad(postl$bD2)/median(postl$bD2),
      s_mu_uc = mad(postl$s_mu)/median(postl$s_mu),
      s_psi_uc = mad(postl$s_psi)/median(postl$s_psi),
      phi_uc = (mad(postl$phi)/median(phi)) )
    
    #PP
    PP3p[[i]]<- list(
      Rcor_pp = sum(postl$Rcor>0)/length(postl$Rcor),
      psi_pp = sum(postl$psi>0)/length(postl$psi),
      Beta_psi_pp = sum(postl$Beta_psi>0)/length(postl$Beta_psi),
      bN1_pp = sum(postl$bN1>0)/length(postl$bN1),
      bN2_pp = sum(postl$bN2<0)/length(postl$bN2),
      bS1_pp = sum(postl$bS1>0)/length(postl$bS1),
      bS2_pp = sum(postl$bS2<0)/length(postl$bS2),
      bD1_pp = sum(postl$bD1<0)/length(postl$bD1),
      bD2_pp = sum(postl$bD2<0)/length(postl$bD2),
      s_mu_pp = sum(postl$s_mu>0)/length(postl$s_mu),
      s_psi_pp = sum(postl$s_psi<0)/length(postl$s_psi),
      phi_pp = sum(postl$phi>0)/length(postl$phi) )
    
    #Rhat
    summary = summary(SAM_m)$summary
    parnames <- c("alpha_0","nu_0","psi_1", "beta_N1", "beta_N2", "beta_S1", "beta_S2", "beta_D1", "beta_D2", 
                  "V_P[1]","V_P[2]", "V_R[1]", "V_R[2]","Rcor[1,2]", "phi")
    
    Rhat3p[[i]] = summary[parnames,"Rhat"]
    
    saveRDS(time3p, "time3p.RDS")
    saveRDS(median3p, "median3p.RDS")
    saveRDS(bias3p, "bias3p.RDS")
    saveRDS(UC3p, "UC3p.RDS")
    saveRDS(PP3p, "PP3p.RDS")
    saveRDS(Rhat3p, "Rhat3p.RDS")
    
    counter <- counter + 1; 
    print(paste(counter/run*100, "% has been processed"))
  }

#####################################################################
#Results
#####################################################################

#set wd
setwd()

#load packages
library(bayestestR)
library(plyr)
library(tidyr)

#####################################################################

#Median

#load PP lists
med1 <- readRDS("median1p.RDS")
med2 <- readRDS("median2p.RDS")
med3 <- readRDS("median3p.RDS")

#unlist to dataframe
med1<-ldply (med1, data.frame)
med2<-ldply (med2, data.frame)
med3<-ldply (med3, data.frame)

#combine to single dataframe
med = rbind(med1,med2,med3)

#summarize
apply(med1,2,median); apply(med1,2,quantile, c(0.05,0.95) )
apply(med2,2,median); apply(med2,2,quantile, c(0.05,0.95) )
apply(med3,2,median); apply(med3,2,quantile, c(0.05,0.95) )

#save
saveRDS(med, "medp.RDS")

#####################################################################

#bias 

#load bias lists
bias1 <- readRDS("bias1p.RDS")
bias2 <- readRDS("bias2p.RDS")
bias3 <- readRDS("bias3p.RDS")

#unlist to dataframe
bias1<-ldply (bias1, data.frame)
bias2<-ldply (bias2, data.frame)
bias3<-ldply (bias3, data.frame)

#combine to single dataframe
bias<-rbind(bias1,bias2,bias3)

#summarize
apply(bias1,2,median); apply(bias1,2,hdi, 0.9, allowSplit=FALSE)
apply(bias2,2,median); apply(bias1,2,hdi, 0.9, allowSplit=FALSE)
apply(bias3,2,median); apply(bias1,2,hdi, 0.9, allowSplit=FALSE)

#save
saveRDS(bias, "biasp.RDS")

#save full results
ci = matrix(unlist(apply(bias1,2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
b100 = data.frame( name = colnames(bias1), med = round(apply(bias1,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )

write.csv(b100, "bias100p.csv")

ci = matrix(unlist(apply(bias2,2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
b200 = data.frame( name = colnames(bias2),med = round(apply(bias2,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(b200, "bias200p.csv")

ci = matrix(unlist(apply(bias3,2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
b300 = data.frame( name = colnames(bias3),  med = round(apply(bias3,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(b300, "bias300p.csv")


#####################################################################

#UC

#load PP lists
UC1 <- readRDS("UC1p.RDS")
UC2 <- readRDS("UC2p.RDS")
UC3 <- readRDS("UC3p.RDS")

#unlist to dataframe
UC1<-ldply (UC1, data.frame)
UC2<-ldply (UC2, data.frame)
UC3<-ldply (UC3, data.frame)

#combine to single data frame
UC<-rbind(UC1, UC2, UC3)

#summarize
apply(UC1,2,median); apply(UC1,2,hdi, 0.9, allowSplit=FALSE)
apply(UC2,2,median); apply(UC1,2,hdi, 0.9, allowSplit=FALSE)
apply(UC3,2,median); apply(UC1,2,hdi, 0.9, allowSplit=FALSE)

#save
saveRDS(UC, "UCp.RDS")

#save full results
ci = matrix(unlist(apply(abs(UC1),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
uc100 = data.frame( name = colnames(UC1), med = round(apply(abs(UC1),2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(uc100, "UC100p.csv")

ci = matrix(unlist(apply(abs(UC2),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
uc200 = data.frame( name = colnames(UC2), med = round(apply(abs(UC2),2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(uc200, "UC200p.csv")

ci = matrix(unlist(apply(abs(UC3),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
uc300 = data.frame( name = colnames(UC3), med = round(apply(abs(UC3),2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(uc300, "UC300p.csv")

#####################################################################

#PP

#load PP lists
PP1 <- readRDS("PP1p.RDS")
PP2 <- readRDS("PP2p.RDS")
PP3 <- readRDS("PP3p.RDS")

#unlist to dataframe
PP1<-ldply (PP1, data.frame)
PP2<-ldply (PP2, data.frame)
PP3<-ldply (PP3, data.frame)

#combine to single data frame
PP<-rbind(PP1, PP2, PP3)

#summarize
apply(PP1,2,median); apply(PP1,2,hdi, 0.9, allowSplit=FALSE)
apply(PP2,2,median); apply(PP1,2,hdi, 0.9, allowSplit=FALSE)
apply(PP3,2,median); apply(PP1,2,hdi, 0.9, allowSplit=FALSE)

#save
saveRDS(PP, "PPp.RDS")

#save full organized results
ci = matrix(unlist(apply(abs(PP1),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
PP100 = data.frame( name = colnames(PP1), med = round(apply(PP1,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(PP100, "PP100p.csv")

ci = matrix(unlist(apply(abs(PP2),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
PP200 = data.frame( name = colnames(PP2), med = round(apply(PP2,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(PP200, "PP200p.csv")

ci = matrix(unlist(apply(abs(PP3),2,hdi, 0.9, allowSplit=FALSE)), ncol = 3, byrow =TRUE)
PP300 = data.frame( name = colnames(PP3), med = round(apply(PP3,2,median),2),
            hdil = paste0(round(ci[,2],2),", ", round(ci[,3],2) ) )
write.csv(PP300, "PP300p.csv")

#####################################################################
#Summarize and plot results
#####################################################################

#set wd
setwd()

#load packages
library(ggplot2); library(tidybayes); library(tidyr); library(cowplot)

#load results
med = readRDS("medp.RDS")
bias = readRDS("biasp.RDS")
UC = readRDS("UCp.RDS")
PP = readRDS("PPp.RDS")

#absolute uncertainty
UC = abs(UC)

#wide to long
med.l = gather(med, parameter, median, factor_key=FALSE)
bias.l = gather(bias, parameter, bias, factor_key=FALSE)
UC.l = gather(UC, parameter, uc, factor_key=FALSE)
PP.l = gather(PP, parameter, pp, factor_key=FALSE)

#add sample size
med.l$size = rep( rep(c(100,200,300), each = (nrow(med))/3), ncol(med))
bias.l$size = rep( rep(c(100,200,300), each = (nrow(bias))/3), ncol(bias))
UC.l$size = rep( rep(c(100,200,300), each = (nrow(UC))/3), ncol(UC))
PP.l$size = rep( rep(c(100,200,300), each = (nrow(PP))/3), ncol(PP))

#subset to relevant parameters and create plots
#####################################################################

med.s = med.l[med.l$parameter %in% 
         paste0(c("psi_","Beta_psi_", "s_mu_","s_psi_"), "med") , ]
bias.s = bias.l[bias.l$parameter %in% #derived from med
         paste0(c("psi_","Beta_psi_", "s_mu_","s_psi_"), "bias") , ]
UC.s = UC.l[UC.l$parameter %in% 
         paste0(c("psi_","Beta_psi_", "s_mu_","s_psi_"), "uc") , ]
PP.s = PP.l[PP.l$parameter %in% 
         paste0(c("psi_","Beta_psi_", "s_mu_","s_psi_"), "pp") , ]

#set factor level order
med.s$parameter =factor(med.s$parameter, levels= unique(med.s$parameter) )
bias.s$parameter =factor(bias.s$parameter, levels= unique(bias.s$parameter) )
UC.s$parameter =factor(UC.s$parameter, levels= unique(UC.s$parameter) )
PP.s$parameter =factor(PP.s$parameter, levels= unique(PP.s$parameter) )

#parameter labels
library(latex2exp)

param_label =
c(
  parse(text=TeX('$\\psi_{1}$')),
  parse(text=TeX('$\\beta_{\\bar{\\psi}^{\\prime}\\psi}$')),
 parse(text=TeX('$s_{\\bar{\\mu}}$')),
  parse(text=TeX('$s_{\\bar{\\psi}}$'))
  ) 


levels(med.s$parameter) = param_label
levels(bias.s$parameter) = param_label
levels(UC.s$parameter) = param_label
levels(PP.s$parameter) = param_label

######################################################################
#bias
bias.plotp = 
ggplot(bias.s, aes(x = bias, y = factor(size,levels=c(300,200,100)),
                   group = factor(size,levels=c(300,200,100)), 
                   color = factor(size,levels=c(300,200,100)),
                   fill = factor(size,levels=c(300,200,100)) ) ) +
stat_pointinterval(alpha = 0.45, point_interval = median_hdci,
                    size = 5, .width=c(0.9))+
  scale_fill_manual(values = c("#b91d73","#e378b3","#f953c6"))+
  facet_wrap(.~factor(parameter,levels=(levels(parameter))), ncol=6,
             labeller=label_parsed)+
  coord_cartesian(xlim = c(-1,1))+
  geom_vline(xintercept = c(-0.2,0.2),  linetype="dashed")+ 
    geom_vline(xintercept = 0, linetype=c("solid"))+
  scale_x_continuous(breaks = c(-1,0,1))+
  scale_color_manual(values = c("#b91d73","#e378b3","#f953c6"))+
  xlab("Bias\n")+
  ylab("\n")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        strip.text.x = element_text(size = 14), 
        axis.line = element_line(size = 1),
        panel.spacing = unit(1.1, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.2,0.05,0.2), "cm"))+
        guides(fill=FALSE, color=FALSE)

#save
save_plot("Fig 3i1.tiff", bias.plotp, compression="lzw",
              dpi=600, base_height=2, base_width=7 )

######################################################################
#uncertainty

UC.plotp = 
ggplot(UC.s, aes(x = uc, y = factor(size,levels=c(300,200,100)),
                   group = factor(size,levels=c(300,200,100)), 
                   color = factor(size,levels=c(300,200,100)),
                   fill = factor(size,levels=c(300,200,100))  
                 ) ) +
stat_pointinterval(alpha = 0.45, point_interval = median_hdci,
                    size = 5, .width=c(0.9))+
  scale_color_manual(values = c("#1751ad","#4e9ed9", "#94d2ff"))+
  facet_wrap(.~factor(parameter,levels=(levels(parameter))), ncol=6,
             labeller=label_parsed)+
  coord_cartesian(xlim = c(-0.1,2.1))+ 
  geom_vline(xintercept = c(0.5),  linetype="dashed")+ 
    geom_vline(xintercept = 0, linetype=c("solid"))+
  scale_x_continuous(breaks = c(0,1,2), expand=c(0,0))+
  scale_fill_manual(values = c("#1751ad","#4e9ed9", "#94d2ff"))+
  xlab("Uncertainty\n")+
  ylab("Sample size\n")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        strip.text.x = element_text(size = 14),
        axis.line = element_line(size = 1),
        panel.spacing = unit(1.1, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.2,0.05,0.2), "cm"))+
        guides(fill=FALSE, color=FALSE)

#save
save_plot("Fig 3i2.tiff", UC.plotp, compression="lzw",
              dpi=600, base_height=2, base_width=7 )

######################################################################
#uncertainty

PP.plotp = 
ggplot(PP.s, aes(x = pp, y = factor(size,levels=c(300,200,100)),
                   group = factor(size,levels=c(300,200,100)), 
                   color = factor(size,levels=c(300,200,100)),
                   fill = factor(size,levels=c(300,200,100)) )) +
  stat_pointinterval(alpha = 0.45, point_interval = median_hdci,
                    size = 5, .width=c(0.9))+
  facet_wrap(.~factor(parameter,levels=(levels(parameter))), ncol=6,
             labeller=label_parsed)+
  coord_cartesian(xlim = c(0.55,1.05))+ 
  geom_vline(xintercept = c(0.95),  linetype="dashed")+ 
  geom_vline(xintercept = c(1),  linetype="solid")+ 
  scale_x_continuous(breaks = c(0.6,0.8,1.0), expand=c(0,0))+
  scale_color_manual(values =c("#0c8708","#5ad656", "#aefcac") )+
  xlab("Power\n")+
  ylab("\n")+
  theme(axis.ticks.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=12,face="bold"),
        axis.title.y=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        strip.text.x = element_text(size = 14),
        axis.line = element_line(size = 1),
        panel.spacing = unit(1.1, "lines"),
        panel.border=element_rect(fill=NA,color="black", size=1, 
                                  linetype="solid"),
        panel.background= element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.05,0.2,0.05,0.2), "cm"))+
        guides(fill=FALSE, color=FALSE)

#save
save_plot("Fig 3i3.tiff", PP.plotp, compression="lzw",
              dpi=600, base_height=2, base_width=7 )

######################################################################
#group plots
library(cowplot)

p.gridp <-
  plot_grid(bias.plotp, UC.plotp, PP.plotp, ncol=1, nrow=3,  align="h")

#save
save_plot("Fig 3i.tiff", p.gridp, compression="lzw",
              dpi=600, base_height=5.5, base_width=6.2)


saveRDS(bias.plotp, "bias_plotp.RDS")
saveRDS(UC.plotp, "uc_plotp.RDS")
saveRDS(PP.plotp, "pp_plotp.RDS")

######################################################################
#group with QG results
bias.plot = readRDS("bias_plot.RDS")
UC.plot = readRDS("uc_plot.RDS")
PP.plot = readRDS("pp_plot.RDS")

p.grid <-
  plot_grid(bias.plot, UC.plot, PP.plot, ncol=1, nrow=3,  align="h")

p.comb <-
  plot_grid(NULL, p.gridp, p.grid, nrow = 2 )

#save
save_plot("Fig 3.tiff", p.comb, compression="lzw",
              dpi=600, base_height=10, base_width=9)







