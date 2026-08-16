# Example usage for single run
# Make sure that func.R is in the working directory.

source("func.R")

# set simulation parameters
n=32              #32, 256, 2048
phi=0.05          #0.05, 0.20
tau2=0            #0, 0.25
rtype="poisson"   #"poisson", "negbin", "zip"
seed=234          #use different seed for independent replicates  

# generate data
data_list <- generate_spatial_count_data(
    n = n,
    phi = phi,
    tau2 = tau2,
    response_type = rtype,
    seed = seed  )
  
train_data <- data_list$train_data
test_data <- data_list$test_data

# fit models
res_spaMM <- fit_predict_spaMM(train_data, test_data, rtype ) 
res_inla <- fit_predict_bru_weak_prior(train_data, test_data, rtype)
res_ensGP <- fit_predict_ensGP(train_data, test_data) 
res_spStack <- fit_predict_spStack(train_data, test_data) 
res_TMB <- fit_predict_glmmTMB(train_data, test_data, rtype) 

# summarize results 

metrics=data.frame( method=c("spaMM","INLA","glmmTMB", "ensGP","spStack"), 
        RMSPE=NA, MAE=NA, time_fit=NA, time_pred=NA)
metrics[metrics[,"method"]=="spaMM",2:5]=with(res_spaMM,c(rmspe,mae,time_fit,time_pred))
metrics[metrics[,"method"]=="INLA",2:5]=with(res_inla,c(rmspe,mae,time_fit,time_pred))
metrics[metrics[,"method"]=="glmmTMB",2:5]=with(res_TMB,c(rmspe,mae,time_fit,time_pred))
metrics[metrics[,"method"]=="ensGP",2:5]=with(res_ensGP,c(rmspe,mae,time_fit,time_pred))
metrics[metrics[,"method"]=="spStack",2:5]=with(res_spStack,c(rmspe,mae,time_fit,time_pred))

print(metrics)
