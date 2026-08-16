# =============================================================================
# 1. Data Generation Function
# =============================================================================

#' Generate Spatial Count Data with Train/Test Split
#' 
#' @param n Total sample size
#' @param beta0 Intercept (default = 1)
#' @param sigma2 Spatial variance (default = 1)
#' @param phi Spatial range parameter
#' @param nu Matern smoothness (default = 0.5, exponential)
#' @param tau2 Nugget variance (default = 0)
#' @param response_type One of "poisson", "negbin", "zip"
#' @param kappa Dispersion parameter for negative binomial (default = 2)
#' @param pi_zero Zero-inflation probability for ZIP (default = 0.3)
#' @param train_prop Proportion for training set (default = 0.8)
#' @param seed Random seed for reproducibility
#' 
#' @return List containing:
#'   - train_data: data.frame with x, y, z, count for training
#'   - test_data: data.frame with x, y, z, count for testing
#'   - true_params: list of true parameter values
#'   - x,y for location coordinates, z =  β0 + w(s) + ε(s), count at(x,y) 

generate_spatial_count_data <- function(n = 1000,
                                        beta0 = 1,
                                        sigma2 = 1,
                                        phi = 0.05,
                                        nu = 0.5,
                                        tau2 = 0,
                                        response_type = "poisson",
                                        kappa = 2,
                                        pi_zero = 0.3,
                                        train_prop = 0.8,
                                        seed = NULL)
 {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Generate random training locations in [0,1]^2
  coords <- cbind(x = runif(n), y = runif(n))

  # Form prediction grid
  grid_size <- 40
  grid_seq <- seq(0, 1, length.out = grid_size)
  pred_grid <- expand.grid(x = grid_seq, y = grid_seq)
  coords_round <- round(coords, 5)
  grid_round <- round(as.matrix(pred_grid), 5)
  overlap_idx <- apply(coords_round, 1, function(pt) {
    any(grid_round[,1] == pt[1] & grid_round[,2] == pt[2])
  })
  coords <- coords[!overlap_idx, ]   #remove coords on grid

  coords_all <- rbind(coords, as.matrix(pred_grid))
  n_all <- nrow(coords_all)
  
  # Use geoR to generate Gaussian random field
  # For nu = 0.5 (exponential covariance)
  # geoR uses "matern" with kappa parameter for smoothness
  
  # Generate spatial random effects using geoR::grf
  # Note: geoR's grf() generates data at irregular locations
  grf_result <- geoR::grf(
    n = n,
    grid = "irreg",
    cov.model = "matern",
    cov.pars = c(sigma2, phi),  # c(sill, range)
    kappa = nu,  # smoothness parameter
    nugget = tau2,
    mean = 0,
    method = "cholesky"  # Use Cholesky decomposition
  )
  
  dist_mat <- as.matrix(dist(coords_all))
  
  # Matern correlation function
  # For nu = 0.5 (exponential): C(h) = sigma2 * exp(-h/phi)
  if (nu == 0.5) {
    cov_mat <- sigma2 * exp(-dist_mat / phi)
  } else {
    # General Matern using geoR's cov.spatial function
    cov_mat <- geoR::cov.spatial(
      obj = dist_mat,
      cov.model = "matern",
      cov.pars = c(sigma2, phi),
      kappa = nu
    )
  }
  
  # Add nugget to diagonal
  cov_mat <- cov_mat + diag(tau2, n_all)
  
  # Generate spatial random effects
  w <- MASS::mvrnorm(n = 1, mu = rep(0, n_all), Sigma = cov_mat)
  
  # Linear predictor
  z <- beta0 + w
  
  # Generate count responses based on type
  if (response_type == "poisson") {
    lambda <- exp(z)
    y <- rpois(n_all, lambda)
    
  } else if (response_type == "negbin") {
    mu <- exp(z)
    # Parameterization: Var = mu + mu^2/kappa
    y <- rnbinom(n_all, size = kappa, mu = mu)
    
  } else if (response_type == "zip") {
    lambda <- exp(z)
    # Zero-inflated: with prob (1-pi_zero) get 0, otherwise Poisson
    u <- runif(n_all)
    y <- ifelse(u < (1 - pi_zero), 0, rpois(n_all, lambda))
    
  } else {
    stop("response_type must be one of: poisson, negbin, zip")
  }
  
  # Create data frame
  data_full <- data.frame(
    x = coords_all[, 1],
    y = coords_all[, 2],
    z = z,
    count = y
  )
  
  n_train <- nrow(coords)
  train_data <- data_full[1:n_train, ]
  test_data  <- data_full[(n_train + 1):n_all, ]
  
  # Save true parameters
  true_params <- list(
    beta0 = beta0,
    sigma2 = sigma2,
    phi = phi,
    nu = nu,
    tau2 = tau2,
    response_type = response_type,
    kappa = kappa,
    pi_zero = pi_zero
  )
  
  return(list(
    train_data = train_data,
    test_data = test_data,
    true_params = true_params
  ))
}


# =============================================================================
# 2. spaMM Fitting and Prediction Function
# =============================================================================

#' Fit and predict using spaMM (PQL approximation)
#' 
#' @param train_data Training data frame
#' @param test_data Test data frame
#' @param family Family for GLM (e.g., "poisson", "negbin", "zip")
#' 
#' @return List with predictions, metrics, and timing

fit_predict_spaMM <- function(train_data, test_data, family ) 
{
  
  # Start timing
  time_start <- Sys.time()
  
  if (family == "poisson") {
    family_spec <- poisson()
  } else if (family == "zip") {
    family_spec <- spaMM::negbin1()
  } else if (family == "negbin") {
    family_spec <- spaMM::negbin2()
  } else {
    stop("family must 'poisson', 'nbinom1', or 'nbinom2'")
  }  
  
  # Fit model
  if(nrow(train_data)>1000)
  {
		fit <- tryCatch({
		  spaMM::fitme(
			count ~ 1 + Matern(1 | x + y),
			data = train_data,
			family = family_spec,
			fixed = list(nu = 0.5),
			control.HLfit=list(algebra="decorr") 
		  )
		}, error = function(e) {
		  message("spaMM fitting error: ", e$message)
		  return(NULL)
		})		
  } else
  {
		fit <- tryCatch({
		  spaMM::fitme(
			count ~ 1 + Matern(1 | x + y),
			data = train_data,
			family = family_spec,
			fixed = list(nu = 0.5)
		  )
		}, error = function(e) {
		  message("spaMM fitting error: ", e$message)
		  return(NULL)
		})
  }
  
  time_fit <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  if (is.null(fit)) {
    return(list(
      method = "spaMM",
      predictions = NA,
      rmspe = NA,
      mae = NA,
      time_fit = NA,
      time_pred = NA,
      time_total = NA,
      error = TRUE
    ))
  }
  
  # Prediction
  time_pred_start <- Sys.time()
  
  pred <- tryCatch({
    predict(fit, newdata = test_data, type = "response",
	variances=list(respVar=TRUE))
  }, error = function(e) {
    message("spaMM prediction error: ", e$message)
    return(NULL)
  })
  
  time_pred <- as.numeric(difftime(Sys.time(), time_pred_start, units = "secs"))
  time_total <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  if (is.null(pred)) {
    return(list(
      method = "spaMM",
      predictions = NA,
      rmspe = NA,
      mae = NA,
      time_fit = time_fit,
      time_pred = NA,
      time_total = time_total,
      error = TRUE
    ))
  }
  
  # Calculate PMSE
  pmse <- mean((pred - test_data$count)^2)
  
  # Calculate MAE
  mae <- mean(abs(pred - test_data$count))
  
  return(list(
    method = "spaMM",
    predictions = pred,
	  estpvar=attr(pred,"predVar"),   
    rmspe = sqrt(pmse),
    mae = mae,
    time_fit = time_fit,
    time_pred = time_pred,
    time_total = time_total,
    error = FALSE
  ))
}


# =============================================================================
# 3. INLA Fitting and Prediction Function
# =============================================================================

#' Fit and predict using INLA (SPDE approach)
#' 
#' @param train_data Training data frame
#' @param test_data Test data frame
#' @param family_name  Family for GLM (e.g., "poisson", "negbin", "zip")
#' 
#' @return List with predictions, metrics, and timing


fit_predict_bru_weak_prior <- function(train_data, test_data, family_name = "poisson") 
{

  time_start <- Sys.time()

  safe_return <- function(time_fit = NA, time_pred = NA, time_total = NA) {
    return(list(
      method = "inlabru",
      predictions = NA,
      estpvar = NA,
      rmspe = NA,
      mae = NA,
      time_fit = time_fit,
      time_pred = time_pred,
      time_total = time_total,
      inla_family = NA,
      error = TRUE
    ))
  }
  
  family_name <- tolower(family_name)
  if (family_name == "poisson") {
    inla_family <- "poisson"
  } else if (family_name == "negbin") {
    inla_family <- "nbinomial"
  } else if (family_name == "zip") {
    inla_family <- "zeroinflatedpoisson2  "
  } else {
    stop("family_name must be one of: poisson, negbin, zip")
  }
  
  
  # ---------------------------
  # Convert to sf
  # ---------------------------
  
  train_sf <- tryCatch(
    st_as_sf(train_data, coords = c("x","y"), crs = NA),
    error = function(e) return(NULL)
  )
  
  test_sf <- tryCatch(
    st_as_sf(test_data, coords = c("x","y"), crs = NA),
    error = function(e) return(NULL)
  )
  
  if (is.null(train_sf) || is.null(test_sf)) {
    return(safe_return())
  }  
  
  coords <- as.matrix(train_data[,c("x","y")])
  
  # ---------------------------
  # Mesh 
  # ---------------------------
  
  mesh <- tryCatch(
    fm_mesh_2d(loc = coords),
    error = function(e) return(NULL)
  )
  
  if (is.null(mesh)) {
    return(safe_return())
  }
  
  
  # ---------------------------
  # weak prior
  # ---------------------------
  
  spde <- tryCatch(
    inla.spde2.pcmatern(
      mesh = mesh,
      alpha = 2,
      prior.range = c(0.01, 0.5),
      prior.sigma = c(10, 0.05)
    ),
    error = function(e) return(NULL)
  )
  
  if (is.null(spde)) {
    return(safe_return())
  }
 
  cmp <- ~ intercept(1) +
    spatial(geometry, model = spde)
  
  # ---------------------------
  # Fit model
  # ---------------------------
  
  time_fit_start <- Sys.time()
  
  fit_final <- tryCatch({    
    bru(
      components = cmp,
      formula = count ~ intercept + spatial,
      data = train_sf,
      family = inla_family
    )    
  }, error = function(e) {    
    return(NULL)    
  })
  
  
  time_fit <- as.numeric(
    difftime(
      Sys.time(),
      time_fit_start,
      units = "secs"
    )
  )
  
  
  if (is.null(fit_final)) {
    return(safe_return(time_fit = time_fit))
  }
  
  
  # ---------------------------
  # Prediction
  # ---------------------------
  
  time_pred_start <- Sys.time()  
  
  pred <- tryCatch({    
    predict(
      fit_final,
      newdata = test_sf,
      formula = ~ exp(intercept + spatial)
    )    
  }, error = function(e) return(NULL))  
  
  time_pred <- as.numeric(
    difftime(
      Sys.time(),
      time_pred_start,
      units = "secs"
    )
  )  
  
  time_total <- as.numeric(
    difftime(
      Sys.time(),
      time_start,
      units = "secs"
    )
  )
    
  if (is.null(pred)) {
    return(safe_return(
      time_fit = time_fit,
      time_pred = time_pred,
      time_total = time_total
    ))
  }
  
  
  # ---------------------------
  # Extract predictions
  # ---------------------------
  
  pred_mean <- pred$mean
  pred_var  <- pred$sd^2  
  
  mpse <- tryCatch(
    mean((pred_mean - test_data$count)^2),
    error = function(e) NA
  )  
  
  mae <- tryCatch(
    mean(abs(pred_mean - test_data$count)),
    error = function(e) NA
  )  
  
  if (length(pred_mean) > 1) {
    pred_mean <- as.matrix(pred_mean)
  }
  
  
  return(list(
    method = "inlabru",
    predictions = pred_mean,
    estpvar = pred_var,
    rmspe = sqrt(mpse),
    mae = mae,
    time_fit = time_fit,
    time_pred = time_pred,
    time_total = time_total,
    inla_family = inla_family,
    error = FALSE
  ))
}


# =============================================================================
# 4. ensGP Fitting and Prediction Function
# =============================================================================

#' Fit and predict using ensemble GP
#' 
#' @param train_data Training data frame
#' @param test_data Test data frame
#' @param n_base Number of base learners (default = 10)
#' 
#' @return List with predictions, metrics, and timing


fit_predict_ensGP <- function(train_data, test_data, n_base = 10) 
{
  
  # Start timing
  time_start <- Sys.time()
  
  # Create a grid of candidate range parameters
  phi_grid <- seq(0.03, 0.3, length.out = n_base)
  
  # Prepare coordinates
  coords_train <- as.matrix(train_data[, c("x", "y")])
  coords_test <- as.matrix(test_data[, c("x", "y")])
  
  # Storage for base learner predictions
  base_predictions <- matrix(NA, nrow = nrow(test_data), ncol = n_base)
  
  # Fit base learners with different fixed phi values
  for (i in 1:n_base) {
    
    base_fit <- tryCatch({
      # Use a simple conjugate model (e.g., Gaussian approximation)
      # In practice, would use more sophisticated base learners
      
      # Compute distance matrices
      dist_train <- as.matrix(dist(coords_train))
      dist_test_train <- as.matrix(fields::rdist(coords_test, coords_train))
      
      # Exponential covariance
      K_train <- exp(-dist_train / phi_grid[i])
      K_test_train <- exp(-dist_test_train / phi_grid[i])
      
      # Simple prediction (Gaussian process prediction)
      # Then transform for count data
      K_inv <- solve(K_train + diag(0.01, nrow(K_train)))
      
      # Use log-transformed counts for Gaussian approximation
      y_trans <- log(train_data$count + 0.5)
      
      pred_latent <- K_test_train %*% K_inv %*% y_trans
      pred_count <- exp(pred_latent)
      
      pred_count
      
    }, error = function(e) {
      rep(mean(train_data$count), nrow(test_data))
    })
    
    base_predictions[, i] <- base_fit
  }
  
  time_fit <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  # Stacking: use equal weights (in practice, would optimize on validation set)
  time_pred_start <- Sys.time()
  
  weights <- rep(1/n_base, n_base)
  pred_mean <- base_predictions %*% weights
  
  time_pred <- as.numeric(difftime(Sys.time(), time_pred_start, units = "secs"))
  time_total <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  # Calculate PMSE
  pmse <- mean((pred_mean - test_data$count)^2)
  
  # Calculate MAE
  mae <- mean(abs(pred_mean - test_data$count))
  
  return(list(
    method = "ensGP",
    predictions = pred_mean,
    rmspe = sqrt(pmse),
    mae = mae,
    time_fit = time_fit,
    time_pred = time_pred,
    time_total = time_total,
    error = FALSE
  ))
}


# =============================================================================
# 5. spStack Fitting and Prediction Function
# =============================================================================

#' Fit and predict using Bayesian Predictive Stacking
#' 
#' @param train_data Training data frame
#' @param test_data Test data frame
#' @param n_base Number of base learners (default = 10)
#' @param family Family for GLM
#' 
#' @return List with predictions, metrics, and timing

fit_predict_spStack <- function(train_data, test_data, n_base = 10) 
{

  ErrObj <- list(
      method = "spStack",
      predictions = NA,
      rmspe = NA,
      mae = NA,
      time_fit = NA,
      time_pred = NA,
      time_total = NA,
      error = TRUE
  )
  
  time_start <- Sys.time()
  
  coords_train <- as.matrix(train_data[, c("x", "y")])
  coords_test  <- as.matrix(test_data[, c("x", "y")])
  
  mod <- tryCatch({
	  spGLMstack(
      formula       = count ~ 1,
      data          = train_data,
      family        = "poisson",
      coords        = coords_train,
      cor.fn        = "matern",           
      params.list   = list(
		    phi= seq(0.03, 0.3, length.out = n_base),   
		    nu = 0.5,                                 
        boundary = 0.5 ),                    
      n.samples     = 1000,                
      loopd.controls = list(method = "CV", CV.K = 10, nMC = 500),
      priors        = list(nu.beta = 5, nu.z = 5),
      parallel      = TRUE,              
      verbose       = FALSE
    )
  }, error = function(e) {
    message("spStack failed", e$message)
    NULL
  })
  
  if(is.null(mod)) return(ErrObj) else
	mod <- recoverGLMscale(mod)  
  
  time_fit <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))  
  
  time_pred_start <- Sys.time()
  
  predictions <- tryCatch({
    if (is.null(mod)) stop("spStack failed")
    
    X_new <- matrix(1, nrow = nrow(test_data), ncol = 1)
    
    mod_pred <- spStack::posteriorPredict(
      mod,                                   
      coords_new = coords_test,
      covars_new = X_new,
      joint      = FALSE                    
    )
    
    
    post_samps <- spStack::stackedSampler(mod_pred)   
    rowMeans(post_samps$y.pred)
  }, error = function(e) {
    message("spStack failed", e$message)
    rep(mean(train_data$count, na.rm = TRUE), nrow(test_data))
  })
  
  time_pred <- as.numeric(difftime(Sys.time(), time_pred_start, units = "secs"))
  time_total <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  if(is.null(predictions)) return(ErrObj)
  
  pmse <- mean((predictions - test_data$count)^2)
  mae  <- mean(abs(predictions - test_data$count))
  

  list(
    method      = "spStack",
    predictions = predictions,
    rmspe        = sqrt(pmse),
    mae         = mae,
    time_fit    = time_fit,
    time_pred   = time_pred,
    time_total  = time_total,
    error       = FALSE
  )
}


#=============================================================================
# 6. glmmTMB Fitting and Prediction Function
# =============================================================================

#' @param train Training data frame with x, y, count columns
#' @param test Test data frame with x, y, count columns
#' @param family Family for count data (default = "poisson")    #can be "poisson", "negbin", "zip"
#' @param spatial_cov Spatial covariance type (default = "mat")
#' 
#' @return List with predictions, metrics, and timing
#' 


run_glmmTMB <- function(train, test, family = "poisson", spatial_cov = "exp", slot=100) 
{
  

  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop(" glmmTMB required: install.packages('glmmTMB')")
  }  

  time_start <- Sys.time()
  

  #train$pos <- numFactor(train$x, train$y)
  #test$pos <- numFactor(test$x, test$y)
  n_train <- nrow(train)
  all_data <- rbind( train, test)
  all_data$pos <- numFactor( all_data$x, all_data$y)
  train <- all_data[1:n_train, ]
  test  <- all_data[(n_train+1):nrow(all_data), ]    
  
  train$group <- factor(1)
  test$group <- factor(1)
  
  
  if (family == "poisson") {
    family_spec <- poisson()
	zifo <- as.formula("~0")
  } else if (family == "zip") {
    family_spec <- poisson()
	zifo <- as.formula("~1")
  } else if (family == "negbin") {
    family_spec <- nbinom2()
	zifo <- as.formula("~0")
  } else {
    stop("family must 'poisson', 'nbinom1', or 'nbinom2'")
  }
  
  if (spatial_cov == "exp") {
    formula_str <- "count ~ 1 + exp(pos + 0 | group)"
  } else if (spatial_cov == "gau") {
    formula_str <- "count ~ 1 + gau(pos + 0 | group)"
  } else if (spatial_cov == "mat") {
    formula_str <- "count ~ 1 + mat(pos + 0 | group)"
  } else {
    stop("spatial_cov must be 'exp', 'gau', or	 'mat'")
  }
  
  if(family == "negbin")
  fit <- tryCatch({
    glmmTMB(
      formula = as.formula(formula_str),
      data = train,
      family = family_spec,
	  ziformula =zifo,
	  control=glmmTMBControl(
	  parallel=list(n=3L, autopar=TRUE)),
      REML = TRUE       
    )
  }, error = function(e) {
    message("glmmTMB error: ", e$message)
    return(NULL)
  })  else
  fit <- tryCatch({
    glmmTMB(
      formula = as.formula(formula_str),
      data = train,
      family = family_spec,
      dispformula = ~0,
	  ziformula =zifo,
	  control=glmmTMBControl(
	  parallel=list(n=3L, autopar=TRUE)),	  
      REML = TRUE       
    )
  }, error = function(e) {
    message("glmmTMB error: ", e$message)
    return(NULL)
  })
  
  time_fit <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  

  if (!is.null(fit)) {
    if (!is.null(fit$sdr) && !fit$sdr$pdHess) {
      warning("Hessian not PSD")
    }
  }
  
  if (is.null(fit)) {
    return(list(
      method = "glmmTMB",
      predictions = NA,
      rmspe = NA,
      mae = NA,
      time_fit = time_fit,
      time_pred = NA,
      time_total = time_fit,
      error = TRUE,
      error_message = "failed"
    ))
  }
  
  # prediction
  time_pred_start <- Sys.time()
  
  pred <- tryCatch({
    predict(fit, newdata = test, type = "response", allow.new.levels = TRUE)
  }, error = function(e) {
    message("glmmTMB prediction error: ", e$message)
    return(NULL)
  })  
  
  time_pred <- as.numeric(difftime(Sys.time(), time_pred_start, units = "secs"))
  time_total <- as.numeric(difftime(Sys.time(), time_start, units = "secs"))
  
  if (is.null(pred)) {
    return(list(
      method = "glmmTMB",
      predictions = NA,
      rmspe = NA,
      mae = NA,
      time_fit = time_fit,
      time_pred = time_pred,
      time_total = time_total,
      error = TRUE,
      error_message = "prediction error"
    ))
  }
  
  
  pred_mean <- as.numeric(pred)
  
  pmse <- mean((pred_mean - test$count)^2)
  mae <- mean(abs(pred_mean - test$count))
  
  return(list(
    method = "glmmTMB",
    predictions = pred_mean,
    rmspe = sqrt(pmse),
    mae = mae,
    time_fit = time_fit,
    time_pred = time_pred,
    time_total = time_total,
    error = FALSE,
    converged = if(!is.null(fit$sdr)) fit$sdr$pdHess else NA
  ))
}

fit_predict_glmmTMB <- function(train_data, test_data, rtype){

  tryCatch(

    callr::r(

      function(fun, train_data, test_data, rtype){

        library(glmmTMB)
        library(TMB)
        library(RcppEigen)
        library(sf)
        library(nlme)

        fun(train_data, test_data, rtype)

      },

      args = list(
        fun = run_glmmTMB,
        train_data = train_data,
        test_data = test_data,
        rtype = rtype
      )

    ),

    error = function(e){

      message(conditionMessage(e))

      list(
        method = "glmmTMB",
        rmpse = NA,
        mae = NA
      )

    }

  )

}