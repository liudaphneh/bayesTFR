
#########################################################
# MCMC sampling for DLpar for UN estimates
#########################################################
tfr.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10, uncertainty=FALSE) {
  if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
    nr_simu <- mcmc$iter
    nr_countries <- mcmc$meta$nr_countries_estimation
    nr_countries_all <- mcmc$meta$nr_countries
    has_extra_countries <- nr_countries != nr_countries_all
    Triangle_c4.low <- mcmc$meta$Triangle_c4.low
    Triangle_c4.up <- mcmc$meta$Triangle_c4.up
    lower_d <- mcmc$meta$d.low
    upper_d <- mcmc$meta$d.up

    id_DL <- mcmc$meta$id_DL[mcmc$meta$id_DL <= nr_countries]
    id_DL_all <- mcmc$meta$id_DL
    id_early <- mcmc$meta$id_early[mcmc$meta$id_early <= nr_countries]
    id_early_all <- mcmc$meta$id_early
    id_notearly <- mcmc$meta$id_notearly[mcmc$meta$id_notearly <= nr_countries]
    id_notearly_all <- mcmc$meta$id_notearly
    nr_notearly <- length(id_notearly)
    nr_DL <- length(id_DL)
    sqrt.nr.dl.plus1 <- sqrt(nr_DL + 1)
    CoverCplus1 <- nr_DL/(nr_DL+1)

    chi0 <- mcmc$meta$chi0
    psi0 <- mcmc$meta$psi0
    nu_psi0 <- mcmc$meta$nu.psi0
    nu_psiD <- nu_psi0 + nr_DL
    prod_psi0 <- nu_psi0*psi0^2
    
    mean_eps_tau_0 <- mcmc$meta$mean.eps.tau0
    sd_eps_tau_0 <- mcmc$meta$sd.eps.tau0
    nu_tau0 <- mcmc$meta$nu.tau0
    nu_tauD <- nu_tau0 + nr_notearly
    prod_tau0 <- nu_tau0*sd_eps_tau_0^2

    delta4_0 <- mcmc$meta$delta4.0
    Triangle4_0  <- mcmc$meta$Triangle4.0
    nu4_D <- mcmc$meta$nu4 +  nr_DL
    prod_delta4_0 <- mcmc$meta$nu4*delta4_0^2
    
    nu_deltaD <- mcmc$meta$nu.delta0 +  nr_DL
    a_delta <- nu_deltaD/2
    prod_delta0 <- mcmc$meta$nu.delta0*mcmc$meta$delta0^2 
    
    suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
    # Daphne edit 20240414 to force suppl.T to be 0 if annual = FALSE 
    # since in my five-year run the suppl.data just duplicates my data...
    if(!mcmc$meta$annual.simulation){
      suppl.T <- 0
    }
    
    # Daphne: load in MCMC object from first stage of estimation
    m.default <- get.tfr.mcmc(mcmc$meta$first.stage.directory, chain.ids = mcmc$meta$first.stage.chains)
    if (mcmc$meta$second.stage.uncertainty){
      m.default.3 <- get.tfr3.mcmc(mcmc$meta$first.stage.directory, chain.ids = mcmc$meta$first.stage.3.chains)
    }
    
    if (uncertainty)
    {
      meta <- mcmc$meta
      nr.countries <- meta$nr.countries
      countries.index <- meta$id_phase3
      ardata <- list()
      Ts <- rep(0, nr.countries)
      for(country in 1:nr.countries) {
        data <- get.observed.tfr(countries.index[country], meta, 'tfr_matrix_all')
        ardata[[country]] <- data[meta$lambda_c[countries.index[country]]:meta$T_end_c[countries.index[country]]]
        Ts[country] <- length(ardata[[country]])
      }
      mcmc$observations <- ardata
      mcmc$length_obs <- Ts
      mcmc$recompute.mu.integral <- TRUE
      mcmc$recompute.rho.integral <- TRUE
    }

    mcenv <- as.environment(mcmc) # Create an environment for the mcmc stuff in order to avoid 
					   			  # copying of the mcmc list 

    mcenv$thin <- thin
    ones <- matrix(1, ncol=nr_DL, nrow=3)
    
  matrix.name <- ifelse(uncertainty, 'tfr_all', 'tfr_matrix')
  
  if(is.null(mcenv$eps_Tc)) 
  {
    if (!is.null(mcenv$meta$ar.phase2) && mcenv$meta$ar.phase2) 
      mcenv$eps_Tc <- get_eps_T_all(mcmc, matrix.name=matrix.name, rho.phase2=mcenv$rho.phase2)
    else mcenv$eps_Tc <- get_eps_T_all(mcmc, matrix.name=matrix.name)
  }
    
	if(has_extra_countries) {
		if (verbose)
			cat('\nCountries', mcmc$meta$regions$country_code[(nr_countries+1):nr_countries_all], 
							'not included in the estimation of world parameters.\n')
	}
    if(length(id_DL_all) < nr_countries_all && verbose) 
        cat('\nCountries', mcmc$meta$regions$country_code[!seq_len(nr_countries_all) %in% id_DL_all],
            'not included in the estimation because in Phase I.\n')
  
    # sd_Tc with sigma0 and sd_tau_eps
    # the non-constant variance is sum of sigma0 and add_to_sd_Tc
    # matrix with each column one country
    mcenv$add_to_sd_Tc <- matrix(NA, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
    mcenv$data.list <- list()
    for (country in 1:nr_countries_all){
    	# could exclude 1:(tau_c-1) here
      this.data <- array(dim = mcenv$meta$T_end_c[country] - 1)
      if (country %in% id_DL_all)
      	this.data[mcenv$meta$start_c[country]:(mcenv$meta$lambda_c[country] - 1)] <- 
      	  get.observed.tfr(country, mcenv$meta, matrix.name=matrix.name)[mcenv$meta$start_c[country]:(mcenv$meta$lambda_c[country] - 1)]
      	
    	mcenv$data.list[[country]] <- this.data
    	mcenv$add_to_sd_Tc[1:(mcenv$meta$T_end_c[country]-1),country] <- (
        			this.data - mcenv$S_sd)*ifelse(this.data > mcenv$S_sd, -mcenv$a_sd, mcenv$b_sd)
    }
    
	mcenv$const_sd_dummie_Tc <- matrix(0, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(mcenv$meta$suppl.data$tfr_matrix) else c(), rownames(mcenv$meta$tfr_matrix)))
	mcenv$const_sd_dummie_Tc[mid.years[-length(mid.years)] < 1975,] <- 1

  	mcenv$mean_eps_Tc <- matrix(0, mcenv$meta$T_end -1+suppl.T, nr_countries_all)
  	
  	# DAPHNE: modified idx.tau_c.id.notearly.all and idx.tau_c.id.notearly
  	#    - tau_c restricts to Phase II country-time pairs
  	#    - but also need to restrict to where we have covariate data (e.g. start_cov_data_c)
  	start.both.phaseII.cov <-	apply(matrix(c(mcenv$meta$start_cov_data_c[id_notearly_all], mcenv$meta$tau_c[id_notearly_all]), ncol = 2), 1, max)
  	idx.tau_c.id.notearly.all <- matrix(c(start.both.phaseII.cov, id_notearly_all), ncol=2)
  	# some of these are NA, so remove from idx.tau_c.id.notearly.all 
  	idx.tau_c.id.notearly.all <- idx.tau_c.id.notearly.all[which(!is.na(idx.tau_c.id.notearly.all[,1])), ]
	mcenv$mean_eps_Tc[idx.tau_c.id.notearly.all] <- mcenv$mean_eps_tau
	
	start.both.phaseII.cov <-	apply(matrix(c(mcenv$meta$start_cov_data_c[id_notearly], mcenv$meta$tau_c[id_notearly]), ncol = 2), 1, max)
	idx.tau_c.id.notearly <- matrix(c(start.both.phaseII.cov, id_notearly), ncol=2)
	# remove NAs
	idx.tau_c.id.notearly <- idx.tau_c.id.notearly[which(!is.na(idx.tau_c.id.notearly[,1])), ]
	
	############ DAPHNE: setting up sampling for covariate term ############

	#### priors for beta
	X.educ <- as.vector(t(mcenv$meta$educ))
	X.fp <- as.vector(t(mcenv$meta$fp))
	idxX.educ <- !is.na(X.educ)  # NA when data is not available
	idxX.fp <- !is.na(X.fp)      # NA when data is not available
	X.g <- as.vector(t(mcenv$meta$gdp))
	idxX.g <- !is.na(X.g)        # NA when data is not available
	SSA.ind <- rep(mcenv$meta$SSA_indicator, times = nrow(mcenv$meta$educ))
	tfr.decr.vector <- as.vector(t(mcenv$meta$tfr_decr))
	idxtfr.decr <- !is.na(tfr.decr.vector) 
	
	# retain non-NA values
	X.educ <- X.educ[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	X.fp <- X.fp[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	X.g <- X.g[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	SSA.ind <- SSA.ind[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	tfr.decr.vector <- tfr.decr.vector[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	
	# prior sd for beta: 0.5*sd(TFR decr)/sd(variable of interest)
	# for sensitivity analysis: wider = *sqrt(2) and narrower = *(1/sqrt(2))
	sdbeta.e <- 0.5*sd(tfr.decr.vector)/sd(X.educ)
	sdbeta.fp <- 0.5*sd(tfr.decr.vector)/sd(X.fp)
	sdbeta.g <- 0.5*sd(tfr.decr.vector)/sd(X.g)
	sdbeta.e.SSA <- 0.5*sd(tfr.decr.vector)/sd(SSA.ind * X.educ)
	sdbeta.fp.SSA <- 0.5*sd(tfr.decr.vector)/sd(SSA.ind * X.fp)
	sdbeta.g.SSA <- 0.5*sd(tfr.decr.vector)/sd(SSA.ind * X.g)
	
	# for beta posterior
	sdbeta.inv.e <- 1/sdbeta.e^2
	sdbeta.inv.fp <- 1/sdbeta.fp^2
	sdbeta.inv.g <- 1/sdbeta.g^2
	sdbeta.inv.e.SSA <- 1/sdbeta.e.SSA^2
	sdbeta.inv.fp.SSA <- 1/sdbeta.fp.SSA^2
	sdbeta.inv.g.SSA <- 1/sdbeta.g.SSA^2

	#### for estimation of between-country correlation (bc_rho)
	# retain non-NA values for clusterUNregion
	cluster <- mcenv$meta$cluster
	cluster <- cluster[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
	
	# find indices when cluster membership changes
	# used to create block diagonal correlation matrix R
	new.clust <- c(1)
	for(i in 1:(length(cluster)-1)){
	  if(cluster[i] != cluster[i+1]){ new.clust <- c(new.clust, i+1) }
	}
	
	# function to obtain n x n block diagonal matrix R^{-1} from given rho
	# used to find posterior of beta
	get.R.from.rho <- function(rho, n, new.clust){
	  R <- diag(n)
	  
	  for(i in 1:(length(new.clust) - 1)){
	    # each of these indices indicates the start of a new block A_j^{-1}
	    # fill in the off-diagonals of that n_j x n_j block with rho
	    
	    # size of block Aj
	    nj <- new.clust[i+1] - new.clust[i]
	    
	    # indices of R corresponding to block A_j^{-1}
	    A.j.idx <- new.clust[i]:(new.clust[i] + nj - 1)
	    
	    # polynomials for entries of A_j^{-1}
	    a.2 <- (-1*rho)/((1-rho) * (1 - rho + nj*rho))
	    a.1 <- 1/(1-rho) + a.2
	    
	    # R[j,k] elements filled in for A_j^{-1}
	    for(j in A.j.idx){
	      for(k in A.j.idx){
	        if(j == k){ R[j,k] <- a.1}  # diagonal
	        if(j != k){ R[j,k] <- a.2 } # off-diagonals
	      }
	    }
	  }
	  return(R)
	}
	
	# function to obtain posterior loglikelihood of rho, conditional on other parameters 
	rho.post.ll <- function(rho, X, Y, beta, sigma, new.clust){
	  # to store loglik
	  ll <- 0
	  
	  # create z = (y- xbeta)/sigma
	  z.vector <- (Y - X %*% beta)/sigma
	  
	  # for each block, obtain contribution to posterior loglik
	  for(j in 1:(length(new.clust) - 1)){
	    
	    # size of block Aj
	    nj <- new.clust[j+1] - new.clust[j]
	    
	    # values of z.vector corresponding to block A_j
	    zj <- z.vector[j:(j + nj)]
	    
	    # add this block's contribution to loglik
	    ll <- ll + rho.post.ll.Aj(rho, nj, zj)
	  }
	  
	  # add indicator function
	  ll <- ll + log((rho >= 0) & (rho <= 1))
	  return(ll)
	}
	
	# function to obtain component for block A_j for posterior log likelihood of rho
	#      nj = size of block A_j
	#      zj = z.vector for block A_j
	rho.post.ll.Aj <- function(rho, nj, zj){
	  # define polynomials that frequently show up
	  a.2 <- (-1*rho)/((1-rho) * (1 - rho + nj*rho))
	  a.1 <- 1/(1-rho) + a.2
	  
	  # determinant part
	  det.Aj <- -(1/2) * (log(1 - rho + nj*rho) + (nj - 1)*log(1 - rho))
	  
	  # create z.2 and z.cross terms for the exponential term
	  # z.2 = sum_{i=1}^{n_j} z_i^2
	  z.2 <- sum(zj^2)
	  # z.cross = sum_{i=1}^{n_j} sum_{k \neq i} z_i z_k
	  z.cross <- 0
	  for(i in 1:nj){
	    for(k in 1:nj){
	      if(i != k){ z.cross <- z.cross + (zj[i] * zj[k]) }
	    }
	  }
	  
	  # exponent part
	  exp.Aj <- -(1/2) * (a.1*z.2 + a.2*z.cross)
	  
	  return(det.Aj + exp.Aj)
	}
	
	# initialize acceptance rate counter variables
	rho.n.iter <- 0
	rho.accept <- 0

	############ DAPHNE: load parameter traces from first.stage.directory ############
	# in second stage, beta and rho_bc are sampled conditionally on the parameters of bayesTFR from the first stage of estimation
	# load in the parameter traces from first.stage.directory, which contains a converged run of default (unconditional) bayesTFR
	mean_eps_tau.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("mean_eps_tau"),
	                                                burnin = mcenv$meta$first.stage.burnin)
	sd_eps_tau.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("sd_eps_tau"),
	                                              burnin = mcenv$meta$first.stage.burnin)
	chi.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("chi"), burnin = mcenv$meta$first.stage.burnin)
	psi.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("psi"), burnin = mcenv$meta$first.stage.burnin)
	Triangle4.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("Triangle4"),
	                                             burnin = mcenv$meta$first.stage.burnin)
	delta4.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("delta4"),
	                                          burnin = mcenv$meta$first.stage.burnin)
	alpha.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("alpha"),
	                                         burnin = mcenv$meta$first.stage.burnin)
	delta.traces <- get.tfr.parameter.traces(m.default$mcmc.list, c("delta"),
	                                         burnin = mcenv$meta$first.stage.burnin)
	
	# if using Phase II AR(1) 
	if(!is.null(mcenv$meta$ar.phase2) && mcenv$meta$ar.phase2){
	  rho_phase2.traces <- get.tfr.parameter.traces(m.default$mcmc.list, "rho_phase2", 
	                                                burnin = mcenv$meta$first.stage.burnin)
	}
	
	# if m.default used uncertainty = TRUE, also get the Phase III traces
	if(mcenv$second.stage.uncertainty){
	  mu.traces <- get.tfr3.parameter.traces(m.default.3$mcmc.list, c("mu"), 
	                                         burnin = mcenv$meta$first.stage.burnin)
	  rho.traces <- get.tfr3.parameter.traces(m.default.3$mcmc.list, c("rho"), 
	                                          burnin = mcenv$meta$first.stage.burnin)
	  sigma.mu.traces <- get.tfr3.parameter.traces(m.default.3$mcmc.list, c("sigma.mu"), 
	                                               burnin = mcenv$meta$first.stage.burnin)
	  sigma.rho.traces <- get.tfr3.parameter.traces(m.default.3$mcmc.list, c("sigma.rho"), 
	                                                burnin = mcenv$meta$first.stage.burnin)
	  sigma.eps.traces <- get.tfr3.parameter.traces(m.default.3$mcmc.list, c("sigma.eps"), 
	                                                burnin = mcenv$meta$first.stage.burnin)
	}

	# sample once from parameter traces for whole iter
	iter.sample <- sample(1:nrow(delta.traces), size = nr_simu-start.iter+1, replace = TRUE)

	# country-specific parameters for Phase II
	# including past TFR if second.stage.uncertainty = TRUE
	# store samples as list, each country is separate element
	d.samples <- list()
	gamma.samples <- list()
	Trianglec4.samples <- list()
	U.samples <- list()
	if(mcenv$second.stage.uncertainty){ tfr.samples <- list() }

	for (country in id_DL_all){
	  cd <- get.country.object(country, meta = mcenv$meta, index = TRUE)$code
	  country.i.obj <- get.country.object(cd, meta = m.default$meta, index = FALSE)
	  country.traces <- get.tfr.parameter.traces.cs(m.default$mcmc.list, country.obj = country.i.obj, par.names = c("d", "gamma", "Triangle_c4"), burnin = mcenv$meta$first.stage.burnin)
	  d.samples[[country]] <- country.traces[iter.sample, 1]
	  gamma.samples[[country]] <- country.traces[iter.sample, c(2,3,4)]
	  Trianglec4.samples[[country]] <- country.traces[iter.sample, 5]
	  U.samples[[country]] <- get.tfr.parameter.traces.cs(m.default$mcmc.list, country.obj = country.i.obj, par.names = "U", burnin = mcenv$meta$first.stage.burnin)[iter.sample]
	  
	  if(mcenv$second.stage.uncertainty){
	    tfr.year <- rownames(m.default$meta$tfr_matrix)
	    tfr.samples[[country]] <- get.tfr.parameter.traces.cs(m.default$mcmc.list, country.obj = country.i.obj, par.names = c("tfr"), burnin = mcenv$meta$first.stage.burnin)[iter.sample, which(tfr.year %in% rownames(mcenv$meta$tfr_matrix))]
	  }
	}

	# country-specific parameters for Phase III 
	# store samples as list, each country is separate element
	if(mcenv$second.stage.uncertainty){
	  mu.c.samples <- list()
	  rho.c.samples <- list()
	  for(country in 1:mcenv$meta$nr.countries){
	    country.i.obj <- get.country.object(m.default.3$meta$id_phase3[country], meta = m.default.3$meta, index = TRUE)
	    country.traces <- get.tfr3.parameter.traces.cs(m.default.3$mcmc.list, country.obj = country.i.obj, par.names = c("mu.c", "rho.c"), burnin = mcenv$meta$first.stage.burnin)
	    mu.c.samples[[country]] <- country.traces[iter.sample, 1]
	    rho.c.samples[[country]] <- country.traces[iter.sample, 2]
	  }
	}

	# MCMC sampling steps for default bayesTFR parameters have all been modified to be conditional sampling based on parameter traces from m.default
  ################################################################### 
  # Start MCMC
	############
	  for (simu in start.iter:nr_simu) {
	    # DAPHNE: index for iter.sample 
	    iter.idx <- simu - start.iter + 1

	    if(verbose.iter > 0 && (simu %% verbose.iter == 0))
        	cat('\nIteration:', simu, '--', date())
        unblock.gtk('bDem.TFRmcmc')
        #################################################################
        ## a_sd, b_sd, f_sd and sigma0
        #################################################################
        # updates sd_Tc
        # start with this to get the right sd_Tc in the next steps!!

        mcmc.update.abSsigma0const(mcenv, idx.tau_c.id.notearly, trace.sample = iter.sample[iter.idx], m.default = m.default, first.stage.burnin = mcenv$meta$first.stage.burnin)

       #################################################################### 
        # 2. mean_eps_tau sd_eps_tau: gibbs step
        ##################################################################
        
        mcenv$mean_eps_tau <- mean_eps_tau.traces[iter.sample[iter.idx]]
        mcenv$sd_eps_tau <- sd_eps_tau.traces[iter.sample[iter.idx]]

 		    #update all not-early countries
		    mcenv$sd_Tc[idx.tau_c.id.notearly.all] <- mcenv$sd_eps_tau
		    mcenv$mean_eps_Tc[idx.tau_c.id.notearly.all] <- mcenv$mean_eps_tau
		    
        #################################################################### 
        # 2. chi's and psi's: gibbs step
        ##################################################################
		    mcenv$chi <- chi.traces[iter.sample[iter.idx]]
		    mcenv$psi <- psi.traces[iter.sample[iter.idx]]

        #################################################################### 
        # 2. Triangle4 and delta4_star: gibbs step
        ##################################################################
        mcenv$Triangle4 <- Triangle4.traces[iter.sample[iter.idx]]
		    mcenv$delta4 <- delta4.traces[iter.sample[iter.idx]]

        if (!is.null(mcenv$meta$ar.phase2) && mcenv$meta$ar.phase2)
        {
          var_prop <- rho_phase2.traces[iter.sample[iter.idx]]
          # update eps_T based on second.stage.uncertainty
          if(mcenv$second.stage.uncertainty){
            suppl.T.tmp <- if(!is.null(mcenv$meta$suppl.data$regions)) mcenv$meta$suppl.data$T_end else 0
            eps_prop <- matrix(NA, mcenv$meta$T_end-1 + suppl.T.tmp, mcenv$meta$nr_countries)
            beta_prop <- c(mcenv$beta_e, mcenv$beta_fp, mcenv$beta_g, 
                           mcenv$beta_e_SSA, mcenv$beta_fp_SSA, mcenv$beta_g_SSA)
            
            for (country in mcenv$meta$id_DL){
              theta_prop <- c((mcenv$U_c[country]-mcenv$Triangle_c4[country])*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), mcenv$Triangle_c4[country], mcenv$d_c[country])
              
              start_idx <- max(mcenv$meta$start_c[country], mcenv$meta$start_cov_data_c[country])
              idx <- start_idx:(mcenv$meta$lambda_c[country]-1)
              raw.outliers <- mcenv$meta$indices.outliers[[as.character(country)]]
              if (!is.null(mcenv$meta$ar.phase2) && mcenv$meta$ar.phase2) 
                raw.outliers <- sort(unique(c(raw.outliers, raw.outliers+1)))
              idx <- setdiff(idx, raw.outliers)

              eps_prop[idx, country] <- get.eps.T(theta_prop, beta_prop, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=var_prop)
            }
          } else{
            eps_prop <- get_eps_T_all(mcenv, matrix.name=matrix.name, rho.phase2=var_prop)
          }
          
          mcenv$rho.phase2 <- var_prop
          mcenv$eps_Tc <- eps_prop
        }
		    
        #################################################################### 
        # country-specific parameters: d_c, gamma's, U_c and Triangle_c4
        ##################################################################
		    beta_prop <- c(mcenv$beta_e, mcenv$beta_fp, mcenv$beta_g, 
		                   mcenv$beta_e_SSA, mcenv$beta_fp_SSA, mcenv$beta_g_SSA)
        for (country in id_DL_all){
          start_idx <- max(mcenv$meta$start_c[country], mcenv$meta$start_cov_data_c[country])
          idx <- start_idx:(mcenv$meta$lambda_c[country]-1)
          
          ##### D
          theta_prop <-  c((mcenv$U_c[country]-mcenv$Triangle_c4[country])*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), mcenv$Triangle_c4[country], mcenv$d_c[country])
          d_prop <- d.samples[[country]][iter.idx]
          if(mcenv$second.stage.uncertainty){
            eps_T_prop <- get.eps.T(c(theta_prop[-5], d_prop), beta_prop, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          } else{
            eps_T_prop <- get.eps.T(c(theta_prop[-5], d_prop), beta_prop, tfr_trace = NULL, country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          }
          
          mcenv$eps_Tc[idx, country] <- eps_T_prop
          mcenv$d_c[country] <- d_prop

          ##### GAMMA
          gamma_prop <- gamma.samples[[country]][iter.idx, ]
          pci_prob <- exp(gamma_prop)/sum(exp(gamma_prop))
          theta_prop <- c(pci_prob*(mcenv$U_c[country] - mcenv$Triangle_c4[country]), 
                          mcenv$Triangle_c4[country], mcenv$d_c[country]) 
          if(mcenv$second.stage.uncertainty){
            eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          } else{
            eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = NULL, country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          }
          
          mcenv$gamma_ci[country,] <- gamma_prop
          mcenv$eps_Tc[idx, country] <- eps_T_prop
          
          ##### TRIANGLEC4
          Triangle_c4_prop <- Trianglec4.samples[[country]][iter.idx]
          theta_prop <-  c((mcenv$U_c[country]-Triangle_c4_prop)*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), Triangle_c4_prop, mcenv$d_c[country])
          if(mcenv$second.stage.uncertainty){
            eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          } else{
            eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = NULL, country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
          }

          mcenv$eps_Tc[idx, country] <- eps_T_prop
          mcenv$Triangle_c4[country] <- Triangle_c4_prop
        
		      ##### U
          theta_prop <- c((mcenv$U_c[country]-mcenv$Triangle_c4[country])*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), mcenv$Triangle_c4[country], mcenv$d_c[country])
          U_prop <- U.samples[[country]][iter.idx]
          theta_prop[1:3] <- theta_prop[1:3]/(mcenv$U_c[country] - mcenv$Triangle_c4[country])*(U_prop - mcenv$Triangle_c4[country])
          if(mcenv$second.stage.uncertainty){
            eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
            } else{
              eps_T_prop <- get.eps.T(theta_prop, beta_prop, tfr_trace = NULL, country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
              }
           
          mcenv$eps_Tc[start_idx:(mcenv$meta$lambda_c[country]-1), country] <- eps_T_prop
          mcenv$U_c[country] <- U_prop
        } 

         ##################################################################
         #mcenv# alpha_i's and delta_i's, with gibbs step
         ##################################################################
         mcenv$alpha <- alpha.traces[iter.sample[iter.idx], ]
		     mcenv$delta <- delta.traces[iter.sample[iter.idx], ]

         ##################################################################
         # DAPHNE: bc_rho, beta
         ##################################################################
         #### compute outcome Y = TFR Decr + DL ####
         # using current TFR to compute DL
         Y.matrix <- matrix(nrow = nrow(mcenv$meta$tfr_matrix), ncol = ncol(mcenv$meta$tfr_matrix))
         rownames(Y.matrix) <- rownames(mcenv$meta$tfr_matrix)
         colnames(Y.matrix) <- colnames(mcenv$meta$tfr_matrix)

         for(country in id_DL_all){
           # vector of observed non-NA TFR values
           # if have NAs, then start of Phase II will be earlier than start of covariate data
           start_idx <- max(mcenv$meta$start_c[country], mcenv$meta$start_cov_data_c[country])
           if(mcenv$second.stage.uncertainty){
             tfr.vctr <- tfr.samples[[country]][iter.idx, start_idx:mcenv$meta$lambda_c[country]]
             names(tfr.vctr) <- names(get.observed.tfr(country, mcenv$meta)[start_idx:mcenv$meta$lambda_c[country]])
           } else{
             tfr.vctr <- get.observed.tfr(country, mcenv$meta)[start_idx:mcenv$meta$lambda_c[country]] # matrix = matrix.name?
           }
           
           ldl <- length(tfr.vctr)-1
           
           # expected decrements for the observed TFR
           theta.pars <- c((mcenv$U_c[country]-mcenv$Triangle_c4[country])*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), mcenv$Triangle_c4[country], mcenv$d_c[country])
           DL.obs <- DLcurve(theta.pars, tfr.vctr[1:ldl], mcenv$meta$dl.p1, mcenv$meta$dl.p2, annual = mcenv$meta$annual.simulation)
           names(DL.obs) <- names(tfr.vctr)[2:(ldl+1)]
           
           # vector of TFR Decr
           tfr.decr.vctr <- tfr.vctr[2:(ldl+1)] - tfr.vctr[1:ldl]
           
           tfr.decr.vctr <- tfr.decr.vctr[!is.na(tfr.decr.vctr)]
           tfr.decr.obs.years <- names(tfr.decr.vctr)
           obs.years <- intersect(names(DL.obs), tfr.decr.obs.years)
           
           # compute outcome Y = TFR Decr + Expected Decr
           Y.matrix[obs.years, country] <- tfr.decr.vctr + DL.obs
         }
         
         #### set up vectors ####
         Y.v <- as.vector(t(Y.matrix))[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
         X.mv <- cbind(X.educ, X.fp, X.g, SSA.ind*X.educ, SSA.ind*X.fp, SSA.ind*X.g)
         # for indices to match up for SD vector, add in the first time period row to mcenv$sd_Tc as all NAs
         sd_Tc_temp <- rbind(rep(NA, ncol(mcenv$sd_Tc)), mcenv$sd_Tc)
         sd.v <- as.vector(t(sd_Tc_temp))[idxX.educ & idxX.fp & idxtfr.decr & idxX.g]
         beta.v <- c(mcenv$beta_e, mcenv$beta_fp, mcenv$beta_g,
                     mcenv$beta_e_SSA, mcenv$beta_fp_SSA, mcenv$beta_g_SSA)
         
         #### bc_rho: Metropolis-Hastings step ####
         # proposal function: truncated normal from truncnorm
         rho.old <- mcenv$bc_rho
         rho.new <- rtruncnorm(1, a=0, b=1, mean=rho.old, sd=0.04)
         
         # acceptance ratio = exp(loglik_new + prior_new - loglik_old - prior_old)
         loglik.new <- rho.post.ll(rho.new, X.mv, Y.v, beta.v, sd.v, new.clust)
         loglik.old <- rho.post.ll(rho.old, X.mv, Y.v, beta.v, sd.v, new.clust)
         prior.new <- log(dtruncnorm(1, a=0, b=1, mean=rho.new, sd=0.04))
         prior.old <- log(dtruncnorm(1, a=0, b=1, mean=rho.old, sd=0.04))
         rho.AR <- exp(loglik.new + prior.new - loglik.old - prior.old)
         
         # acceptance rate = fraction of proposed samples that is accepted
         rho.n.iter <- rho.n.iter + 1
         if(runif(1, 0, 1) <= rho.AR){
           mcenv$bc_rho <- rho.new
           rho.accept <- rho.accept + 1
         } # else, reject (mcenv$bc_rho does not chage)
         
         # check acceptance rate at end of each chain
         # if(rho.n.iter == (nr_simu-1)){
         #   print(paste("acceptance rate =", rho.accept/rho.n.iter))
         # }
         
         #### beta: Gibbs step ####
         # multivariate posterior of beta, using correlation matrix R 
         
         # get matrix R^{-1} from rho
         R.inv <- get.R.from.rho(mcenv$bc_rho, length(sd.v), new.clust)
         # convert to sparse matrix using Matrix library
         R.inv <- as(R.inv, "sparseMatrix")
         diag.sd <- as(diag(1/sd.v), "sparseMatrix")
         
         # create matrix \Sigma^{-1}
         Sigma.inv <- diag.sd %*% R.inv %*% diag.sd
         # convert back to normal matrix
         Sigma.inv <- as.matrix(Sigma.inv)
         
         # posterior mean and SD of beta
         # var = (X^T \Sigma^{-1} X + \Omega^{-1})^{-1}
         # mean = (var) X^T \Sigma^{-1} Y 
         X.Sigma.inv <- crossprod(X.mv, Sigma.inv) # = t(X.mv) %*% Sigma.inv
         sXY.mv <- X.Sigma.inv %*% Y.v
         XX.mv <- X.Sigma.inv %*% X.mv
         denom.mv <- solve(XX.mv + diag(c(sdbeta.inv.e, sdbeta.inv.fp, sdbeta.inv.g, 
                                          sdbeta.inv.e.SSA, 
                                          sdbeta.inv.fp.SSA, sdbeta.inv.g.SSA)))

         mv.sample <- rmvnorm(1, mean=crossprod(sXY.mv, denom.mv), sigma = denom.mv)
         mcenv$beta_e <- mv.sample[1]
         mcenv$beta_fp <- mv.sample[2]
         mcenv$beta_g <- mv.sample[3]
         mcenv$beta_e_SSA <- mv.sample[4]
         mcenv$beta_fp_SSA <- mv.sample[5]
         mcenv$beta_g_SSA <- mv.sample[6]
         
         # save sampled iterations of first stage's TFR as part of mcenv
         mcenv$sampled_iter <- ifelse(mcenv$second.stage.uncertainty, 
                                      iter.sample[iter.idx], 0)
         
         # update eps_T after all parameters are estimated for this iter
         beta_for_eps_T <- c(mcenv$beta_e, mcenv$beta_fp, mcenv$beta_g, 
                             mcenv$beta_e_SSA, mcenv$beta_fp_SSA, mcenv$beta_g_SSA)
         for(country in id_DL_all){
           theta_for_eps_T <-  c((mcenv$U_c[country]-mcenv$Triangle_c4[country])*exp(mcenv$gamma_ci[country,])/sum(exp(mcenv$gamma_ci[country,])), mcenv$Triangle_c4[country], mcenv$d_c[country])
           start_idx <- max(mcenv$meta$start_c[country], mcenv$meta$start_cov_data_c[country])
           if(mcenv$second.stage.uncertainty){
             temp_eps <- get.eps.T(theta_for_eps_T, beta_for_eps_T, tfr_trace = tfr.samples[[country]][iter.idx, ], country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
           } else{
             temp_eps <- get.eps.T(theta_for_eps_T, beta_for_eps_T, tfr_trace = NULL, country, mcenv$meta, rho.phase2=mcenv$rho.phase2)
           }
           mcenv$eps_Tc[start_idx:(mcenv$meta$lambda_c[country]-1), country] <- temp_eps
         }
         

         ##################################################################
         # Phase III
         ##################################################################
         
         if (uncertainty)
         {
           one.step.mcmc3.sampling(mcenv)
           
           mcmc.update.tfr.year(mcenv)
         }

         # Daphne: if using second.stage.uncertainty, need to update Phase III here
         if (mcenv$second.stage.uncertainty){
           mcenv$mu <- mu.traces[iter.sample[iter.idx]]
           mcenv$rho <- rho.traces[iter.sample[iter.idx]]
           mcenv$sigma.mu <- sigma.mu.traces[iter.sample[iter.idx]]
           mcenv$sigma.rho <- sigma.rho.traces[iter.sample[iter.idx]]
           mcenv$sigma.eps <- sigma.eps.traces[iter.sample[iter.idx]]
           
           for(country in 1:mcenv$meta$nr.countries){
             mcenv$mu.c[country] <- mu.c.samples[[country]][iter.idx]
             mcenv$rho.c[country] <- rho.c.samples[[country]][iter.idx]
           }
           
           # mcmc.update.tfr.year(mcenv)
         }
         ################################################################### 
         # write samples simu/thin to disk
         ##################################################################
         mcenv$finished.iter <- simu
         mcenv$rng.state <- .Random.seed
         
         if (simu %% thin == 0){
           mcenv$length <- mcenv$length + 1
           flush.buffer <- FALSE
           if (simu + thin > nr_simu) flush.buffer <- TRUE
           store.mcmc(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
           # Daphne: edited to be uncertainty OR second.stage.uncertainty
           if (uncertainty | mcenv$second.stage.uncertainty) store.mcmc3(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
         }
	}       # end simu loop MCMC
	.cleanup.mcmc(mcenv)
	resmc <- as.list(mcenv)
	store.bayesTFR.meta.object(mcenv$meta, mcenv$meta$output.dir)
	
	class(resmc) <- class(mcmc)
    return(resmc)
}

one.step.mcmc3.sampling <- function(mcmc)
{
  meta <- mcmc$meta
  niter <- mcmc$iter
  nr.countries <- meta$nr.countries
  countries.index <- meta$id_phase3
  ardata <- mcmc$observations
  Ts <- mcmc$length_obs
  gamma.mu.low <- 1/(meta$sigma.mu.prior.range[2])^2
  gamma.mu.up <- if (meta$sigma.mu.prior.range[1] == 0) NA else 1/(meta$sigma.mu.prior.range[1])^2
  gamma.rho.low <- 1/(meta$sigma.rho.prior.range[2])^2
  gamma.rho.up <- if (meta$sigma.rho.prior.range[1] == 0) NA else 1/(meta$sigma.rho.prior.range[1])^2
  gamma.eps.low <- 1/(meta$sigma.eps.prior.range[2])^2
  gamma.eps.up <- if (meta$sigma.eps.prior.range[1] == 0) NA else 1/(meta$sigma.eps.prior.range[1])^2
  
  #### Start MCMC
  
  mu.integral.to.mC <- (mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0))^(-nr.countries)
  prop.mu <- proposal.mu.rho(mcmc[['mu.c']], mcmc[['sigma.mu']], nr.countries, 
                             meta[['mu.prior.range']][1], meta[['mu.prior.range']][2])
  accept.prob <- min(((mu.rho.integral(prop.mu, mcmc[['sigma.mu']], low=0))^(-nr.countries))/mu.integral.to.mC, 1)
  if (runif(1) < accept.prob) {
    mcmc[['mu']] <- prop.mu
    mcmc[['recompute.mu.integral']] <- TRUE
  } else mcmc[['recompute.mu.integral']] <- FALSE
  
  # Metropolis-Hastings for sigma_mu=1/sqrt(lambda_mu)
  S <- sum((mcmc[['mu.c']]-mcmc[['mu']])^2)
  if(mcmc[['recompute.mu.integral']]) mu.integral.to.mC <- mu.rho.integral(mcmc[['mu']], mcmc[['sigma.mu']], low=0)^(-nr.countries)
  prop.lambda.mu <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.mu.low, high=gamma.mu.up)
  accept.prob <- min(((mu.rho.integral(mcmc[['mu']], 1/prop.lambda.mu, low=0))^(-nr.countries))/mu.integral.to.mC, 1)
  if (runif(1) < accept.prob) {
    mcmc[['sigma.mu']] <- 1/sqrt(prop.lambda.mu)
    recompute.mu.integral <- TRUE
  } else recompute.mu.integral <- FALSE	
  
  # Slice sampling for rho
  mcmc[['rho']] <- slice.sampling(mcmc[['rho']], logdensity.mu.rho, 1, 
                                  low=meta[['rho.prior.range']][1], up=meta[['rho.prior.range']][2], 
                                  par.c=mcmc[['rho.c']], sd=mcmc[['sigma.rho']], 
                                  c.low=0, c.up=1)
  # Metropolis-Hastings for rho
  # rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
  # prop.rho <- proposal.mu.rho(mcmc[['rho.c']], mcmc[['sigma.rho']], nr.countries, 
  # meta[['rho.prior.range']][1], meta[['rho.prior.range']][2])
  # accept.prob <- min(((mu.rho.integral(prop.rho, mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
  # if (runif(1) < accept.prob) {
  # mcmc[['rho']] <- prop.rho
  # recompute.rho.integral <- TRUE
  # } else recompute.rho.integral <- FALSE
  
  mcmc[['sigma.rho']] <- slice.sampling(mcmc[['sigma.rho']], logdensity.sigma.mu.rho, 1, 
                                        low=meta$sigma.rho.prior.range[1], up=meta$sigma.rho.prior.range[2], 
                                        par.c=mcmc[['rho.c']], mean=mcmc[['rho']],
                                        c.low=0, c.up=1)
  
  # Metropolis-Hastings for sigma_rho=1/sqrt(lambda_rho)
  # S <- sum((mcmc[['rho.c']]-mcmc[['rho']])^2)
  # if(recompute.rho.integral) 
  # rho.integral.to.mC <- (mu.rho.integral(mcmc[['rho']], mcmc[['sigma.rho']], low=0, up=1))^(-nr.countries)
  # prop.lambda.rho <- rgamma.trunc((nr.countries-1)/2, S/2, low=gamma.rho.low, high=gamma.rho.up)
  # accept.prob <- min(((mu.rho.integral(mcmc[['rho']], 1/prop.lambda.rho, low=0, up=1))^(-nr.countries))/rho.integral.to.mC, 1)
  # if (runif(1) < accept.prob) {
  # mcmc[['sigma.rho']] <- 1/sqrt(prop.lambda.rho)
  # recompute.rho.integral <- TRUE
  # } else recompute.rho.integral <- FALSE
  
  sigma.eps.sq <- mcmc$sigma.eps^2
  sigma.mu.sq <- mcmc$sigma.mu^2
  sigma.mu.sq.inv <- 1/sigma.mu.sq
  mu.over.sigma.sq <- mcmc$mu/sigma.mu.sq
  sigma.rho.sq <- mcmc$sigma.rho^2
  sigma.rho.sq.inv <- 1/sigma.rho.sq
  rho.over.sigma.sq <- mcmc$rho/sigma.rho.sq
  S.eps <- STn <- 0
  one.minus.rho <- 1-mcmc$rho.c
  one.minus.rho.sq <- one.minus.rho^2
  W <- one.minus.rho.sq/sigma.eps.sq
  
  # country-specific parameters - Gibbs sampler
  for(country in 1:nr.countries) {
    f.ct <- ardata[[country]][2:Ts[country]]
    f.ctm1 <- ardata[[country]][1:(Ts[country]-1)]
    # mu.c
    s <- sum((f.ct - mcmc$rho.c[country]*f.ctm1)/one.minus.rho[country])
    nomin <- W[country] * s + mu.over.sigma.sq
    denom <- (Ts[country]-1) * W[country] + sigma.mu.sq.inv
    mcmc$mu.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0)	
    
    # rho.c
    d1 <- f.ctm1 - mcmc$mu.c[country]
    a <- sum(d1^2)/sigma.eps.sq
    b <- sum(d1*(f.ct - mcmc$mu.c[country]))/sigma.eps.sq
    nomin <- b + rho.over.sigma.sq
    denom <- a + sigma.rho.sq.inv
    mcmc$rho.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0, #high=1-10*.Machine$double.xmin
                                       high=0.999999)
    S.eps <- S.eps + sum((f.ct - mcmc$mu.c[country] - mcmc$rho.c[country]*d1)^2)
    STn <- STn + Ts[country]-1
  }
  # Gibbs for sigma.eps
  mcmc$sigma.eps <- 1/sqrt(rgamma.trunc((STn-1)/2, S.eps/2, low=gamma.eps.low, high=gamma.eps.up))
  
}

.cleanup.mcmc <- function(mcmc) {
	if(is.environment(mcmc)) {
		rm(list=mcmc$dontsave[mcmc$dontsave != 'meta' & mcmc$dontsave %in% ls(mcmc)], envir=mcmc)
		return(NULL)
	}
	for(rmitem in mcmc$dontsave[mcmc$dontsave != 'meta'  & mcmc$dontsave %in% names(mcmc)]) mcmc[[rmitem]] <- NULL
	return(mcmc)
}

unblock.gtk <- function(option, options.list=NULL) {
	if(!getOption(option, default=FALSE)) return()
	if(!is.null(options.list)) options(options.list)
	# This is to unblock the GUI, if the run is invoked from bayesDem
	# In such a case the gtk libraries are already loaded
	while(do.call('gtkEventsPending', list()))
		do.call('gtkMainIteration', list())

}

tfr.mcmc.sampling.extra <- function(mcmc, mcmc.list, countries, posterior.sample,
											 iter=NULL, thin = 1, burnin=2000, verbose=FALSE, verbose.iter=100, uncertainty=FALSE) {
	#run mcmc sampling for countries given by the index 'countries'
  nr_simu <- iter
	if (is.null(iter))
    	nr_simu <- mcmc$length
    nr_countries_all <- mcmc$meta$nr_countries
    nr_countries <- length(countries)
    Triangle_c4.low <- mcmc$meta$Triangle_c4.low
    Triangle_c4.up <- mcmc$meta$Triangle_c4.up
    lower_d <- mcmc$meta$d.low
    upper_d <- mcmc$meta$d.up

    id_DL <- mcmc$meta$id_DL[is.element(mcmc$meta$id_DL, countries)]
    id_early <- mcmc$meta$id_early[is.element(mcmc$meta$id_early, countries)]
    id_notearly <- mcmc$meta$id_notearly[is.element(mcmc$meta$id_notearly, countries)]
	    
	suppl.T <- if(!is.null(mcmc$meta$suppl.data$regions)) mcmc$meta$suppl.data$T_end else 0
	const_sd_dummie_Tc_extra <- matrix(0, mcmc$meta$T_end-1+suppl.T, nr_countries)
	mid.years <- as.integer(c(if(suppl.T > 0) rownames(mcmc$meta$suppl.data$tfr_matrix) else c(), rownames(mcmc$meta$tfr_matrix)))
	const_sd_dummie_Tc_extra[mid.years[-length(mid.years)] < 1975,] <- 1
	
	
	# get values of the hyperparameters (sample from the posterior)
  hyperparameter.names <- tfr.parameter.names(trans=FALSE)
  hyperparameters <- list()
  sampled.index <- sample(posterior.sample, nr_simu, replace=TRUE)
  th.burnin <- get.thinned.burnin(mcmc, burnin)
  if (uncertainty)
  {
    mcmc3.format <- list()
    for (i in 1:length(mcmc.list))
    {
      mcmc3.format[[i]] <- list()
      mcmc3.format[[i]]$meta$phase <- 3
      mcmc3.format[[i]]$compression.type <- mcmc.list[[i]]$compression.type
      mcmc3.format[[i]]$meta$compression.type <- mcmc.list[[i]]$meta$compression.type
      mcmc3.format[[i]]$meta$output.dir <- file.path(mcmc.list[[i]]$meta$output.dir, 'phaseIII')
      mcmc3.format[[i]]$output.dir <- mcmc.list[[i]]$output.dir
      mcmc3.format[[i]]$traces <- mcmc.list[[i]]$traces
      mcmc3.format[[i]] <- structure(mcmc3.format[[i]], class='bayesTFR.mcmc')
    }
    hyperparameter3.names <- tfr3.parameter.names()
    for (par in hyperparameter3.names) {
      hyperparameters[[par]] <- c()
      for(mc in mcmc3.format) {
        if (no.traces.loaded(mc)  || th.burnin < mc$traces.burnin) {
          traces <- bdem.parameter.traces(mc, par, burnin=th.burnin)
        } else {
          traces <- get.burned.tfr.traces(mc, par, th.burnin)
        }
        hyperparameters[[par]] <- rbind(hyperparameters[[par]], traces)
      }
      hyperparameters[[par]] <- hyperparameters[[par]][sampled.index,]
    }
  }
  if (!is.null(mcmc$meta$ar.phase2) && mcmc$meta$ar.phase2) hyperparameter.names <- c(hyperparameter.names, 'rho_phase2')
  else mcmc$rho.phase2 <- NULL
  for (par in hyperparameter.names) {
    hyperparameters[[par]] <- c()
    for(mc in mcmc.list) {
    	if (no.traces.loaded(mc)  || th.burnin < mc$traces.burnin) {
    		traces <- bdem.parameter.traces(mc, par, burnin=th.burnin)
        } else {
        		traces <- get.burned.tfr.traces(mc, par, th.burnin)
       	}
       	hyperparameters[[par]] <- rbind(hyperparameters[[par]], traces)
      }
    hyperparameters[[par]] <- hyperparameters[[par]][sampled.index,]
  }
  
	mcmc.orig <- mcmc
	if (uncertainty)
	{
	  meta <- mcmc$meta
	  nr.countries <- meta$nr.countries
	  if (nr.countries > 0)
	  {
	    countries.index <- meta$id_phase3
	    ardata <- list()
	    Ts <- rep(0, nr.countries)
	    for(country in 1:nr.countries) {
	      data <- get.observed.tfr(countries.index[country], meta, 'tfr_matrix_all')
	      ardata[[country]] <- data[meta$lambda_c[countries.index[country]]:meta$T_end_c[countries.index[country]]]
	      Ts[country] <- length(ardata[[country]])
	    }
	    mcmc$observations <- ardata
	    mcmc$length_obs <- Ts
	  }
	}
	matrix.name <- ifelse(uncertainty, 'tfr_all', 'tfr_matrix')
	
	mcenv <- as.environment(mcmc) # Create an environment for the mcmc stuff in order to avoid 
					              # copying of the mcmc list 
  
	mcenv$const_sd_dummie_Tc <- matrix(0, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
	mcenv$const_sd_dummie_Tc[mid.years[-length(mid.years)] < 1975,] <- 1
	
	updated.var.names <- c('gamma_ci', 'd_c', 'Triangle_c4', 'U_c')
	if (uncertainty) updated.var.names <- c(updated.var.names, 'rho.c', 'mu.c')
	idx.tau_c.id.notearly <- matrix(c(mcmc$meta$tau_c[id_notearly], id_notearly), ncol=2)

    ################################################################### 
    # Start MCMC
	############
    for (simu in 1:nr_simu) {
        if(verbose.iter > 0 && (simu %% verbose.iter == 0))
			cat('\nIteration:', simu, '--', date())
			unblock.gtk('bDem.TFRmcmcExtra')
        # set hyperparameters for this iteration
        for (par in hyperparameter.names) {
        	if(is.null(dim(hyperparameters[[par]]))) {
        		mcenv[[par]] <- hyperparameters[[par]][simu]
        	} else {
        		mcenv[[par]] <- hyperparameters[[par]][simu,]
        	}
        }
			  
        # compute eps_T, mean_eps_t and sd_Tc
        if(is.null(mcenv$eps_Tc)) mcenv$eps_Tc <- get_eps_T_all(mcenv)
         	
        add_to_sd_Tc_extra <- matrix(NA, mcenv$meta$T_end-1 + suppl.T, nr_countries)
        mcenv$data.list <- list()
    	for (icountry in 1:length(countries)){
    		country <- countries[icountry]
    		this.data <- get.observed.tfr(country, mcenv$meta, matrix.name = matrix.name)
    		this.data <- this.data[1:(mcenv$meta$T_end_c[country]-1)]
    		mcenv$data.list[[country]] <- this.data
			add_to_sd_Tc_extra[1:(mcenv$meta$T_end_c[country]-1),icountry] <- (
						this.data - mcenv$S_sd)*ifelse(this.data > mcenv$S_sd, -mcenv$a_sd, mcenv$b_sd)
    	}
      mcenv$add_to_sd_Tc <- matrix(0, mcenv$meta$T_end-1+suppl.T, nr_countries_all)
        
      mcenv$add_to_sd_Tc[,countries] <- add_to_sd_Tc_extra
		  mcenv$mean_eps_Tc <- matrix(0, mcenv$meta$T_end -1 + suppl.T, nr_countries_all)
        mcenv$sd_Tc <- matrix(NA, mcenv$meta$T_end -1 + suppl.T, nr_countries_all)
       	mcenv$sd_Tc[,countries] <- ifelse(const_sd_dummie_Tc_extra==1, 
         						mcenv$const_sd, 1)*
            			   ifelse((mcenv$sigma0 + add_to_sd_Tc_extra)>0, mcenv$sigma0 + add_to_sd_Tc_extra, 
            						mcenv$meta$sigma0.min)
         	
 		mcenv$sd_Tc[idx.tau_c.id.notearly] <- mcenv$sd_eps_tau
 		mcenv$mean_eps_Tc[idx.tau_c.id.notearly] <- mcenv$mean_eps_tau
        #################################################################### 
        # country-specific parameters: d_c, gamma's, U_c and Triangle_c4
        ##################################################################
        for (country in id_DL) {
        	mcmc.update.d(country, mcenv, matrix.name=matrix.name, rho.phase2=mcenv$rho.phase2)
          mcmc.update.gamma(country, mcenv, matrix.name=matrix.name, rho.phase2=mcenv$rho.phase2)
          mcmc.update.Triangle_c4(country, mcenv, matrix.name=matrix.name, rho.phase2=mcenv$rho.phase2)
        } 
 
         # U_c updated only for countries with early decline
         for (country in id_early){
                mcmc.update.U(country, mcenv, matrix.name=matrix.name, rho.phase2=mcenv$rho.phase2)
         } 
 		
 		    if (uncertainty && (nr.countries > 0))
 		    {
 		      for (par in hyperparameter3.names) {
 		        if(is.null(dim(hyperparameters[[par]]))) {
 		          mcenv[[par]] <- hyperparameters[[par]][simu]
 		        } else {
 		          mcenv[[par]] <- hyperparameters[[par]][simu,]
 		        }
 		      }
 		      sigma.eps.sq <- mcenv$sigma.eps^2
 		      sigma.mu.sq <- mcenv$sigma.mu^2
 		      sigma.mu.sq.inv <- 1/sigma.mu.sq
 		      mu.over.sigma.sq <- mcenv$mu/sigma.mu.sq
 		      sigma.rho.sq <- mcenv$sigma.rho^2
 		      sigma.rho.sq.inv <- 1/sigma.rho.sq
 		      rho.over.sigma.sq <- mcenv$rho/sigma.rho.sq
 		      S.eps <- STn <- 0
 		      one.minus.rho <- 1-mcenv$rho.c
 		      one.minus.rho.sq <- one.minus.rho^2
 		      W <- one.minus.rho.sq/sigma.eps.sq
 		      for(country in 1:nr.countries) {
 		        f.ct <- ardata[[country]][2:Ts[country]]
 		        f.ctm1 <- ardata[[country]][1:(Ts[country]-1)]
 		        # mu.c
 		        s <- sum((f.ct - mcenv$rho.c[country]*f.ctm1)/one.minus.rho[country])
 		        nomin <- W[country] * s + mu.over.sigma.sq
 		        denom <- (Ts[country]-1) * W[country] + sigma.mu.sq.inv
 		        mcenv$mu.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0)	
 		        
 		        # rho.c
 		        d1 <- f.ctm1 - mcenv$mu.c[country]
 		        a <- sum(d1^2)/sigma.eps.sq
 		        b <- sum(d1*(f.ct - mcenv$mu.c[country]))/sigma.eps.sq
 		        nomin <- b + rho.over.sigma.sq
 		        denom <- a + sigma.rho.sq.inv
 		        mcenv$rho.c[country] <- rnorm.trunc(mean=nomin/denom, sd=1/sqrt(denom), low=0, #high=1-10*.Machine$double.xmin
 		                                           high=0.999999)
 		        S.eps <- S.eps + sum((f.ct - mcenv$mu.c[country] - mcenv$rho.c[country]*d1)^2)
 		        STn <- STn + Ts[country]-1
 		      }
 		    }
 		    
 		    if (uncertainty)
 		    {
 		      if (length(countries) > 3) mcmc.update.tfr.year(mcenv, countries)
 		      else
 		      {
 		        for (country in countries)
 		        {
 		          mcmc.update.tfr(country, mcenv)
 		        }
 		      }
 		    }

         ################################################################### 
         # write samples to disk
         ##################################################################
         #update the original mcmc with the new values
         for(var in updated.var.names) {
         	mcmc.orig[[var]] <- mcenv[[var]]
         }
 		    if (uncertainty)
 		    {
 		      mcmc.orig[['meta']][['tfr_all']][, countries] <- mcenv[['meta']][['tfr_all']][, countries]
 		    }
         flush.buffer <- FALSE
         append <- TRUE
		 if (simu <= thin) {
			append <- FALSE
			flush.buffer <- TRUE
		 } else {
			if (simu + thin > nr_simu) flush.buffer <- TRUE
		 }
         if (simu %% thin == 0){
         store.mcmc(mcmc.orig, append=append, flush.buffer=flush.buffer, countries=countries, 
         				verbose=verbose)
         if (uncertainty && (length(mcmc.orig$meta$id_phase3)>0)) 
           store.mcmc3(mcmc.orig, append=append, flush.buffer=flush.buffer, countries=1:length(mcmc.orig$meta$id_phase3), 
                                      verbose=verbose)}
         
	}       # end simu loop MCMC
	mcmc.orig <- .cleanup.mcmc(mcmc.orig)
    return(mcmc.orig)
}


