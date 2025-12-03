#####################################################################
####### CODE TO RECREATE FIGURE 3 FROM FALLA & CUNNIFFE PAPER #######
#####################################################################

## This script runs code to produce and save a pdf of Figure 3 to the current 
## directory

rm(list = ls())

library(deSolve)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(numDeriv)
library(cowplot)
library(latex2exp)

#### FUNCTIONS ####
transient_colonising_ode <- function(times, y, parms, death_linked_v = F) {
  
  ## model ODE function
  
  # parameters
  k <- parms[["k"]]
  t <- parms[["t"]]
  H <- parms[["H"]]
  b <- parms[["b"]]
  a <- parms[["a"]]
  gamma <- parms[["gamma"]]
  v <- parms[["v"]]
  rho <- parms[["rho"]]
  q <- parms[["q"]] # transient aphid emigration probability
  q_c <- parms[["q_c"]] # colonising aphid emigration probability
  lambda <- parms[["lambda"]] 
  psi <- parms[["psi"]] # proportion of immigrants that are transient
  pi <- parms[["pi"]]
  omega <- parms[["omega"]]
  epsilon <- parms[["epsilon"]]
  eta <- parms[["eta"]]
  sigma <- parms[["sigma"]] # colonising aphid growth rate
  theta <- parms[["theta"]] # colonising aphid death rate
  kappa <- parms[["kappa"]] # colonising carrying capacity per plant
  mu <- parms[["mu"]] # severity of plant death-linked attractiveness
  mu2 <- parms[["mu2"]] # severity of reduced inoculation-linked attractiveness
  
  # states
  I <- y[["I"]]
  X_t <- y[["X_t"]]
  Z_t <- y[["Z_t"]]
  X_c <- y[["X_c"]]
  Z_c <- y[["Z_c"]]
  
  S <- H - I
  
  # aphid flight rates (c = colonising, t = transient)
  phi_t <- 1 / (k + t)
  phi_c <- (S + v*I) / 
    (omega*eta*(S + v*epsilon*I) + (t + k)*(S + v*I))
  # phi_c <- phi_t + 
  #   (S + v*I) / 
  # (omega*eta*(S + v*epsilon*I))
  
  # aphid infectivity loss rates (c = colonising, t = transient)
  tau_t <- (1 - q)*phi_t*rho*(S + (1 - a)*v*I) / (S + v*I)
  tau_c <- (1 - q_c)*phi_c*
    (S*(rho + (1 - rho)*omega) + 
       v*I*(rho*(1 - a*(1 - epsilon*omega)) + (1 - rho)*epsilon*omega)) / (S + v*I)
  
  # PLANTS EQUATION
  # if (death_linked_v) { plant_death <- gamma*v*I }
  # else { plant_death <- gamma*I } 
  plant_death <- gamma*(1 + mu*(v-1)^2)
  inoculation <- b*(1/(1 + mu2*(v-1)^2))
  
  dI <- ((1 - q)*phi_t*Z_t + (1 - q_c)*phi_c*Z_c)*inoculation*S/(S + v*I) - 
    plant_death*I
  
  # TRANSIENT APHID EQUATIONS
  dX_t <- (1 - pi)*lambda*psi - # immigration
    (1 - q)*phi_t*a*X_t*v*I/(S + v*I) + # acquisition
    tau_t*Z_t - # infectivity loss
    q*phi_t*X_t # emigration
  
  dZ_t <- pi*lambda*psi + # immigration
    (1 - q)*phi_t*a*X_t*v*I/(S + v*I) - # acquisition
    tau_t*Z_t - # infectivity loss
    q*phi_t*Z_t # emigration
  
  # COLONISING APHID EQUATIONS
  dX_c <- sigma*(X_c + Z_c)*(1 - (X_c + Z_c)/(H*kappa)) + # birth
    (1 - pi)*lambda*(1 - psi) +  # immigration
    tau_c*Z_c - # infectivity loss
    (1 - q_c)*phi_c*X_c*a*(1 - epsilon*omega)*v*I/(S + v*I) - # acquisition
    q_c*phi_c*X_c - # emigration
    theta*X_c # death
  
  dZ_c <- pi*lambda*(1 - psi) + # immigration
    (1 - q_c)*phi_c*X_c*a*(1 - epsilon*omega)*v*I/(S + v*I) - # acquisition
    tau_c*Z_c - # infectivity loss
    q_c*phi_c*Z_c - # emigration
    theta*Z_c # death
  
  
  return(list(c(dI, dX_t, dZ_t, dX_c, dZ_c)))
}

run_combined_ode <- function(times, init_states, parms, death_linked_v = F) {
  run <- data.frame(deSolve::ode(y = init_states, 
                                 times = times,
                                 parms = parms,
                                 func = transient_colonising_ode, 
                                 death_linked_v = death_linked_v))
  S <- parms[["H"]] - run$I
  phi_t <- 1 / (parms[["k"]] + parms[["t"]])
  phi_c <-  (S + parms[["v"]]*run$I) / 
    (parms[["omega"]]*parms[["eta"]]*(S + parms[["v"]]*parms[["epsilon"]]*run$I) +
       (parms[["k"]] + parms[["t"]])*(S + parms[["v"]]*run$I))
  
  run$Xc_em <- parms[["q_c"]]*phi_c*run$X_c
  run$Zc_em <- parms[["q_c"]]*phi_c*run$Z_c
  run$Xt_em <- parms[["q"]]*phi_t*run$X_t
  run$Zt_em <- parms[["q"]]*phi_t*run$Z_t
  run$phi_c <- phi_c
  
  return(run)
}

run_combined_ode_v_vals <- function(v_vals, times, init_states, parms) {
  out <- data.frame()
  for (v in v_vals) {
    parms[["v"]] <- v
    out <- rbind(out, data.frame(run_combined_ode(times, init_states, parms),
                                 v_val = v))
  }
  return(out)
}

calculate_coloni_eqX <- function(parms, include_pi = F) {
  # finds X_c at disease-free equilibrium (coefficients found by setting dX_c/dt=0, see lab book).
  # used in calculate_combi_R0()
  
  coeff0 <- -(1 - parms[["psi"]])*parms[["lambda"]]
  if (include_pi) { coeff0 <- coeff0*(1 - parms[["pi"]])}
  
  coeff1 <- parms[["q_c"]]/(parms[["k"]] + parms[["t"]] + parms[["omega"]]*parms[["eta"]]) + 
    parms[["theta"]] - parms[["sigma"]]
  coeff2 <- parms[["sigma"]]/(parms[["H"]]*parms[["kappa"]])
  eqX <- polyroot(c(coeff0, coeff1, coeff2))
  
  if (!all(round(Im(eqX),6) == 0)) browser()
  eqX <- Re(eqX) # remove imaginary part (is 0)
  
  if (length(eqX) != 2 | sum(eqX > 0) != 1)  {
    if (parms[["psi"]] == 1) {
      
      if (max(eqX) != 0) browser()
      warning("psi=0: DFE Xc = 0, death + emigration exceeds birth")
      
      return(max(eqX)) # will be 0
    }
    else { browser() }
  }
  
  return(eqX[eqX > 0])
}

calculate_inf_loss <- function(parms, combi_run) {
  ## for given parameters and run of the combined (transient + colonising aphid)
  ## model, returns the transient and colonising aphid infectivity loss rates,
  ## tau_t and tau_c, for each value of I
  
  I <- combi_run$I
  S <- parms[["H"]] - I
  
  with(as.list(parms), {
    phi_t <- 1 / (k + t)
    phi_c <- (S + v*I) / 
      (omega*eta*(S + v*epsilon*I) + (t + k)*(S + v*I))
    tau_t <- (1 - q)*phi_t*rho*(S + (1 - a)*v*I) / (S + v*I)
    tau_c <- (1 - q_c)*phi_c*
      (S*(rho + (1 - rho)*omega) + 
         v*I*(rho*(1 - a*(1 - epsilon*omega)) + (1 - rho)*epsilon*omega)) / (S + v*I)
    return(data.frame(t_inf_loss = tau_t, c_inf_loss = tau_c))
  })
}

calculate_transi_acq <- function(parms, I_vals, X_vals) {
  # for given value(s) of I and X and parameters, calculates aphid virus
  # acquisition rate (i.e. attract-influenced acquisition) for transient aphid model
  
  (1 - parms[["q"]]) * (1/(parms[["k"]] + parms[["t"]])) *
    parms[["a"]] * X_vals * parms[["v"]] * I_vals / 
    (parms[["H"]] - I_vals + parms[["v"]] * I_vals)
}

calculate_coloni_acq <- function(parms, I_vals, X_vals) {
  # for given value(s) of I and X and parameters, calculates aphid virus
  # acquisition rate (i.e. attract-influenced acquisition) for transient aphid model
  phi_c <- (parms[["H"]] - I_vals + parms[["v"]] * I_vals) / 
    (parms[["omega"]] * parms[["eta"]] *
       (parms[["H"]] - I_vals + parms[["v"]] * parms[["epsilon"]] * I_vals) + 
       (parms[["t"]] + parms[["k"]]) * (parms[["H"]] - I_vals + parms[["v"]] * I_vals))
  
  return((1 - parms[["q_c"]]) * phi_c * (1 - parms[["epsilon"]]*parms[["omega"]]) *
           parms[["a"]] * X_vals * parms[["v"]] * I_vals / 
           (parms[["H"]] - I_vals + parms[["v"]] * I_vals))
}

calculate_combi_endIXZ <- function(parms, death_linked_v = F) {
  
  run_combined_ode(times = c(0,100000000000),
                   init_states = c(I=1, X_t=1, Z_t=1, X_c=1, Z_c=1),
                   parms = parms,
                   death_linked_v = death_linked_v)[2,-1]
}

calculate_transi_R0 <- function(parms) {
  plant_death <- parms[["gamma"]]*(1 + parms[["mu"]]*(parms[["v"]] - 1)^2)
  inoculation <- parms[["b"]]*(1 / (1 + parms[["mu2"]]*(parms[["v"]] - 1)^2))
  
  return(
    parms[["lambda"]] * (1 - parms[["q"]])^2 * parms[["a"]] * inoculation * parms[["v"]] /
      (parms[["q"]] * ((1 - parms[["q"]])*parms[["rho"]] + parms[["q"]]) * plant_death * parms[["H"]])
  )
}

calculate_combi_R0 <- function(parms, calc_parts = F, include_pi = F) {
  plant_death <- parms[["gamma"]]*(1 + parms[["mu"]]*(parms[["v"]] - 1)^2)
  inoculation <- parms[["b"]]*(1 / (1 + parms[["mu2"]]*(parms[["v"]] - 1)^2))
  
  flight <- parms[["omega"]]*parms[["eta"]] + parms[["t"]] + parms[["k"]]
  
  colonising_R0 <- (1 - parms[["q_c"]])^2 * parms[["a"]] * inoculation * parms[["v"]] *
    (1 - parms[["epsilon"]]*parms[["omega"]]) * calculate_coloni_eqX(parms, include_pi) /
    (flight * parms[["H"]] * plant_death * 
       ((1 - parms[["q_c"]])*(parms[["rho"]] + (1 - parms[["rho"]])*parms[["omega"]]) + 
          parms[["q_c"]] + parms[["theta"]]*flight)
    )
  
  if (calc_parts) {
    return(c("total" = calculate_transi_R0(parms)*parms[["psi"]] + colonising_R0,
             "transient" = calculate_transi_R0(parms)*parms[["psi"]], 
             "colonising" = colonising_R0))
  }
  return(calculate_transi_R0(parms)*parms[["psi"]] + colonising_R0)
}

calculate_R0_over_v <- function(parms, v_vals) {
  v_R0_df <- data.frame()
  for (v in v_vals) {
    parms[["v"]] <- v
    v_R0_df <- rbind(v_R0_df,
                     cbind(t(calculate_combi_R0(parms, calc_parts = T)),
                           v = v))
  }
  return(v_R0_df)
}

find_equilibria_over_v <- function(parms, 
                                   parm_name, 
                                   parm_vals, 
                                   v_vals = seq(1, 10, length.out = 50),
                                   silent = T) {
  
  results_df <- data.frame()
  for (parm_val in parm_vals) {
    
    for (v_val in v_vals) {
      if (!silent) print(paste("v:", v_val, "other parm:", parm_val))
      parms[["v"]] <- v_val
      parms[[parm_name]] <- parm_val
      
      results_df <- rbind(results_df,
                          data.frame(calculate_combi_endIXZ(parms),#[1:5],
                                     v_val = v_val,
                                     parm_val = parm_val,
                                     parm_name = parm_name
                          ))
    }
  }
  return(results_df)
}

set_text_size <- function(plot, size=15) {
  if (inherits(plot, "ggplot")) {
    plot + theme(text = element_text(size = size))
  } else {
    plot  # leave grobs (e.g. legend) unchanged
  }
}


### INVASION FITNESS FUNCTIONS
# for one-year (equilibrium) invasion fitness
mutant_combined_ode <- function(y, parms, resident_eq) {
  
  # mutant part of 2 strain combined model with mutant and resident virus strains (set up for 1 timepoint only),
  # where the values of S, I_r, X_t and X_c are determined by a vector 'resident_eq', and v_m and v_r are
  # mutant and resident values of v respectively. 
  # - returns vector of mutant rates of infection: dI_m, dZt_m, dZc_m
  # - y is the mutant initial conditions: states of I_m, Zt_m, Zc_m
  
  with(as.list(c(y, parms)), {
    S <- H - resident_eq[["I"]]
    X_c <- resident_eq[["X_c"]]
    X_t <- resident_eq[["X_t"]]
    
    # v and I for res and mutant
    v_vals <- c(v_r, v_m)
    I_i <- c(resident_eq[["I"]], I_m)
    
    # assume mutant immigrates proportional to resident equilibrium, like assuming the mutant has evolved in the lanscape +
    # field at the same time
    Zt_last_m <- Zt_m / (Zt_m + resident_eq[["Z_t"]])
    Zc_last_m <- Zc_m / (Zc_m + resident_eq[["Z_c"]])
    
    H_hat <- S + sum(v_vals*I_i)
    
    # aphid flight rates (c = colonising, t = transient)
    phi_t <- 1 / (k + t)
    phi_c <- H_hat / 
      (omega*eta*(S + epsilon*sum(v_vals*I_i)) + (t + k)*H_hat)
    
    # aphid infectivity loss rates (c = colonising, t = transient)
    tau_t <- (1 - q)*phi_t*rho*(S + (1 - a)*sum(v_vals*I_i)) / H_hat
    tau_c <- (1 - q_c)*phi_c*
      (S*(rho + (1 - rho)*omega) + 
         sum(v_vals*I_i)*(rho*(1 - a*(1 - epsilon*omega)) + (1 - rho)*epsilon*omega)) / H_hat
    
    # PLANTS EQUATION
    plant_death_m <- gamma*(1 + mu*(v_m - 1)^2)
    inoculation_m <- b*(1 / (1 + mu2*(v_m - 1)^2))
    dI_m <- ((1 - q)*phi_t*Zt_m + (1 - q_c)*phi_c*Zc_m)*inoculation_m*S/H_hat - plant_death_m*I_m
    
    # TRANSIENT APHID EQUATIONS
    dZt_m <- pi*lambda*psi*Zt_last_m + # immigration = 0
      (1 - q)*phi_t*a*X_t*v_m*I_m / H_hat - # acquisition
      tau_t*Zt_m - # infectivity loss
      q*phi_t*Zt_m # emigration
    
    # COLONISING APHID EQUATIONS
    dZc_m <- pi*lambda*(1 - psi)*Zc_last_m + # immigration = 0
      (1 - q_c)*phi_c*X_c*a*(1 - epsilon*omega)*v_m*I_m / H_hat - # acquisition
      tau_c*Zc_m - # infectivity loss
      q_c*phi_c*Zc_m - # emigration
      theta*Zc_m # death
    
    return(c(dI_m, dZt_m, dZc_m))
  })
}

calculate_invasion_fitness <- function(v_r, v_m, parms, tmax = "equilibrium", silent=T) {
  
  if (!silent) { print(paste("v_r:", v_r, "v_m:", v_m)) }
  
  #if (parms[["pi"]] != 0) { warning("inv fitness: setting pi to 0"); parms[["pi"]] <- 0 }
  if (!is.numeric(tmax) & tmax != "equilibrium") {
    stop("tmax must be 'equilibrium' or a number of days")
  }
  
  parms[["v"]] <- v_r # for calculating resident eq
  parms_mut <- c(parms, v_r = v_r, v_m = v_m)
  
  y <- c(I_m = 1e-7, Zt_m = 1e-7, Zc_m = 1e-7) # mutant initial conditions, very low levels
  
  if (tmax == "equilibrium") {
    resident_eq <- calculate_combi_endIXZ(parms) %>% select(c("I", "X_t", "X_c", "Z_t", "Z_c"))
  }
  else {
    resident_eq <- run_combined_ode(times = c(0, tmax),
                                    init_states = c(I=1, X_t=1, Z_t=1, X_c=1, Z_c=1),
                                    parms = parms)[2,] %>% select(c("I", "X_t", "X_c", "Z_t", "Z_c"))
  }
  # calculate jacobian of mutant model at low mutant level (y) feeding in resident equilibrium values (resident_eq)
  J <- jacobian(func = mutant_combined_ode, 
                x = y, 
                parms = parms_mut, 
                resident_eq = resident_eq, 
                method.args=list(eps = 1e-9)) # dictates resolution of result, smaller=higher
  
  return(max(Re(eigen(J)$values))) # return largest eigenvalue
}

pairwise_invasibility_plot <- function(parms, type = c("one_year", "seasonal"), v_range = c(0, 5), plot_grad=F,
                                       v_lengthout = 100, tmax = "equilibrium", return_plot = T, diagnostic_plots = F) {
  # if return_plot=T, diagnostic_plots=F (default), returns pairwise invasibility plot
  # if return_plot=T, diagnostic_plots=T, returns list of that plot + one with invasion fitness +
  # resident equilibrium I.
  # else will return data
  # see calculate_invasion_fitness() for explanation of tmax
  
  type <- match.arg(type) # one_year or seasonal (multiyear)
  if (type == "seasonal") {
    if (tmax == "equilibrium") stop("type=='seasonal' so set tmax")
    if (plot_grad & v_lengthout > 5) warning(paste("plot_grad=T: will display", v_lengthout^2, "plots, are you sure?"))
    inv_fitness_func <- calculate_seasonal_invasion_fitness
  } else {
    inv_fitness_func <- calculate_invasion_fitness
  }
  
  v_vals <- seq(v_range[1], v_range[2], length.out = v_lengthout)
  
  out <- expand.grid(v_r = v_vals, v_m = v_vals)
  
  if (type == "seasonal") {
    out[, (ncol(out)+1):(ncol(out)+2)] <- 
      t(apply(out, 1, function(row) {
        calculate_seasonal_invasion_fitness(v_r = row[[1]], 
                                            v_m = row[[2]], 
                                            parms = parms, 
                                            tmax = tmax,
                                            plot_grad = plot_grad) }))
    names(out) <- c("v_r", "v_m", "Zt_inv_fitness", "Zc_inv_fitness")
    out$Zt_will_invade <- ifelse(out$Zt_inv_fitness > 0, "yes", "no") 
    out$Zc_will_invade <- ifelse(out$Zc_inv_fitness > 0, "yes", "no") 
    
    if (return_plot) {
      Zt_invas_plot <- ggplot(out, aes(x = v_r, y = v_m)) +
        geom_tile(aes(fill = Zt_will_invade)) +
        geom_contour(aes(z = Zt_inv_fitness), breaks = 0, colour = "black", size = 0.6) +
        scale_fill_manual(values = c("yes" = "grey40", "no" = "white")) +
        coord_equal() +
        geom_hline(yintercept=1, col = "red") +
        geom_vline(xintercept=1, col = "red") +
        labs(title = "Zt invasibility plot")
      Zc_invas_plot <- ggplot(out, aes(x = v_r, y = v_m)) +
        geom_tile(aes(fill = Zc_will_invade)) +
        geom_contour(aes(z = Zc_inv_fitness), breaks = 0, colour = "black", size = 0.6) +
        scale_fill_manual(values = c("yes" = "grey40", "no" = "white")) +
        coord_equal() +
        geom_hline(yintercept=1, col = "red") +
        geom_vline(xintercept=1, col = "red") +
        labs(title = "Zc invasibility plot")
      return(list(Zt_invas_plot, Zc_invas_plot))
    }
  } 
  else {
    out$inv_fitness <- 
      apply(out, 1, function(row) {
        calculate_invasion_fitness(v_r = row[[1]], v_m = row[[2]], parms = parms, tmax = tmax) })
    out$will_invade <- ifelse(out$inv_fitness > 0, "yes", "no") 
    out$res_eq <- sapply(v_vals, function(v) {
      parms[["v"]] <- v
      round(calculate_combi_endIXZ(parms)[["I"]],7) })
    
    if (return_plot) {
      if (diagnostic_plots) {
        invas_plot <- ggplot(out, aes(x = v_r, y = v_m)) +
          geom_tile(aes(fill = will_invade)) +
          geom_contour(aes(z = inv_fitness), breaks = 0, colour = "black", size = 0.6) +
          scale_fill_manual(values = c("yes" = "grey40", "no" = "white")) +
          coord_equal() +
          geom_hline(yintercept=1, col = "red") +
          geom_vline(xintercept=1, col = "red")
        inv_fitness_plot <- ggplot(out, aes(x = v_r, y = v_m)) +
          geom_tile(aes(fill = inv_fitness)) +
          geom_contour(aes(z = inv_fitness), breaks = 0, colour = "black", size = 0.6) +
          coord_equal() +
          geom_hline(yintercept=1, col = "red") +
          geom_vline(xintercept=1, col = "red") +
          labs(fill = "Invasion\nfitness")
        resI_plot <- ggplot(out, aes(x = v_r, y = v_m)) +
          geom_tile(aes(fill = res_eq)) +
          geom_contour(aes(z = inv_fitness), breaks = 0, colour = "black", size = 0.6) +
          coord_equal() +
          geom_hline(yintercept=1, col = "red") +
          geom_vline(xintercept=1, col = "red") +
          labs(fill = "Resident\nequilibrium I")
        return(list(invas_plot, inv_fitness_plot, resI_plot))
      }
      return(ggplot(out, aes(x = v_r, y = v_m)) +
               geom_tile(aes(fill = will_invade)) +
               geom_contour(aes(z = inv_fitness), breaks = 0, colour = "black", size = 0.6) +
               scale_fill_manual(values = c("yes" = "grey40", "no" = "white")) +
               coord_equal() +
               geom_hline(yintercept=1, col = "red") +
               geom_vline(xintercept=1, col = "red"))
    }
  }
  return(out) # else return the data
}

# for seasonal invasion fitness - for reference, not used in this figure
mutant_ode_res_background <- function(time, y, parms, Z_last, I_interp, Xt_interp, Xc_interp, Zt_interp, Zc_interp) {
  with(as.list(c(y, parms)), {
    S <- H - I_interp(time)
    X_c <- Xc_interp(time)
    X_t <- Xt_interp(time)
    
    # v and I for res and mutant
    v_vals <- c(v_r, v_m)
    I_i <- c(I_interp(time), I_m)
    
    # proportion of mutant immigrants given by Z_last - proportion of Z emigration last season
    Zt_last_m <- Z_last[["Zt"]]
    Zc_last_m <- Z_last[["Zc"]]
    
    H_hat <- S + sum(v_vals*I_i)
    
    # aphid flight rates (c = colonising, t = transient)
    phi_t <- 1 / (k + t)
    phi_c <- H_hat / 
      (omega*eta*(S + epsilon*sum(v_vals*I_i)) + (t + k)*H_hat)
    
    # aphid infectivity loss rates (c = colonising, t = transient)
    tau_t <- (1 - q)*phi_t*rho*(S + (1 - a)*sum(v_vals*I_i)) / H_hat
    tau_c <- (1 - q_c)*phi_c*
      (S*(rho + (1 - rho)*omega) + 
         sum(v_vals*I_i)*(rho*(1 - a*(1 - epsilon*omega)) + (1 - rho)*epsilon*omega)) / H_hat
    
    # PLANTS EQUATION
    plant_death_m <- gamma*(1 + mu*(v_m - 1)^2)
    inoculation_m <- b*(1 / (1 + mu2*(v_m - 1)^2))
    dI_m <- ((1 - q)*phi_t*Zt_m + (1 - q_c)*phi_c*Zc_m)*inoculation_m*S/H_hat - plant_death_m*I_m
    
    # TRANSIENT APHID EQUATIONS
    dZt_m <- pi*lambda*psi*Zt_last_m + # immigration
      (1 - q)*phi_t*a*X_t*v_m*I_m / H_hat - # acquisition
      tau_t*Zt_m - # infectivity loss
      q*phi_t*Zt_m # emigration
    
    # COLONISING APHID EQUATIONS
    dZc_m <- pi*lambda*(1 - psi)*Zc_last_m + # immigration 
      (1 - q_c)*phi_c*X_c*a*(1 - epsilon*omega)*v_m*I_m / H_hat - # acquisition
      tau_c*Zc_m - # infectivity loss
      q_c*phi_c*Zc_m - # emigration
      theta*Zc_m # death
    
    return(list(c(dI_m, dZt_m, dZc_m)))
  })
}

run_years_mutant_res_background <- function(v_r, v_m, num_years, tmax, parms) {
  
  v_vals <- c(v_r, v_m)
  num_states <- 2
  out_df <- data.frame(year = rep(1:num_years, each = num_states),
                       state = rep(c("Zt", "Zc"), times = num_years),
                       tmax_val = NA,
                       tmax_prop = NA) # proportion compared to resident value
  
  # initialise parms, times
  times <- seq(0, tmax, length.out = 100)
  
  parms_res <- parms
  parms_res[["v"]] <- v_r
  
  parms_mut <- c(parms, v_r = v_r, v_m = v_m)
  #parms_mut[["pi"]] <- 0 # assume mutant appears in field
  
  # run epidemic of resident for tmax value
  init_states_res <- c(I = 0, # invading scenario
                       X_t = 0,
                       Z_t = 0,
                       X_c = 0,
                       Z_c = 0)
  res_out <- data.frame(ode(y = init_states_res,
                            times = times,
                            parms = parms_res,
                            func = transient_colonising_ode))
  res_endI <- res_out[res_out$time == tmax, "I"]
  res_endZc <- res_out[res_out$time == tmax, "Z_c"]
  res_endZt <- res_out[res_out$time == tmax, "Z_t"]
  
  # create interpolation function for S, Xt, and Xc for the resident strain
  interp_Ifunc <- approxfun(x = res_out$time, y = res_out$I, rule = 2)
  interp_Xtfunc <- approxfun(x = res_out$time, y = res_out$X_t, rule = 2)
  interp_Xcfunc <- approxfun(x = res_out$time, y = res_out$X_c, rule = 2)
  interp_Ztfunc <- approxfun(x = res_out$time, y = res_out$Z_t, rule = 2)
  interp_Zcfunc <- approxfun(x = res_out$time, y = res_out$Z_c, rule = 2)
  
  init_mut_prop <- 1e-5
  # Z_last <- c(Zt = init_mut / (init_mut + Zt_interp(time)),
  #                  Zc = init_mut / (init_mut + Zc_interp(time)))
  Z_last <- c(Zt = init_mut_prop, Zc = init_mut_prop) # so they're the same proportion
  
  for (yr in 1:num_years) {
    run <- data.frame(ode(y = c(I_m = 0, Zt_m = 0, Zc_m = 0),
                          times = seq(0, tmax, length.out = 100),
                          parms = parms_mut,
                          func = mutant_ode_res_background,
                          Z_last = Z_last,
                          I_interp = interp_Ifunc,
                          Xt_interp = interp_Xtfunc,
                          Xc_interp = interp_Xcfunc,
                          Zt_interp = interp_Ztfunc,
                          Zc_interp = interp_Zcfunc))
    
    endZ <- run[nrow(run),-c(1,2)]
    endI <- run[nrow(run), "I_m"]
    I_vals <- c(res_endI, endI)
    H_hat_mut <- parms[["H"]] - res_endI + sum(v_vals*I_vals)
    H_hat_res <- parms[["H"]] - res_endI + v_r*res_endI
    
    phi_mut <- H_hat_mut / (parms[["omega"]]*parms[["eta"]]*(parms[["H"]]-res_endI + parms[["epsilon"]]*sum(v_vals*I_vals)) + 
                              (parms[["t"]] + parms[["k"]])*H_hat_mut)
    phi_res <- H_hat_res / (parms[["omega"]]*parms[["eta"]]*(parms[["H"]]-res_endI + parms[["epsilon"]]*v_r*res_endI) + 
                              (parms[["t"]] + parms[["k"]])*H_hat_res)
    
    Z_last <- c(Zt = endZ[["Zt_m"]] / (endZ[["Zt_m"]] + res_endZt),
                #Zc = endZ[["Zc_m"]] / (endZ[["Zc_m"]] + interp_Zcfunc(tmax))#,
                Zc = phi_mut*endZ[["Zc_m"]] / (phi_mut*endZ[["Zc_m"]] + phi_res*res_endZc)
    )
    out_df[((yr - 1)*num_states + 1):(yr*num_states), "tmax_prop"] <- Z_last
    out_df[((yr - 1)*num_states + 1):(yr*num_states), "tmax_val"] <- unlist(endZ)
    
  }
  return(out_df)
}

calculate_seasonal_invasion_fitness <- function(v_r, v_m, parms, tmax, lower_bound=1, upper_bound=3, plot_grad=F, silent=T) {
  if (!silent) print(paste("v_r:", v_r, "v_m:", v_m))
  
  num_years <- 25
  seasonal_out <- run_years_mutant_res_background(v_r = v_r, 
                                                  v_m = v_m, 
                                                  num_years = num_years, 
                                                  tmax = tmax, 
                                                  parms = parms)
  if (plot_grad) { # plot gradient used as invasion fitness and boundaries
    plot <- ggplot(seasonal_out %>% filter(year < 100), 
                   aes(x = year, y = log(tmax_prop), col = state)) + 
      geom_line() +
      geom_vline(xintercept = lower_bound) +
      geom_vline(xintercept = upper_bound)
    grid.arrange(plot)
  }
  Zc_model <- lm(log(tmax_prop) ~ year, data = (seasonal_out %>% filter(year >= lower_bound &
                                                                          year <= upper_bound &
                                                                          state == "Zc")))
  Zt_model <- lm(log(tmax_prop) ~ year, data = (seasonal_out %>% filter(year >= lower_bound &
                                                                          year <= upper_bound & 
                                                                          state == "Zt")))
  Zc_gradient <- Zc_model$coefficients[["year"]]
  Zt_gradient <- Zt_model$coefficients[["year"]]
  
  return(c(Zt_gradient, Zc_gradient))
  
}


#### PARAMETERS AND INITIAL STATES ####

parms_default <- c(
  k = 0.007, # average probing time (days), equivalent to 10.8 mins
  t = 0.0035, # average time flying between 2 plants (days), equivalent to 5.04 mins
  H = 100, # number of plants
  b = 0.3, # virus inoculation probability
  a = 0.3, # virus acquisition probability
  gamma = 0.02, # infected plant death rate per day
  v = 1, # degree of virus manipulation of aphid landing (1 = no manipulation)
  rho = 0.7, # probability of virus loss from aphid probing
  q = 0.3, # q_t: transient aphid per-flight emigration probability
  lambda = 3, # number of immigrating aphids per day
  pi = 0.05, # proportion of immigrating aphids that are infective
  mu = 0, # degree of cost to manipulation: plant death--linked attractiveness
  mu2 = 0,  # degree of cost to manipulation: reduced inoculation--linked attractiveness
  psi = 0.5, # proportion of immigrants that are transient aphids
  
  # colonising aphids only
  omega = 0.5, # probability of a colonising aphid feeding after probing
  epsilon = 1, # degree of virus manipulation of aphid feeding (1 = no manipulation)
  eta = 2, # length of colonising aphid feed (days)
  q_c = 0.05, # prob of emigration for colonising aphids
  sigma = 0.35, # colonising aphid growth rate (per day) 
  theta = 0.3, # colonising aphid death rate
  kappa = 25 # colonising aphid carrying capacity per plant
)

init_states <- c(I = 0,
                 X_t = 0,
                 Z_t = 0,
                 X_c = 0,
                 Z_c = 0)

v_vals <- c(1, 1.5, 3.5)

#### PANELS A, B, C - TIME COURSES ####
trajecs_v_vals <- run_combined_ode_v_vals(v_vals = v_vals, 
                                          times = 0:900, 
                                          init_states = init_states, 
                                          parms = parms_default)

# calculate acquisition and infectivity loss rates
trajecs_v_vals$C_acq_rate <- calculate_coloni_acq(parms_default, trajecs_v_vals$I, trajecs_v_vals$X_c)
trajecs_v_vals$T_acq_rate <- calculate_transi_acq(parms_default, trajecs_v_vals$I, trajecs_v_vals$X_t)
trajecs_v_vals <- cbind(trajecs_v_vals, calculate_inf_loss(parms_default, trajecs_v_vals))

# PANEL A #
Iplot <- ggplot(trajecs_v_vals, aes(x = time, y = I, col = as.factor(v_val))) +
  geom_line() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Time (days)",
       y = "Number of infected plants (I)")

# PANEL B #
Zplot_w_leg <- Iplot + aes(y = Z_c, lty = "colonising") +
  geom_line(aes(y = Z_t, lty = "transient")) +
  labs(lty = "Aphid type",
       y = expression("Number of infective aphids ("*Z[C]*" or "*Z[T]*")"),
       col = "Infected plant attractiveness (v)",
       title = "(b)") +
  theme(legend.position = "right")

# extract legend
fig3_legend <- cowplot::get_legend(Zplot_w_leg + theme(text = element_text(size = 15)))

Zplot <- Zplot_w_leg + theme(legend.position = "none")

Z_inset <- ggplot(trajecs_v_vals, aes(x = time, y = Z_t, col = as.factor(v_val))) +
  geom_line(lty = 2)  +
  #ylim(0, 0.022) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(legend.position = "none")

combined_Zplot <- cowplot::ggdraw() +
  cowplot::draw_plot(Zplot + theme(text = element_text(size = 15))) +
  cowplot::draw_plot(Z_inset,
                     x = 0.42, y = 0.18,
                     width = 0.5, height = 0.35)

# PANEL C #
Zem_plot <- Iplot + aes(y = Zc_em, lty = "colonising") +
  geom_line(aes(y = Zt_em, lty = "transient")) +
  labs(y = "Rate of infective aphid emigration (/day)")


#### PANELS D, E, F - EQUILIBRIUM VERSUS v ####
equilibrium_over_v <- find_equilibria_over_v(parms = parms_default, 
                                             parm_name = "a", 
                                             parm_vals = parms_default[["a"]],
                                             v_vals = seq(1, 5, length.out = 50))
equilibrium_over_v_points <- find_equilibria_over_v(parms = parms_default, 
                                                    parm_name = "a", 
                                                    parm_vals = parms_default[["a"]],
                                                    v_vals = c(1, 1.5, 3.5))

# PANEL D #
Ieq_plot <- ggplot(equilibrium_over_v, aes(x = v_val, y = I)) +
  geom_line() +
  geom_point(data = equilibrium_over_v_points, aes(col = as.factor(v_val)),
             show.legend = F) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Infected plant attractiveness (v)",
       y = "Equilibrium infected plants (I)") +
  ylim(0, 50)

# PANEL E #
Z_eq_plot <- Ieq_plot + aes(y = Z_c, lty = "colonising") +
  geom_line(aes(y = Z_t, lty = "transient")) +
  geom_point(data = equilibrium_over_v_points, aes(col = as.factor(v_val)),
             show.legend = F) +
  geom_point(data = equilibrium_over_v_points, aes(y = Z_t, col = as.factor(v_val)),
             show.legend = F) +
  labs(y = expression("Equilibrium infective aphids ("*Z[C]*" or "*Z[T]*")"),
       title = "(e)") +
  ylim(0, 10)

Zt_eq_plot <- ggplot(equilibrium_over_v, aes(x = v_val, y = Z_t)) +
  geom_line(lty = 2) +
  geom_point(data = equilibrium_over_v_points, aes(y = Z_t, col = as.factor(v_val)),
             show.legend = F) +
  ylim(0, 0.012) +
  labs(x = NULL, y = NULL) +
  theme_bw()

combined_Zeq_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(Z_eq_plot + theme(text = element_text(size = 15))) +
  cowplot::draw_plot(Zt_eq_plot,
                     x = 0.42, y = 0.18,
                     width = 0.5, height = 0.35)

# PANEL F #
Zem_eq_plot <- Ieq_plot + aes(y = Zc_em, lty = "colonising") +
  geom_line(aes(y = Zt_em, lty = "transient")) +
  geom_point(data = equilibrium_over_v_points, aes(col = as.factor(v_val)),
             show.legend = F) +
  geom_point(data = equilibrium_over_v_points, aes(y = Zt_em, col = as.factor(v_val)),
             show.legend = F) +
  ylim(0, 0.5) +
  labs(y = "Equil. emigrating infective aphids (/day)")


#### PANEL G - v versus R0 ####

v_R0_df <- calculate_R0_over_v(parms_default, seq(0.01, 5, length.out = 100))
v_R0_points_df <- calculate_R0_over_v(parms_default,
                                      v_vals = v_vals)

v_R0_plot <- ggplot(v_R0_df, aes(x = v, y = transient)) +
  geom_line(lty = 2) +
  geom_line(aes(y = colonising)) +
  geom_point(data = v_R0_points_df, aes(col = as.factor(v))) +
  geom_point(data = v_R0_points_df, aes(y = colonising, col = as.factor(v))) +
  labs(x = "Infected plant attractiveness (v)",
       y = expression(""*R[0]*"")) +
  theme_bw() +
  theme(legend.position = "none")

#### PANEL H - PAIRWISE INVASIBILITY PLOT ####
# note: for run time purposes, the resolution parameter v_lengthout 
# has been set to a lower value - the figure in the paper uses v_lengthout=50
pairwise_inv_mu0 <- pairwise_invasibility_plot(parms_default, 
                                               v_range = c(0.001, 4), 
                                               v_lengthout = 20,#50, # figure uses 50
                                               return_plot = T, 
                                               diagnostic_plots = F)
pairwise_inv_mu0 <- pairwise_inv_mu0 + 
  labs(x = expression("Resident strain attractiveness ("*v[r]*")"),
       y = expression("Mutant strain attractiveness ("*v[m]*")"),
       fill = "Mutant can invade?") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(colour = "black"))) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# get legend
pairwise_inv_mu0_leg <- cowplot::get_legend(pairwise_inv_mu0 + 
                                              theme(text = element_text(size = 15)))

pairwise_inv_mu0_no_leg <- pairwise_inv_mu0 + 
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA))
pairwise_inv_mu0_no_leg$layers[[4]] <- NULL # remove red lines at v=1 for final fig
pairwise_inv_mu0_no_leg$layers[[3]] <- NULL


#### COMPILE FIGURE ####
layout_matrix <- matrix(c(1,2,3,
                          1,2,3,
                          4,5,6,
                          4,5,6,
                          7,8,9,
                          7,8,10),
                        byrow=T,nrow=6)

init_fig2_panels <- list(Iplot + ggtitle("(a)"),
                         combined_Zplot,
                         Zem_plot + ggtitle("(c)"),
                         Ieq_plot + ggtitle("(d)"),
                         combined_Zeq_plot,
                         Zem_eq_plot + ggtitle("(f)"),
                         v_R0_plot + ggtitle("(g)"),
                         pairwise_inv_mu0_no_leg + ggtitle("(h)"),
                         fig3_legend,
                         pairwise_inv_mu0_leg)
final_fig2_panels <- lapply(init_fig2_panels, set_text_size) # set text size to 15

pdf("figure3.pdf", height = 12, width = 11)
do.call("grid.arrange", c(final_fig2_panels, 
                          list(layout_matrix = layout_matrix,
                               heights = c(0.5, 0.5, 0.5, 0.5, 0.7, 0.3))))
dev.off()