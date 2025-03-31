# ---------------------------
# LOAD LIBRARIES
# ---------------------------
library(ggplot2)
library(MASS)
library(dplyr)

# ---------------------------
# HELPER FUNCTIONS
# ---------------------------

# Simulation update for drift estimation
simulate_network_kalman <- function(n_timepoints, initial_state, input_signals, 
                                    F, arousal, W_ff, H) {
  # F: state-transition matrix used for prediction
  # H: observation matrix that extracts position
  # W_ff: 2x1 feedforward weight vector (applied elementwise to the prediction error)
  # Arousal scales the feedforward correction
  state_history <- matrix(NA, nrow = n_timepoints, ncol = length(initial_state))
  state <- initial_state
  for (t in 1:n_timepoints) {
    s_pred <- F %*% state                  # Predict next state
    z_t <- input_signals[t, ]              # 1D observation (position)
    e_t <- z_t - H %*% s_pred              # Prediction error
    correction <- arousal * (W_ff * as.numeric(e_t))  # Compute correction
    state <- s_pred + correction           # Update state
    state_history[t, ] <- state
  }
  return(state_history)
}

# Generate a true state trajectory using a drift-diffusion model
generate_correct_states <- function(n_timepoints, F, initial_state, process_noise_cov) {
  true_states <- matrix(nrow = n_timepoints, ncol = length(initial_state))
  true_states[1, ] <- initial_state
  for (t in 2:n_timepoints) {
    true_states[t, ] <- mvrnorm(n = 1, mu = F %*% true_states[t - 1, ], Sigma = process_noise_cov)
  }
  return(true_states)
}

# Generate 1D noisy input signals from the true state (using only the position)
generate_input_signals_1d <- function(true_states, measurement_noise_sd) {
  n_timepoints <- nrow(true_states)
  input <- numeric(n_timepoints)
  for (t in 1:n_timepoints) {
    input[t] <- rnorm(1, mean = true_states[t, 1], sd = measurement_noise_sd)
  }
  return(matrix(input, ncol = 1))
}

# ---------------------------
# Simulation Function for Adaptive Drift Error Analysis
# ---------------------------

simulate_adaptive_error <- function(arousal_levels, volatility_levels, n_timepoints, n_runs,
                                    F, initial_true_state, measurement_noise_sd,
                                    network_initial, network_F, W_ff, H) {
  results <- data.frame()
  total_runs <- length(volatility_levels) * length(arousal_levels) * n_runs
  run_counter <- 0
  
  for (vol in volatility_levels) {
    cat(sprintf("Processing volatility: %s\n", vol))
    # Process noise covariance: diagonal with 'vol' in both dimensions
    process_noise_cov <- matrix(c(vol, 0, 0, vol), nrow = 2, byrow = TRUE)
    
    for (alpha in arousal_levels) {
      error_sum <- 0
      for (r in 1:n_runs) {
        # Generate a new true state trajectory for each run
        true_states <- generate_correct_states(n_timepoints, F, initial_true_state, process_noise_cov)
        # Generate noisy 1D observations of position
        input_signals <- generate_input_signals_1d(true_states, measurement_noise_sd)
        # Simulate the network using the Kalman-like update
        est_states <- simulate_network_kalman(n_timepoints, network_initial, input_signals,
                                              network_F, alpha, W_ff, H)
        # Compute mean absolute error in position (first component)
        mae <- mean(abs(true_states[,1] - est_states[,1]))
        error_sum <- error_sum + mae
        run_counter <- run_counter + 1
        
        # Provide progress updates (every 5% of total runs)
        if (run_counter %% max(1, round(total_runs/20)) == 0) {
          cat(sprintf("Completed %d/%d runs\n", run_counter, total_runs))
        }
      }
      avg_error <- error_sum / n_runs
      temp <- data.frame(Arousal = alpha, Volatility = vol, MAE = avg_error, Count = n_runs)
      results <- rbind(results, temp)
    }
  }
  return(results)
}

# ---------------------------
# Batch Combination Function
# ---------------------------

combine_error_batches <- function(agg1, agg2) {
  agg1 <- agg1 %>% mutate(
    Arousal = as.numeric(Arousal),
    Volatility = as.character(Volatility),
    Count = as.numeric(Count)
  )
  agg2 <- agg2 %>% mutate(
    Arousal = as.numeric(Arousal),
    Volatility = as.character(Volatility),
    Count = as.numeric(Count)
  )
  
  combined <- dplyr::bind_rows(agg1, agg2) %>%
    group_by(Arousal, Volatility) %>%
    summarize(
      SumError = sum(MAE * Count, na.rm = TRUE),
      TotalCount = sum(Count, na.rm = TRUE),
      MAE = SumError / TotalCount,
      .groups = "drop"
    ) %>%
    dplyr::select(Arousal, Volatility, MAE, TotalCount)
  
  names(combined)[names(combined) == "TotalCount"] <- "Count"
  return(as.data.frame(combined))
}

# ---------------------------
# PARAMETERS
# ---------------------------

n_timepoints <- 500
n_runs=1000

# Define the state-transition matrix F for drift-diffusion dynamics
F <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)

# Set the initial true state
initial_true_state <- c(0, 0)

# Define measurement noise standard deviation
measurement_noise_sd <- 50000

# Network initial state is [0, 0]
network_initial <- c(0, 0)

# For the adaptive drift network, the internal model matches the true dynamics
network_F <- F

# Define the feedforward weight vector for the adaptive condition (applies correction to both position and drift)
W_ff <- matrix(c(1, 1), nrow = 2, ncol = 1)

# Observation matrix H extracts position
H <- matrix(c(1, 0), nrow = 1)

# Define a range of arousal levels
arousal_levels <- seq(0, 0.002, by = 0.00001)

# Define process noise (volatility) levels
volatility_levels <- c(10, 100, 1000)

# ---------------------------
# RUN SIMULATIONS
# ---------------------------

set.seed(42)
# Run batches
results_batch1 <- simulate_adaptive_error(arousal_levels, volatility_levels, n_timepoints, n_runs,
                                          F, initial_true_state, measurement_noise_sd,
                                          network_initial, network_F, W_ff, H)

results_batch2 <- simulate_adaptive_error(arousal_levels, volatility_levels, n_timepoints, n_runs,
                                          F, initial_true_state, measurement_noise_sd,
                                          network_initial, network_F, W_ff, H)

# Combine batches
combined_results <- combine_error_batches(results_batch1, results_batch1)

# ---------------------------
# PLOTTING RESULTS
# ---------------------------

p_error <- ggplot(combined_results, aes(x = Arousal, y = MAE, color = factor(Volatility))) +
  geom_line(size = 1) +
  labs(title = "Mean Absolute Error vs. Arousal (Adaptive Drift Network)",
       x = "Arousal (Gain)",
       y = "Mean Absolute Error in Position",
       color = "Volatility") +
  ylim(0,30000)+
  xlim(0,.001) +
  theme_minimal()

print(p_error)
