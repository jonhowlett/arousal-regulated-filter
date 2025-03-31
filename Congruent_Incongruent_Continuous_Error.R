### LOAD LIBRARIES
library(ggplot2)
library(MASS)
library(dplyr)

### HELPER FUNCTIONS

# Simulate continuous network dynamics
simulate_network <- function(n_timepoints, initial_activity, input_signals, 
                             feedforward_weights, recurrent_weights, arousal) {
  activity_history <- matrix(nrow = n_timepoints, ncol = length(initial_activity))
  activity <- initial_activity
  for (t in 1:n_timepoints) {
    activity <- (arousal * feedforward_weights %*% input_signals[t, ] + 
                   recurrent_weights %*% activity) / (arousal + 1)
    activity_history[t, ] <- activity
  }
  return(activity_history)
}

# Generate true states using state-transition model plus process noise
generate_correct_states <- function(n_timepoints, state_transition, 
                                    initial_true_state, process_noise_cov) {
  correct_states <- matrix(nrow = n_timepoints, ncol = length(initial_true_state))
  correct_states[1, ] <- initial_true_state
  for (t in 2:n_timepoints) {
    correct_states[t, ] <- mvrnorm(n = 1, mu = state_transition %*% correct_states[t - 1, ], 
                                   Sigma = process_noise_cov)
  }
  return(correct_states)
}

# Generate noisy input signals given the true states and measurement noise
generate_input_signals <- function(correct_states, measurement_noise_cov) {
  n_timepoints <- nrow(correct_states)
  n_units <- ncol(correct_states)
  input_signals <- matrix(0, nrow = n_timepoints, ncol = n_units)
  for (t in 1:n_timepoints) {
    input_signals[t, ] <- mvrnorm(n = 1, mu = correct_states[t, ], Sigma = measurement_noise_cov)
  }
  return(input_signals)
}

# Simulation batch function for error analysis
simulate_continuous_error <- function(arousal_levels, noise_levels, volatility_levels, 
                                      n_timepoints, n_runs, network_type,
                                      state_transition, initial_true_state, network_initial,
                                      feedforward_weights, recurrent_weights) {
  error_results <- data.frame()
  
  for (noise_level in noise_levels) {
    for (volatility in volatility_levels) {
      cat(sprintf("Processing: Noise Level = %s, Volatility = %s\n", noise_level, volatility))
      error_means <- numeric(length(arousal_levels))
      
      for (i in seq_along(arousal_levels)) {
        arousal <- arousal_levels[i]
        error_sum <- 0
        
        for (r in 1:n_runs) {
          # Process noise: 
          process_noise_cov <- matrix(c(volatility, 0, 0, 0), nrow = 2, byrow = TRUE)
          # Generate true state trajectory for this run
          correct_states <- generate_correct_states(n_timepoints, state_transition, 
                                                    initial_true_state, process_noise_cov)
          # Measurement noise: 
          measurement_noise_cov <- matrix(c(noise_level, 0, 0, 0), nrow = 2, byrow = TRUE)
          input_signals <- generate_input_signals(correct_states, measurement_noise_cov)
          
          # Run network simulation with network_initial (which differs by condition)
          activity_history <- simulate_network(n_timepoints, network_initial, input_signals, 
                                               feedforward_weights, recurrent_weights, arousal)
          
          # Compute error: mean absolute difference on the position (first component)
          error_sum <- error_sum + mean(abs(correct_states[, 1] - activity_history[, 1]))
        }
        error_means[i] <- error_sum / n_runs
        
        if (i %% max(1, round(length(arousal_levels) / 10)) == 0) {
          cat(sprintf("  Completed %d / %d arousal levels\n", i, length(arousal_levels)))
        }
      }
      temp_df <- data.frame(Arousal = arousal_levels,
                            ErrorMean = error_means,
                            Noise = noise_level,
                            Volatility = volatility,
                            NetworkType = network_type,
                            Count = n_runs)
      error_results <- rbind(error_results, temp_df)
    }
  }
  return(error_results)
}

# Function to combine batches (weighted average by count)
combine_error_batches <- function(agg1, agg2) {
  agg1 <- agg1 %>% mutate(
    Arousal = as.numeric(Arousal),
    Noise = as.character(Noise),
    Volatility = as.character(Volatility),
    NetworkType = as.character(NetworkType),
    Count = as.numeric(Count)
  )
  agg2 <- agg2 %>% mutate(
    Arousal = as.numeric(Arousal),
    Noise = as.character(Noise),
    Volatility = as.character(Volatility),
    NetworkType = as.character(NetworkType),
    Count = as.numeric(Count)
  )
  
  combined <- dplyr::bind_rows(agg1, agg2) %>%
    group_by(Arousal, Noise, Volatility, NetworkType) %>%
    summarize(
      SumError = sum(ErrorMean * Count, na.rm = TRUE),
      TotalCount = sum(Count, na.rm = TRUE),
      ErrorMean = SumError / TotalCount,
      .groups = "drop"
    ) %>%
    dplyr::select(Arousal, Noise, Volatility, NetworkType, ErrorMean, TotalCount)
  
  names(combined)[names(combined) == "TotalCount"] <- "Count"
  return(as.data.frame(combined))
}

### PARAMETERS

n_timepoints <- 500
n_runs <- 10

# One volatility and one noise value
volatility_levels <- c(1000)  
noise_levels <- c(50000)      

# Define arousal levels
arousal_levels <- seq(0, 5, by = .01)

# State-transition matrix F for drift diffusion
F <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)

# For drift process, true state is generated with drift
drift_value <- 100
# Both conditions use the same true state trajectory (generated per run):
initial_true_state <- c(0, drift_value)

# For the **Congruent** condition (network internal model matches true drift):
network_initial_congruent <- c(0, drift_value)
feedforward_weights_congruent <- diag(2)  # full update of both components
recurrent_weights_congruent <- F

# For the **Incongruent** condition (network does not incorporate drift):
network_initial_incongruent <- c(0, 0)   
# Feedforward matrix only updates the position
feedforward_weights_incongruent <- matrix(c(1, 0, 0, 0), nrow = 2, byrow = TRUE)
recurrent_weights_incongruent <- F  # Recurrent dynamics remain the same

### RUN SIMULATIONS FOR ERROR ANALYSIS

#Run Congruent condition (network is congruent with drift process)
error_results_congruent1 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels,
                                                     n_timepoints, n_runs, "Congruent",
                                                     F, initial_true_state, network_initial_congruent,
                                                     feedforward_weights_congruent, recurrent_weights_congruent)
# Run Incongruent condition (network ignores drift)
error_results_incongruent1 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels,
                                                       n_timepoints, n_runs, "Incongruent",
                                                       F, initial_true_state, network_initial_incongruent,
                                                       feedforward_weights_incongruent, recurrent_weights_incongruent)

# Combine the results
error_results1 <- rbind(error_results_congruent1, error_results_incongruent1)



# Run Congruent condition (network is congruent with drift process)
error_results_congruent2 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels,
                                                     n_timepoints, n_runs, "Congruent",
                                                     F, initial_true_state, network_initial_congruent,
                                                     feedforward_weights_congruent, recurrent_weights_congruent)
# Run Incongruent condition (network ignores drift)
error_results_incongruent2 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels,
                                                       n_timepoints, n_runs, "Incongruent",
                                                       F, initial_true_state, network_initial_incongruent,
                                                       feedforward_weights_incongruent, recurrent_weights_incongruent)

# Combine the results
error_results2 <- rbind(error_results_congruent2, error_results_incongruent2)

error_results=combine_error_batches(error_results1,error_results2)

### PLOTTING

# Plot mean absolute error vs. arousal, color by NetworkType.
p_error <- ggplot(error_results, aes(x = Arousal, y = ErrorMean, color = NetworkType)) +
  geom_line(size = 1) +
  labs(title = "Mean Absolute Error vs. Arousal", 
       x = "Arousal", 
       y = "Mean Absolute Error", 
       color = "Network Type") +
  ylim(0,200) +
  theme_minimal()

print(p_error)
