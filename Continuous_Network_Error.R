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
    # Update rule
    activity <- (arousal * feedforward_weights %*% input_signals[t, ] + 
                   recurrent_weights %*% activity) / (arousal + 1)
    activity_history[t, ] <- activity
  }
  return(activity_history)
}

# Generate true states using state-transition model plus process noise
generate_correct_states <- function(n_timepoints, state_transition, 
                                    initial_correct_state, process_noise_cov) {
  correct_states <- matrix(nrow = n_timepoints, ncol = length(initial_correct_state))
  correct_states[1, ] <- initial_correct_state
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

### SIMULATION BATCH FUNCTION FOR ERROR ANALYSIS

simulate_continuous_error <- function(arousal_levels, noise_levels, volatility_levels, 
                                      n_timepoints, n_runs, network_type = "Diffusion",
                                      state_transition, initial_correct_state, 
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
          # Process noise covariance:
          process_noise_cov <- matrix(c(volatility, 0, 0, 0), nrow = 2)
          correct_states <- generate_correct_states(n_timepoints, state_transition,
                                                    initial_correct_state, process_noise_cov)
          # Measurement noise covariance:
          measurement_noise_cov <- matrix(c(noise_level, 0, 0, 0), nrow = 2)
          input_signals <- generate_input_signals(correct_states, measurement_noise_cov)
          
          # For congruent network, initial activity is set to the initial true state
          activity_history <- simulate_network(n_timepoints, initial_correct_state, 
                                               input_signals, feedforward_weights, recurrent_weights, 
                                               arousal)
          # Compute error as the mean absolute difference on the first state
          error_sum <- error_sum + mean(abs(correct_states[, 1] - activity_history[, 1]))
        }
        error_means[i] <- error_sum / n_runs
        
        # Provide a progress update every 10% of the arousal levels
        if (i %% max(1, round(length(arousal_levels) / 10)) == 0) {
          cat(sprintf("  Completed %d / %d arousal levels\n", i, length(arousal_levels)))
        }
      }
      temp_df <- data.frame(Arousal = arousal_levels,
                            ErrorMean = error_means,
                            Noise = noise_level,
                            Volatility = volatility,
                            NetworkType = ifelse(network_type == "Diffusion", "Congruent", "Incongruent"),
                            Count = n_runs)
      error_results <- rbind(error_results, temp_df)
    }
  }
  return(error_results)
}

### FUNCTION TO COMBINE BATCHES
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


### ERROR ANALYSIS

# Parameters for the congruent (diffusion) case
state_transition <- matrix(c(1, 0, 1, 1), nrow = 2, byrow = TRUE)
initial_correct_state <- c(0, 0)  # For diffusion, drift is zero
feedforward_weights <- diag(2)
recurrent_weights <- state_transition

# Define parameters
arousal_levels <- seq(0, 0.2, by = 0.005)
n_timepoints <- 500
noise_levels <- seq(0, 80000, by=20000)
volatility_levels <- c(0, 10, 100, 1000)
n_runs_batch1 <- 50  
n_runs_batch2 <- 50  

# Run batches
batch1 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels, 
                                    n_timepoints, n_runs_batch1, network_type = "Diffusion",
                                    state_transition, initial_correct_state, 
                                    feedforward_weights, recurrent_weights)

batch2 <- simulate_continuous_error(arousal_levels, noise_levels, volatility_levels, 
                                    n_timepoints, n_runs_batch2, network_type = "Diffusion",
                                    state_transition, initial_correct_state, 
                                    feedforward_weights, recurrent_weights)

combined_results <- combine_error_batches(batch1, batch2)

error_results = combined_results

# Plot Error Results:
p_error <- ggplot(error_results, aes(x = Arousal, y = ErrorMean, color = factor(Noise))) +
  geom_line(size = 1) +
  facet_wrap(~ Volatility, scales = "fixed") +
  labs(title = "Mean Absolute Error vs. Arousal (Congruent Diffusion)",
       x = "Arousal",
       y = "Mean Absolute Error",
       color = "Noise Level") +
  ylim(0,100)+
  theme_minimal()
print(p_error)
