library(ggplot2)
library(viridis)
library(dplyr)

#---------------------------
# Simulation Functions
#---------------------------

simulate_network <- function(n_trials, initial_state, input_signals, feedforward_weights, recurrent_weights) {
  state_history <- numeric(n_trials)
  state <- initial_state
  for (t in 1:n_trials) {
    input_signal <- input_signals[t, ]
    net_input <- feedforward_weights %*% input_signal + recurrent_weights[state, ]
    state <- which.max(net_input)
    state_history[t] <- state
  }
  state_history
}

generate_input_signals <- function(n_time_points, correct_states, noise_level) {
  input_signals <- matrix(rnorm(2 * n_time_points, mean = 0, sd = noise_level), 
                          nrow = n_time_points, ncol = 2)
  for (t in 1:n_time_points) {
    # Boost the signal for the correct state 
    input_signals[t, correct_states[t]] <- input_signals[t, correct_states[t]] + 1
  }
  input_signals
}

#---------------------------
# Aggregated Simulation Function
#---------------------------

run_simulation_batch_aggregated <- function(arousal_levels, noise_levels, n_runs, n_time_points,
                                            feedforward_weights, recurrent_weights) {
  aggregated_results <- data.frame()
  for (noise in noise_levels) {
    for (arousal in arousal_levels) {
      scaled_feedforward_weights <- feedforward_weights * arousal
      sum_states <- numeric(n_time_points)
      for (run in 1:n_runs) {
        correct_states <- rep(2, n_time_points)
        input_signals <- generate_input_signals(n_time_points, correct_states, noise)
        state_history <- simulate_network(n_time_points, initial_state = 1,
                                          input_signals, scaled_feedforward_weights, recurrent_weights)
        sum_states <- sum_states + state_history
      }
      batch_df <- data.frame(Time = 1:n_time_points,
                             Arousal = arousal,
                             Noise = noise,
                             SumState = sum_states,
                             Count = n_runs)
      aggregated_results <- rbind(aggregated_results, batch_df)
    }
  }
  aggregated_results
}

#---------------------------
# Function to Combine Aggregated Results
#---------------------------

combine_aggregated_results <- function(agg1, agg2) {
  combined <- full_join(agg1, agg2, by = c("Time", "Arousal", "Noise"), suffix = c(".1", ".2"))
  combined <- combined %>%
    mutate(SumState = coalesce(SumState.1, 0) + coalesce(SumState.2, 0),
           Count = coalesce(Count.1, 0) + coalesce(Count.2, 0),
           MeanState = SumState / Count) %>%
    dplyr::select(Time, Arousal, Noise, SumState, Count, MeanState)
  combined
}

#---------------------------
# Parameters and Weight Matrices
#---------------------------

feedforward_weights <- matrix(c(1, 0,
                                0, 1), nrow = 2, byrow = TRUE)
recurrent_weights <- matrix(c(1, -1,
                              -1, 1), nrow = 2, byrow = TRUE)

# Simulation parameters
arousal_levels <- c(0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1)
noise_levels <- c(5)
n_time_points <- 750

#---------------------------
# Run Simulation in Batches
#---------------------------

n_runs_batch1 <- 50
agg_batch1 <- run_simulation_batch_aggregated(arousal_levels, noise_levels, n_runs_batch1, n_time_points,
                                              feedforward_weights, recurrent_weights)


n_runs_batch2 <- 50
agg_batch2 <- run_simulation_batch_aggregated(arousal_levels, noise_levels, n_runs_batch2, n_time_points,
                                              feedforward_weights, recurrent_weights)

combined_agg <- combine_aggregated_results(agg_batch1, agg_batch2)

#---------------------------
# Plotting
#---------------------------

plot_agg <- ggplot(combined_agg, aes(x = Time, y = MeanState-1, color = factor(Arousal), group = Arousal)) +
  geom_line(size = 1) +
  theme_minimal(base_size = 14) +
  scale_color_viridis(discrete = TRUE) +
  labs(title = "Average Network State Over Time After a Change",
       x = "Time",
       y = "Mean Network State",
       color = "Arousal Level")

print(plot_agg)
