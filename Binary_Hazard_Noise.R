library(ggplot2)
library(viridis)
library(dplyr)

# ---------------------------
# Simulation & Metric Functions
# ---------------------------

simulate_network <- function(n_trials, init_state, inputs, W_ff, W_rec) {
  state_hist <- numeric(n_trials)
  state <- init_state
  for (t in 1:n_trials) {
    net_input <- W_ff %*% inputs[t, ] + W_rec[state, ]
    state <- which.max(net_input)
    state_hist[t] <- state
  }
  return(state_hist)
}

calculate_switch_frequency <- function(state_hist) {
  return(sum(diff(state_hist) != 0) / length(state_hist))
}

calculate_correct_state_percent <- function(state_hist, true_states) {
  return(sum(state_hist == true_states) / length(state_hist))
}

# ---------------------------
# True State and Input Generators
# ---------------------------

generate_correct_states <- function(n_pts, hazard_rate) {
  states <- numeric(n_pts)
  states[1] <- 1
  for (t in 2:n_pts) {
    states[t] <- ifelse(runif(1) < hazard_rate, 3 - states[t - 1], states[t - 1])
  }
  return(states)
}

generate_input_signals <- function(n_pts, true_states, noise_level) {
  inputs <- matrix(rnorm(2 * n_pts, mean = 0, sd = noise_level), nrow = n_pts, ncol = 2)
  for (t in 1:n_pts) {
    inputs[t, true_states[t]] <- inputs[t, true_states[t]] + 1
  }
  return(inputs)
}

# ---------------------------
# Weight Matrices
# ---------------------------

W_ff <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
W_rec <- matrix(c(1, -1, -1, 1), nrow = 2, byrow = TRUE)

# ---------------------------
# Batch Aggregation Functions with Reduced Progress Updates
# ---------------------------

run_simulation_batch_aggregated <- function(arousal_levels, noise_levels, hazard_rates, n_pts, n_runs, W_ff, W_rec) {
  batch_results <- data.frame()
  
  for (noise in noise_levels) {
    for (hazard in hazard_rates) {
      # Print one progress update for each outer combination
      cat(sprintf("Processing: Noise = %s, Hazard = %s\n", noise, hazard))
      
      for (a in arousal_levels) {
        W_ff_scaled <- W_ff * a
        sum_switch <- 0
        sum_correct <- 0
        
        for (r in 1:n_runs) {
          true_states <- generate_correct_states(n_pts, hazard)
          inputs <- generate_input_signals(n_pts, true_states, noise)
          init_state <- sample(1:2, 1)
          state_hist <- simulate_network(n_pts, init_state, inputs, W_ff_scaled, W_rec)
          sum_switch <- sum_switch + calculate_switch_frequency(state_hist)
          sum_correct <- sum_correct + calculate_correct_state_percent(state_hist, true_states)
        }
        temp <- data.frame(
          Arousal = a,
          Noise = paste("Noise Level", formatC(noise, format = "f")),
          Hazard = paste("Hazard Rate", formatC(hazard, format = "f")),
          SumSwitch = sum_switch,
          SumCorrect = sum_correct,
          Count = n_runs
        )
        batch_results <- rbind(batch_results, temp)
      }
    }
  }
  return(batch_results)
}

combine_aggregated_results <- function(agg1, agg2) {
  combined <- full_join(agg1, agg2, by = c("Arousal", "Noise", "Hazard"), suffix = c(".1", ".2"))
  combined <- combined %>%
    mutate(
      SumSwitch = coalesce(SumSwitch.1, 0) + coalesce(SumSwitch.2, 0),
      SumCorrect = coalesce(SumCorrect.1, 0) + coalesce(SumCorrect.2, 0),
      Count = coalesce(Count.1, 0) + coalesce(Count.2, 0),
      SwitchFrequency = SumSwitch / Count,
      CorrectStatePercent = SumCorrect / Count
    ) %>%
    dplyr::select(Arousal, Noise, Hazard, SumSwitch, SumCorrect, Count, SwitchFrequency, CorrectStatePercent)
  return(combined)
}

# ---------------------------
# Parameters
# ---------------------------

arousal_levels <- seq(0, 0.75, by = 0.01)
n_pts <- 500
noise_levels <- seq(1, 5, by = 1)
hazard_rates <- seq(0, 0.05, by = 0.01)
n_runs_batch1 <- 50
n_runs_batch2 <- 50

# ---------------------------
# Run Batch Simulations & Combine Batches
# ---------------------------

batch1 <- run_simulation_batch_aggregated(arousal_levels, noise_levels, hazard_rates, n_pts, n_runs_batch1, W_ff, W_rec)

batch2 <- run_simulation_batch_aggregated(arousal_levels, noise_levels, hazard_rates, n_pts, n_runs_batch2, W_ff, W_rec)

combined_results <- combine_aggregated_results(batch1, batch2)

# ---------------------------
# Plotting Aggregated Results
# ---------------------------

switch_frequency_plot <- ggplot(combined_results[combined_results$Hazard=='Hazard Rate 0.0000',], aes(x = Arousal, y = SwitchFrequency, color = Noise)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Switch Frequency vs. Arousal Level",
       x = "Arousal Level", y = "Switch Frequency") +
  scale_color_viridis(discrete = TRUE, name = "Noise Level")


correct_state_percent_plot <- ggplot(combined_results, aes(x = Arousal, y = CorrectStatePercent, color = Noise)) +
  geom_line() +
  facet_wrap(~ Hazard) +
  theme_minimal() +
  labs(title = "Correct State Percent vs. Arousal Level",
       x = "Arousal Level", y = "Correct State Percent") +
  scale_color_viridis(discrete = TRUE, name = "Noise Level")

print(switch_frequency_plot)
print(correct_state_percent_plot)
