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
                    recurrent_weights %*% activity)/(arousal + 1)
    activity_history[t, ] <- activity
  }
  return(activity_history)
}

# Generate a fixed true state trajectory
generate_correct_states <- function(n_timepoints, F, initial_state, process_noise_cov) {
  true_states <- matrix(nrow = n_timepoints, ncol = length(initial_state))
  true_states[1, ] <- initial_state
  for (t in 2:n_timepoints) {
    true_states[t, ] <- mvrnorm(n = 1, mu = F %*% true_states[t - 1, ], Sigma = process_noise_cov)
  }
  return(true_states)
}

# Generate noisy input signals
generate_input_signals <- function(true_states, measurement_noise_cov) {
  n_timepoints <- nrow(true_states)
  n_units <- ncol(true_states)
  input_signals <- matrix(0, nrow = n_timepoints, ncol = n_units)
  for (t in 1:n_timepoints) {
    input_signals[t, ] <- mvrnorm(n = 1, mu = true_states[t, ], Sigma = measurement_noise_cov)
  }
  return(input_signals)
}

# Bias–variance simulation function with fixed true state trajectory
simulate_bias_variance_condition_fixed <- function(true_states, n_timepoints, n_runs_bias,
                                                   network_initial, feedforward_weights, recurrent_weights, 
                                                   arousal_value, measurement_noise_cov) {
  bias_results <- list()
  for (r in 1:n_runs_bias) {
    input_signals <- generate_input_signals(true_states, measurement_noise_cov)
    activity_history <- simulate_network(n_timepoints, network_initial, input_signals, 
                                         feedforward_weights, recurrent_weights, arousal_value)
    bias_results[[r]] <- data.frame(Timepoint = 1:n_timepoints,
                                    Activation = activity_history[, 1])
    if (r %% max(1, round(n_runs_bias/10)) == 0) {
      cat(sprintf("Arousal = %.3f: Completed run %d of %d\n", arousal_value, r, n_runs_bias))
    }
  }
  combined_df <- do.call(rbind, bias_results)
  summary_df <- combined_df %>%
    group_by(Timepoint) %>%
    summarize(MeanActivation = mean(Activation),
              SDActivation = sd(Activation),
              .groups = "drop")
  summary_df$Arousal <- arousal_value
  return(summary_df)
}

### COMMON PARAMETERS
n_timepoints <- 5000
n_runs_bias <- 50

# State-transition matrix F for drift diffusion:
F <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)

# Process noise covariance
selected_volatility <- 100
process_noise_cov <- matrix(c(selected_volatility, 0, 0, 0), nrow = 2, byrow = TRUE)

# Measurement noise covariance
selected_noise <- 50000000
measurement_noise_cov <- matrix(c(selected_noise, 0, 0, 0), nrow = 2)

# Feedforward weights for congruent conditions
feedforward_weights <- diag(2)
# For internal model, recurrent weights equal F
recurrent_weights <- F

### CONDITIONS
set.seed(250)
## Condition 1: Diffusion network – Diffusion process (Congruent Diffusion)
true_state_diffusion <- generate_correct_states(n_timepoints, F, c(0, 0), process_noise_cov)
network_initial_diffusion <- c(0, 0)

## For Conditions 2 & 3, true process is drift
drift_value <- 1
# Generate one true state trajectory for the drift process
true_state_drift <- generate_correct_states(n_timepoints, F, c(0, drift_value), process_noise_cov)
# Use the same true state for both Incongruent and Congruent Drift conditions
true_state_incongruent <- true_state_drift  

## Condition 2: Diffusion network – Drift process (Incongruent)
# Network initial state
network_initial_incongruent <- c(0, 0)
# Use feedforward matrix that only updates the position
feedforward_weights_incongruent <- matrix(c(1, 0, 0, 0), nrow = 2, byrow = TRUE)

## Condition 3: Drift network – Drift process (Congruent Drift)
# Network initial state (network internal model matches true drift)
network_initial_drift <- c(0, drift_value)


# Choose two arousal levels
arousal_low <- 0.0001
arousal_high <- .01

### RUN BIAS–VARIANCE SIMULATIONS

# 1. Congruent Diffusion
bias_diffusion_low <- simulate_bias_variance_condition_fixed(true_state_diffusion, n_timepoints, n_runs_bias,
                                                             network_initial_diffusion, feedforward_weights,
                                                             recurrent_weights, arousal_low, measurement_noise_cov)
bias_diffusion_low$Condition <- "Congruent Diffusion (Low Arousal)"
bias_diffusion_high <- simulate_bias_variance_condition_fixed(true_state_diffusion, n_timepoints, n_runs_bias,
                                                              network_initial_diffusion, feedforward_weights,
                                                              recurrent_weights, arousal_high, measurement_noise_cov)
bias_diffusion_high$Condition <- "Congruent Diffusion (High Arousal)"

# 2. Incongruent (Diffusion network tracking Drift process)
# Use modified feedforward weights that only update position
bias_incongruent_low <- simulate_bias_variance_condition_fixed(true_state_incongruent, n_timepoints, n_runs_bias,
                                                               network_initial_incongruent, feedforward_weights_incongruent,
                                                               recurrent_weights, arousal_low, measurement_noise_cov)
bias_incongruent_low$Condition <- "Incongruent (Low Arousal)"
bias_incongruent_high <- simulate_bias_variance_condition_fixed(true_state_incongruent, n_timepoints, n_runs_bias,
                                                                network_initial_incongruent, feedforward_weights_incongruent,
                                                                recurrent_weights, arousal_high, measurement_noise_cov)
bias_incongruent_high$Condition <- "Incongruent (High Arousal)"

# 3. Congruent Drift
bias_drift_low <- simulate_bias_variance_condition_fixed(true_state_drift, n_timepoints, n_runs_bias,
                                                         network_initial_drift, feedforward_weights,
                                                         recurrent_weights, arousal_low, measurement_noise_cov)
bias_drift_low$Condition <- "Congruent Drift (Low Arousal)"
bias_drift_high <- simulate_bias_variance_condition_fixed(true_state_drift, n_timepoints, n_runs_bias,
                                                          network_initial_drift, feedforward_weights,
                                                          recurrent_weights, arousal_high, measurement_noise_cov)
bias_drift_high$Condition <- "Congruent Drift (High Arousal)"

### COMBINE BIAS–VARIANCE DATA FOR PLOTTING
bias_diffusion <- rbind(bias_diffusion_low, bias_diffusion_high)
bias_incongruent <- rbind(bias_incongruent_low, bias_incongruent_high)
bias_drift <- rbind(bias_drift_low, bias_drift_high)

### OBTAIN TRUE STATE TRAJECTORIES (first component)
true_diffusion <- true_state_diffusion[, 1]
true_incongruent <- true_state_incongruent[, 1]
true_drift <- true_state_drift[, 1]

true_state_df_diffusion <- data.frame(Timepoint = 1:n_timepoints, TrueState = true_diffusion, Condition = "Congruent Diffusion")
true_state_df_incongruent <- data.frame(Timepoint = 1:n_timepoints, TrueState = true_incongruent, Condition = "Incongruent")
true_state_df_drift <- data.frame(Timepoint = 1:n_timepoints, TrueState = true_drift, Condition = "Congruent Drift")

### PLOTTING

# Plot for Congruent Diffusion
# Combine the simulation output (bias_diffusion) with the true state
plot_data_diffusion <- rbind(
  bias_diffusion,
  data.frame(Timepoint = 1:n_timepoints, 
             MeanActivation = true_diffusion, 
             SDActivation = 0, 
             Arousal = NA, 
             Condition = "True State")
)

p_diffusion <- ggplot(plot_data_diffusion) +
  geom_ribbon(aes(x = Timepoint, ymin = MeanActivation - SDActivation, ymax = MeanActivation + SDActivation, fill = Condition),
              alpha = 0.5) +
  geom_line(aes(x = Timepoint, y = MeanActivation, color = Condition),
            size = 1) +
  labs(title = "Bias–Variance: Congruent Diffusion", x = "Time", y = "Position") +
  scale_color_manual(values = c("Congruent Diffusion (Low Arousal)" = "steelblue",
                                "Congruent Diffusion (High Arousal)" = "red",
                                "True State" = "black"),
                     name = "Condition") +
  scale_fill_manual(values = c("Congruent Diffusion (Low Arousal)" = "steelblue",
                               "Congruent Diffusion (High Arousal)" = "red",
                               "True State" = "black"),
                    name = "Condition") +
  theme_minimal()

print(p_diffusion)

# Plot for Incongruent.
plot_data_incongruent <- rbind(
  bias_incongruent, 
  data.frame(Timepoint = 1:n_timepoints, 
             MeanActivation = true_incongruent, 
             SDActivation = 0, 
             Arousal = NA, 
             Condition = "True State")
)

p_incongruent <- ggplot(plot_data_incongruent) +
  geom_ribbon(aes(x = Timepoint, ymin = MeanActivation - SDActivation, ymax = MeanActivation + SDActivation, fill = Condition),
              alpha = 0.5) +
  geom_line(aes(x = Timepoint, y = MeanActivation, color = Condition),
            size = 1) +
  labs(title = "Bias–Variance: Incongruent (Diffusion Network on Drift Process)", x = "Time", y = "Position") +
  scale_color_manual(values = c("Incongruent (Low Arousal)" = "steelblue",
                                "Incongruent (High Arousal)" = "red",
                                "True State" = "black"),
                     name = "Condition") +
  scale_fill_manual(values = c("Incongruent (Low Arousal)" = "steelblue",
                               "Incongruent (High Arousal)" = "red",
                               "True State" = "black"),
                    name = "Condition") +
  theme_minimal()

print(p_incongruent)


# Plot for Congruent Drift.
plot_data_drift <- rbind(
  bias_drift, 
  data.frame(Timepoint = 1:n_timepoints, 
             MeanActivation = true_drift, 
             SDActivation = 0, 
             Arousal = NA, 
             Condition = "True State")
)

p_drift <- ggplot(plot_data_drift) +
  geom_ribbon(aes(x = Timepoint, ymin = MeanActivation - SDActivation, ymax = MeanActivation + SDActivation, fill = Condition),
              alpha = 0.5) +
  geom_line(aes(x = Timepoint, y = MeanActivation, color = Condition),
            size = 1) +
  labs(title = "Bias–Variance: Congruent Drift", x = "Time", y = "Position") +
  scale_color_manual(values = c("Congruent Drift (Low Arousal)" = "steelblue",
                                "Congruent Drift (High Arousal)" = "red",
                                "True State" = "black"),
                     name = "Condition") +
  scale_fill_manual(values = c("Congruent Drift (Low Arousal)" = "steelblue",
                               "Congruent Drift (High Arousal)" = "red",
                               "True State" = "black"),
                    name = "Condition") +
  theme_minimal()

print(p_drift)
