### LOAD LIBRARIES
library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(gridExtra)

### HELPER FUNCTIONS

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
    correction <- arousal * (W_ff * as.numeric(e_t))  # Correction
    state <- s_pred + correction           # Update state
    state_history[t, ] <- state
  }
  return(state_history)
}

# Generate fixed true state trajectory
generate_correct_states <- function(n_timepoints, F, initial_state, process_noise_cov) {
  true_states <- matrix(nrow = n_timepoints, ncol = length(initial_state))
  true_states[1, ] <- initial_state
  for (t in 2:n_timepoints) {
    true_states[t, ] <- mvrnorm(n = 1, 
                                mu = F %*% true_states[t - 1, ], 
                                Sigma = process_noise_cov)
  }
  return(true_states)
}

# Generate 1D input signals from true state (using only the position)
generate_input_signals_1d <- function(true_states, measurement_noise_sd) {
  n_timepoints <- nrow(true_states)
  input <- numeric(n_timepoints)
  for (t in 1:n_timepoints) {
    input[t] <- rnorm(1, mean = true_states[t, 1], sd = measurement_noise_sd)
  }
  return(matrix(input, ncol = 1))
}

# Bias–variance simulation function for drift estimation that records both Position and Drift
simulate_bias_variance_drift_est2 <- function(n_timepoints, n_runs_bias, F, true_states,
                                              measurement_noise_sd,
                                              network_initial, network_F, arousal,
                                              W_ff, H) {
  bias_results <- list()
  for (r in 1:n_runs_bias) {
    input_signals <- generate_input_signals_1d(true_states, measurement_noise_sd)
    activity_history <- simulate_network_kalman(n_timepoints, network_initial, 
                                                input_signals, network_F, arousal, W_ff, H)
    bias_results[[r]] <- data.frame(Timepoint = 1:n_timepoints,
                                    Position = activity_history[, 1],
                                    Drift = activity_history[, 2])
    if (r %% max(1, round(n_runs_bias / 10)) == 0) {
      cat(sprintf("Arousal = %.3f: Completed run %d of %d\n", arousal, r, n_runs_bias))
    }
  }
  combined_df <- do.call(rbind, bias_results)
  
  # Summarize the Position component
  summary_position <- combined_df %>%
    group_by(Timepoint) %>%
    summarize(Mean = mean(Position),
              SD = sd(Position),
              .groups = "drop")
  summary_position$Component <- "Position"
  
  # Summarize the Drift component
  summary_drift <- combined_df %>%
    group_by(Timepoint) %>%
    summarize(Mean = mean(Drift),
              SD = sd(Drift),
              .groups = "drop")
  summary_drift$Component <- "Drift"
  
  summary_df <- rbind(summary_position, summary_drift)
  summary_df$arousal <- arousal
  return(list(summary = summary_df, true_states = true_states))
}

### PARAMETERS FOR DRIFT ESTIMATION BIAS–VARIANCE ANALYSIS

n_timepoints <- 500
n_runs_bias <- 50

# True process: drift diffusion state-transition matrix F
F <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)

# True process initial state: drift = 0
drift_value <- 0
initial_true_state <- c(0, drift_value)

# Process noise covariance (diagonal) applied to both components
selected_volatility <- 100
process_noise_cov <- matrix(c(selected_volatility, 0, 0, selected_volatility), nrow = 2)

# Measurement noise for 1D input
measurement_noise_sd <- 1000

# Define feedforward weight vectors:
# Adaptive condition: use nonzero gain for drift correction
W_ff_adaptive <- matrix(c(1, 1), nrow = 2, ncol = 1)
# Fixed condition: ignore drift (feed 0 for drift correction).
W_ff_fixed    <- matrix(c(1, 0), nrow = 2, ncol = 1)

# Observation matrix H (1x2) that extracts position
H <- matrix(c(1, 0), nrow = 1)

# Network initial state (for both conditions): start at [0, 0]
network_initial <- c(0, 0)

# Define network internal models:
# Adaptive condition: network's internal model equals F
network_F_adaptive <- F
# Fixed condition: force drift to remain 0 by using a model that ignores the off-diagonal (drift) term
network_F_fixed <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)

# Choose arousal value
arousal_value <- .01

set.seed(42)

# Generate one fixed true state trajectory using the true process F
true_states <- generate_correct_states(n_timepoints, F, initial_true_state, process_noise_cov)

### RUN SIMULATIONS

# Condition A: Adaptive Drift Estimation
adaptive_results <- simulate_bias_variance_drift_est2(n_timepoints, n_runs_bias, F, true_states,
                                                      measurement_noise_sd,
                                                      network_initial, network_F_adaptive,
                                                      arousal_value, W_ff_adaptive, H)
adaptive_results$summary$Condition <- "Adaptive Drift"

# Condition B: Fixed Drift (Network ignores drift; drift remains 0)
fixed_results <- simulate_bias_variance_drift_est2(n_timepoints, n_runs_bias, F, true_states,
                                                   measurement_noise_sd,
                                                   network_initial, network_F_fixed,
                                                   arousal_value, W_ff_fixed, H)
fixed_results$summary$Condition <- "Fixed Drift"

# Combine simulation summaries
bias_plot_data <- rbind(adaptive_results$summary, fixed_results$summary)

# Create true state data frames for both Position and Drift
true_state_df_position <- data.frame(
  Timepoint = 1:n_timepoints,
  Mean = true_states[, 1],
  SD = 0,
  Component = "Position",
  Condition = "True State"
)
true_state_df_drift <- data.frame(
  Timepoint = 1:n_timepoints,
  Mean = true_states[, 2],
  SD = 0,
  Component = "Drift",
  Condition = "True State"
)
true_state_df <- rbind(true_state_df_position, true_state_df_drift)

true_state_df$arousal <- NA
true_state_df <- true_state_df[, c("Timepoint", "Mean", "SD", "Component", "arousal", "Condition")]

# Combine simulation data with true state
plot_data_all <- rbind(bias_plot_data, true_state_df)

### PLOTTING

p_bias_drift_est <- ggplot(plot_data_all, aes(x = Timepoint, y = Mean, color = Condition, fill = Condition)) +
  geom_ribbon(aes(ymin = Mean - SD, ymax = Mean + SD), alpha = 0.3, color = NA) +
  geom_line(size = 1) +
  facet_wrap(~ Component, scales = "free_y", ncol = 1) +
  labs(title = "Bias–Variance Analysis: Drift Estimation",
       x = "Time", y = "Activation",
       color = "Condition", fill = "Condition") +
  scale_color_manual(values = c("Adaptive Drift" = "red",
                                "Fixed Drift" = "blue",
                                "True State" = "black")) +
  scale_fill_manual(values = c("Adaptive Drift" = "red",
                               "Fixed Drift" = "blue",
                               "True State" = "black")) +
  theme_minimal()

print(p_bias_drift_est)

