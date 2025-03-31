library(ggplot2)
library(viridis)
library(dplyr)

# ---- Helper Functions ----

# Map state index to activity vector for each network type
state_index_to_activity_vector <- function(state_index, network_type) {
  if (network_type == "TwoPop") {  # Congruent
    if (state_index == 1) {
      return(c(rep(1, 6), rep(0, 6)))
    } else {
      return(c(rep(0, 6), rep(1, 6)))
    }
  } else if (network_type == "TwoPopIncong") {  # Incongruent
    if (state_index == 1) {
      return(c(rep(1, 3), rep(0, 3), rep(1, 3), rep(0, 3)))
    } else {
      return(c(rep(0, 3), rep(1, 3), rep(0, 3), rep(1, 3)))
    }
  }
}

# Return the number of states
num_states <- function(network_type) {
  if (network_type %in% c("TwoPop", "TwoPopIncong")) {
    return(2)
  }
}

# Simulate network dynamics
simulate_network <- function(n_timepoints, initial_activity, input_signals, feedforward_weights, recurrent_weights) {
  activity_history <- matrix(nrow = n_timepoints, ncol = length(initial_activity))
  activity <- initial_activity
  for (t in 1:n_timepoints) {
    # Threshold the sum of feedforward and recurrent inputs at 0.5
    activity <- as.numeric((feedforward_weights %*% input_signals[t, ] + recurrent_weights %*% activity) > 0.5)
    activity_history[t, ] <- activity
  }
  return(activity_history)
}

# Compute mean Hamming distance error (averaged over time)
calculate_mean_hamming_distance <- function(activity_history, correct_states, network_type) {
  n_timepoints <- nrow(activity_history)
  hamming_sum <- 0
  for (t in 1:n_timepoints) {
    activity <- activity_history[t, ]
    correct_activity <- state_index_to_activity_vector(correct_states[t], network_type)
    hamming_sum <- hamming_sum + sum(abs(activity - correct_activity))
  }
  return(hamming_sum / n_timepoints)
}

# Generate correct state sequence based on hazard rate
generate_correct_states <- function(n_timepoints, hazard_rate, n_states, initial_correct_state) {
  correct_states <- numeric(n_timepoints)
  correct_states[1] <- initial_correct_state
  for (t in 2:n_timepoints) {
    if (runif(1) < hazard_rate) {
      possible_states <- setdiff(1:n_states, correct_states[t - 1])
      correct_states[t] <- sample(possible_states, 1)
    } else {
      correct_states[t] <- correct_states[t - 1]
    }
  }
  return(correct_states)
}

# Generate input signals using the true state sequence
generate_input_signals <- function(n_timepoints, correct_states, noise_level, network_type, n_input_units) {
  input_signals <- matrix(0, nrow = n_timepoints, ncol = n_input_units)
  for (t in 1:n_timepoints) {
    input_mean <- state_index_to_activity_vector(correct_states[t], network_type)
    input_signals[t, ] <- input_mean + rnorm(n_input_units, mean = 0, sd = noise_level)
  }
  return(input_signals)
}

# Change internal network_type label
network_label <- function(network_type) {
  if (network_type == "TwoPop") {
    return("Congruent")
  } else if (network_type == "TwoPopIncong") {
    return("Incongruent")
  } else {
    return(network_type)
  }
}

# ---- Batch Simulation Aggregation Functions ----

# Run one batch of simulations for congruent/incongruent networks
run_simulation_batch_congruence <- function(arousal_levels, noise_levels, hazard_rates, n_timepoints, n_runs, network_types) {
  batch_results <- data.frame()
  
  for (network_type in network_types) {
    # Set up weight matrices
    # For both types, use a 12x12 feedforward matrix (identity) and recurrent weights with:
    # intra-population excitation and inter-population inhibition
    feedforward_weights <- diag(12)
    recurrent_weights <- matrix(1, nrow = 12, ncol = 12)
    recurrent_weights[1:6, 7:12] <- -1
    recurrent_weights[7:12, 1:6] <- -1
    
    for (noise_level in noise_levels) {
      cat(sprintf("Processing %s - Noise Level %s\n", network_type, noise_level))
      for (hazard_rate in hazard_rates) {
        cat(sprintf("  Hazard Rate %s\n", hazard_rate))
        hamming_means <- numeric(length(arousal_levels))
        
        for (i in seq_along(arousal_levels)) {
          arousal <- arousal_levels[i]
          scaled_feedforward_weights <- feedforward_weights * arousal
          hamming_sum <- 0
          
          for (r in 1:n_runs) {
            true_states <- generate_correct_states(n_timepoints, hazard_rate, 
                                                   num_states(network_type), 
                                                   sample(1:num_states(network_type), 1))
            input_signals <- generate_input_signals(n_timepoints, true_states, noise_level, network_type, 12)
            randomized_initial_activity <- sample(0:1, 12, replace = TRUE)
            activity_history <- simulate_network(n_timepoints, randomized_initial_activity, 
                                                 input_signals, scaled_feedforward_weights, recurrent_weights)
            hamming_sum <- hamming_sum + calculate_mean_hamming_distance(activity_history, true_states, network_type)
          }
          hamming_means[i] <- hamming_sum / n_runs
        }
        temp_results <- data.frame(Arousal = arousal_levels,
                                   HammingMean = hamming_means,
                                   Noise = paste("Noise Level", formatC(noise_level, format = "f")),
                                   Hazard = paste("Hazard Rate", formatC(hazard_rate, format = "f")),
                                   NetworkType = network_label(network_type),
                                   Count = n_runs)
        batch_results <- rbind(batch_results, temp_results)
      }
    }
  }
  return(batch_results)
}

# Combine two aggregated batches using weighted average based on count
combine_aggregated_results <- function(agg1, agg2) {
  agg1 <- agg1 %>% mutate(
    Arousal = as.numeric(Arousal),
    Noise = as.character(Noise),
    Hazard = as.character(Hazard),
    NetworkType = as.character(NetworkType),
    Count = as.numeric(Count)
  )
  agg2 <- agg2 %>% mutate(
    Arousal = as.numeric(Arousal),
    Noise = as.character(Noise),
    Hazard = as.character(Hazard),
    NetworkType = as.character(NetworkType),
    Count = as.numeric(Count)
  )
  
  combined <- dplyr::bind_rows(agg1, agg2) %>%
    group_by(Arousal, Noise, Hazard, NetworkType) %>%
    summarize(
      HammingMean = sum(HammingMean * Count, na.rm = TRUE) / sum(Count, na.rm = TRUE),
      Count = sum(Count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::select(Arousal, Noise, Hazard, NetworkType, HammingMean, Count)
  
  return(as.data.frame(combined))
}


# ---- Parameters ----

arousal_levels <- seq(0, 7.5, by = 0.05)
n_timepoints <- 500
noise_levels <- 1
n_runs <- 20000
hazard_rates <- 0
network_types <- c("TwoPop", "TwoPopIncong")

# ---- Run Batch Simulations ----

batch_next <- run_simulation_batch_congruence(arousal_levels, noise_levels, hazard_rates, n_timepoints, n_runs, network_types)

if(exists("results") && nrow(results) > 0) {
  results <- combine_aggregated_results(results, batch_next)
} else {
  results <- batch_next
}

# ---- Plotting ----

# Facet by Noise and Hazard, color by NetworkType
ggplot(results, aes(x = Arousal, y = HammingMean, color = NetworkType)) +
  geom_line() +
  facet_wrap(~ Noise + Hazard, scales = "fixed") +
  scale_color_manual(values = c("Congruent" = "steelblue", "Incongruent" = "firebrick"), name = "Network Type") +
  labs(title = "Simulation Results",
       x = "Arousal",
       y = "Mean Hamming Distance") +
  theme_minimal()

