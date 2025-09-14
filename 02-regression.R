library(tidyverse)
library(rstan)
library(rstanarm)

options(mc.cores = parallel::detectCores(logical = FALSE))

d1 <- read_csv('reg_dfs/reg_df_Q12.csv') |> janitor::clean_names()

fit <- stan_glmer(q12 ~ (1 | cnty) + (1 | race) + (1 | edu) + sex + (1 | sex:race) +
                    (1 | edu:age) + (1 | edu:race) + biden_vote_pct + region,
                  family = binomial(link = "logit"),
                  data = d1,
                  prior = normal(0, 1, autoscale = TRUE),
                  prior_covariance = decov(scale = 0.50),
                  adapt_delta = 0.99,
                  refresh = 0,
                  seed = 1010)

print(fit)

poststrat_df <- read_csv('poststrat_df.csv')

epred_mat <- posterior_epred(fit, newdata = poststrat_df, draws = 1000)
mrp_estimates_vector <- epred_mat %*% poststrat_df$count / sum(poststrat_df$count)
mrp_estimate <- c(mean = mean(mrp_estimates_vector), sd = sd(mrp_estimates_vector))
cat("MRP estimate mean, sd [putative]: ", round(mrp_estimate, 3))

counties <- c(1,  4,  7,  9, 10, 12, 13, 15, 16, 19, 20, 21, 24, 28, 30, 31, 33,
              34, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 48, 49, 50, 54, 56, 57)

county_df <- data.frame(
  cnty = counties,
  mrp_estimate = NA,
  mrp_estimate_se = NA
)

for(i in 1:nrow(county_df)) {
  filtering_condition <- which(poststrat_df$cnty == county_df$cnty[i])
  
  # Filtering matrix epred_mat with filtering_condition
  state_epred_mat <- epred_mat[ ,filtering_condition]
  
  # Filtering poststratification table with filtering_condition
  k_filtered <- poststrat_df[filtering_condition, ]$count
  
  # Poststratification step
  mrp_estimates_vector_sub <- state_epred_mat %*% k_filtered / sum(k_filtered)
  
  # MRP estimate for state in row i 
  county_df$mrp_estimate[i] <- mean(mrp_estimates_vector_sub)
  county_df$mrp_estimate_se[i] <- sd(mrp_estimates_vector_sub)
  
}

county_df <- county_df |> mutate(
  mrp_estimate_app = 1 - mrp_estimate
) |> rename(
  mrp_estimate_disapp = mrp_estimate
)

county_key <- read_csv('cnty_num_key.csv')

county_df <- left_join(county_df, county_key, by='cnty')

county_df <- county_df |> mutate(
  net = mrp_estimate_app - mrp_estimate_disapp,
  net = net * 100,
  mrp_estimate_app = mrp_estimate_app * 100,
  mrp_estimate_disapp = mrp_estimate_disapp * 100,
  mrp_estimate_se = mrp_estimate_se * 100
)

write_csv(county_df, 'elem_mrp_biden_app_ca_county.csv')

glimpse(county_df)