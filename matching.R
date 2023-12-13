#load packages
library(DOS2)
library(optmatch)
library(RItools)
library(ggplot2)
library(dplyr)
source('utility.R')
library(knitr)
library(sensitivitymv)
fludata <- read.csv("fludata.txt", sep="")

# Check covariate balance
test = plot(xBalance(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd - 1, data=fludata))

png(filename="matching/balance_pre_any_matching.png")
print(test)
dev.off()
share_received = mean(fludata$receive)
n_receive = count(fludata, receive)

# TO-DO: Other checks of how covariates are distributed for treated and untreated

# Define some functions
plot_discrete = function(data, treat, var_of_interest) {
  treat = sym(treat)
  var_of_interest = sym(var_of_interest)
  
  freqs = data %>% 
    mutate(!!var_of_interest := as.factor(!!var_of_interest)) %>% 
    group_by(!!treat, !!var_of_interest) %>%
    summarize(n = n()) %>%
    mutate(pct = n/sum(n))
  
  ggplot(freqs, aes(x = !!var_of_interest, y = pct, fill = !!treat)) +
    geom_bar(stat = "identity", position = "dodge")  +
    ggtitle(paste0("Relative Frequency of ", var_of_interest, " by treatment status")) +
    labs(fill = "Treatment Status")
}

plot_rel_freq <- function(var, treat, data) {
  
  var = sym(var)
  treat = sym(treat)
  
  ggplot(data, aes(x = !!var, fill = !!treat)) +
    geom_density(alpha = 0.5) +
    ggtitle(paste0("Density of ", var, " by treatment status")) +
    labs(fill = "Treatment Status")
}

# Compute prop scores and covariate distances
fludata$prop <- glm(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd, family=binomial, data=fludata)$fitted.values

# Plot propensity score by treatment status
ggplot(fludata, aes(x = prop, fill = as.factor(receive))) +
  geom_density(alpha = 0.5) +
  ggtitle(paste0("Density of propensity score by treatment status")) +
  labs(fill = "Treatment Status")
ggsave("matching/propensity_score_by_treatment.png")

# Plot propensity score by COPD status
ggplot(fludata, aes(x = prop, fill = as.factor(copd))) +
  geom_density(alpha = 0.5) +
  ggtitle(paste0("Density of propensity score by COPD")) +
  labs(fill = "COPD")
ggsave("matching/propensity_score_by_copd.png")

# Compute covariate distances
fludata= fludata%>% mutate(z = receive)
fludata= fludata%>% mutate(id = row_number())

# Generate 1:1 dataset
dist1 <- match_on(z ~ assign + age + copd + dm + heartd + race + renal + sex + liverd, method = "mahalanobis", data=fludata)
match1 <- pairmatch(dist1, data=fludata)
match1_summary = summarize.match(fludata, match1)

# 1:1 with caliper
dist1_caliper <- addcaliper(dist1, z=fludata$receive, p=fludata$prop, caliper=0.1)
match1_caliper = pairmatch(dist1_caliper, data=fludata)

match1_caliper_summary = summarize.match(fludata, match1_caliper)
match1_caliper_mean_max = mean_and_max(fludata, match1_caliper)
match1_caliper_dist_mean = round(match1_caliper_mean_max[["mean"]],4)
match1_caliper_dist_max = round(match1_caliper_mean_max[["max"]],4)

match1_caliper_diffs = match1_caliper_mean_max[["all"]] %>% select(id.0, id.1, prop.1, prop.0) %>% mutate(diff = prop.1 - prop.0)

# Plot covariate balance before and after 1:1 matching
png(filename="matching/balance_pre_post_1_1_calip.png")
plot(xBalance(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd -1, strata=list(unstrat=NULL, match1_calip=~match1_caliper), data=fludata))
title('Covariate Balance Before and After 1:1 matching with caliper')
legend(
  "topright",
  legend = c("unstrat", "1:1 matching with caliper"),
  inset = 0.01,
  pch = c(15, 16),
  bg = "white"
)
dev.off()

# Generate 1:k dataset
match1_1_k = pairmatch(dist1_caliper, controls = 2, data=fludata)
match1_1_k_summary = summarize.match(fludata, match1_1_k)

# Plot covariate balance before and after 1:k matching
png(filename="matching/balance_pre_post_1_k.png")
plot(xBalance(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd -1, strata=list(unstrat=NULL, match1_1_k=~match1_1_k), data=fludata))
title('Covariate Balance Before and After Matching with k = 2')
legend(
  "topright",
  legend = c("unstrat", "1:k matching"),
  inset = 0.01,
  pch = c(15, 16),
  bg = "white"
)
dev.off()

# Generate new dataset forcing balance on COPD
dist1_k_penalty <- addalmostexact(dist1_caliper, z=fludata$z, f=fludata$copd)
match1_1_k_exact_copd = pairmatch(dist1_k_penalty, controls = 2, data=fludata)
match1_1_k_exact_copd_summary = summarize.match(fludata, match1_1_k_exact_copd)

png(filename="matching/balance_pre_post_1_k_copd_exact.png")
plot(xBalance(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd -1, strata=list(unstrat=NULL, match1_1_k_exact_copd=~match1_1_k_exact_copd), data=fludata))
title('Covariate Balance Before and After 1:k matching with copd exact')
legend(
  "topright",
  legend = c("unstrat", "1:1 matching with copd exact"),
  inset = 0.01,
  pch = c(15, 16),
  bg = "white"
)
dev.off()


##Plot covariate imbalance before and after matching
all = xBalance(receive ~ assign + age + copd + dm + heartd + race + renal + sex + liverd -1,
                strata=list(unstrat=NULL, match1_caliper=~match1_caliper,  match1_1_k=~match1_1_k,
                            match1_1_k_exact_copd=~match1_1_k_exact_copd),
                data=fludata)

plot_diff <- function(data, strata_labels) {
  # Combine list data into a data frame
  data_long <- reshape2::melt(data, id.vars = c("vars", "strata"), variable.name = "stat") %>% filter(stat == "std.diff") %>% rename(diff = value, dataset = strata)
  
  # Plot using ggplot with switched x and y aesthetics
  p <- ggplot(data_long, aes(x = diff, y = vars, shape = dataset, color = dataset)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # Add dashed line at x = 0
    theme_bw() +
    labs(x = "Standard Difference", y = NULL, 
         title = "Standard Difference by Variable and Dataset") +
    scale_shape_manual(values = c(16, 17, 18, 19)) +  # Use different symbols for each stratum
    
    # Adjust legend position
    theme(
      legend.position = c(1, 1),  # Top right
      legend.justification = c(1, 1),  # Inside the frame
      legend.direction = "vertical"
    )
  
  # Print the plot
  print(p)
}

plot_diff(all)
#save plot as a 800 x 800 pixels image 
ggsave(filename="matching/balance_pre_post_all.png", width=1700, height=1200, units = "px")


match1_mean_max = mean_and_max(fludata, match1)
match1_dist_mean = round(match1_mean_max[["mean"]],4)
match1_dist_max = round(match1_mean_max[["max"]],4)

match1_diffs = match1_mean_max[["all"]] %>% select(id.0, id.1, prop.1, prop.0) %>% mutate(diff = prop.1 - prop.0)

# Estimate treatment effect
## Estimate ATT for 1:1, 1:1 caliper, 1:k, 1:k with exact COPD
att = match1_summary %>% mutate(t_i = outcome.1 - outcome.0)  %>% summarize(att = mean(t_i)) %>% pull(att)
att_caliper = match1_caliper_summary %>% mutate(t_i = outcome.1 - outcome.0)  %>% summarize(att = mean(t_i)) %>% pull(att)
att_1_k = match1_1_k_summary %>% mutate(t_i = outcome.1 - outcome.0)  %>% summarize(att = mean(t_i)) %>% pull(att)
att_1_k_exact_copd = match1_1_k_exact_copd_summary %>% mutate(t_i = outcome.1 - outcome.0)  %>% summarize(att = mean(t_i)) %>% pull(att)

# Compute bias-corrected ATT
## Compute different bias-corrected ATTs for different samples
adjustment_1_to_1 = bias_adjusment(data = fludata, match_data = match1_summary)
adjustment_1_to_1_calip = bias_adjusment(data = fludata, match_data = match1_caliper_summary)
adjustment_1_to_k = bias_adjusment(data = fludata, match_data = match1_1_k_summary)
adjustment_1_to_k_exact_copd = bias_adjusment(data = fludata, match_data = match1_1_k_exact_copd_summary)

# Show unadjusted and adjusted ATT for the different datasets
adj_att_1_to_1 = adjustment_1_to_1[["adjusted_att"]]
unadj_att_1_to_1 = adjustment_1_to_1[["unadjusted_att"]]
ci_1_to_1 = adjustment_1_to_1[["ci"]]

adj_att_1_to_1_calip = adjustment_1_to_1_calip[["adjusted_att"]]
unadj_att_1_to_1_calip = adjustment_1_to_1_calip[["unadjusted_att"]]
ci_1_to_1_calip = adjustment_1_to_1_calip[["ci"]]

adj_att_1_to_k = adjustment_1_to_k[["adjusted_att"]]
unadj_att_1_to_k = adjustment_1_to_k[["unadjusted_att"]]
ci_1_to_k = adjustment_1_to_k[["ci"]]

adj_att_1_to_k_exact_copd = adjustment_1_to_k_exact_copd[["adjusted_att"]]
unadj_att_1_to_k_exact_copd = adjustment_1_to_k_exact_copd[["unadjusted_att"]]
ci_1_to_k_exact_copd = adjustment_1_to_k_exact_copd[["ci"]]

# Implement Rosenbaum's sensitivity analysis
## As explained in the documentation, we need to do - (outcome.1 - outcome.0) to test against the alternative that the treatment decreases hospitalizations.

match1_to_1_sens_data = match1_summary %>% mutate(y = -(outcome.1 - outcome.0)) %>% pull(y)
match1_to_1_sens = senmv(match1_to_1_sens_data,gamma=3,trim=1)
match1_to_1_sense_pv = match1_to_1_sens[["pval"]]

match1_to_1_calip_sens_data = match1_caliper_summary %>% mutate(y = -(outcome.1 - outcome.0)) %>% pull(y)
match1_to_1_calip_sens = senmv(match1_to_1_calip_sens_data,gamma=3,trim=1)
match1_to_1_calip_sense_pv = match1_to_1_calip_sens[["pval"]]

match1_to_k_sens_data = match1_1_k_summary %>% mutate(y = -(outcome.1 - outcome.0)) %>% pull(y)
match1_to_k_sens = senmv(match1_to_k_sens_data,gamma=3,trim=1)
match1_to_k_sense_pv = match1_to_k_sens[["pval"]]

match1_to_k_exact_copd_sens_data = match1_1_k_exact_copd_summary %>% mutate(y = -(outcome.1 - outcome.0)) %>% pull(y)
match1_to_k_exact_copd_sens = senmv(match1_to_k_exact_copd_sens_data,gamma=3,trim=1)
match1_to_k_exact_copd_sense_pv = match1_to_k_exact_copd_sens[["pval"]]
