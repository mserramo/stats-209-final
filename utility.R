## puts the results of a pair match in a nice form
## Usage: summarize.match(dataset,pairmatch_output)
summarize.match <- function(dat, ms, ps.name="prop", keep.mset=FALSE) {
  adat <- dat
  adat$mset <- ms
  adat <- adat[!is.na(adat$mset),]
  adat.treat <- adat[adat$z==1, ]
  adat.ctrl <- adat[adat$z==0, ]
  
  adat.m <- merge(adat.treat, adat.ctrl, by="mset", suffixes=c(".1", ".0"))
  
  if(!keep.mset) {
    adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0", "mset"))]
  } else {
    adat.m <- adat.m[, -which(names(adat.m) %in% c("z.1", "z.0"))]        
  }
  adat.m <- adat.m[, sort(names(adat.m), index.return=TRUE)$ix]
  
  p0.name <- paste0(ps.name,".", 0)
  p1.name <- paste0(ps.name,".",1)
  
  adat.m.tmp.1 <- adat.m[, -which(names(adat.m) %in% c(p0.name, p1.name))]
  adat.m.tmp.2 <- adat.m[, c(p0.name, p1.name)]
  
  adat.m <- cbind(adat.m.tmp.1, adat.m.tmp.2)
  
  return(adat.m)
}

#Function to compute mean and max absolute difference in propensity scores.
mean_and_max = function(original, data) {
  
  sum = summarize.match(original, data)
  
  dist_sum = sum %>%
    summarise(prop_dist = abs(prop.1 - prop.0))
  
  mean_dist = mean(dist_sum$prop_dist)
  max_dist = max(dist_sum$prop_dist)
  
  return(list(mean = mean_dist, max = max_dist, all=sum))
}

#Custom-made function to compute the bias-adjusted ATT (note it only works with specific variable names and formulas, not meant to be general)
bias_adjusment = function(data, match_data){
  #Compute bias-corrected estimate of the average treatment effect on the treated.
  #Step 1: Fit a linear regression used to compute mu_0 for those in the control group and mu_1 for the treatment units
  
  reg_only_control = lm(outcome ~ age + copd + dm + heartd + race + renal + sex + liverd,
                        data = data%>% filter(receive==0)) 
  
  reg_only_treat = lm(outcome ~ age + copd + dm + heartd + race + renal + sex + liverd,
                      data = data%>% filter(receive==1))
  
  #Step 2: Compute fitted value for those in the control group using the parameters estimated in the linear regression
  control_sample = data%>% filter(receive==0)
  treat_sample = data%>% filter(receive==1)
  
  control_sample$pred_y = predict(reg_only_control, newdata = control_sample)
  treat_sample$pred_y = predict(reg_only_treat, newdata = treat_sample) 
  
  control_sample = control_sample %>% select(id, pred_y)
  treat_sample = treat_sample %>% select(id, pred_y)
  pred_y = bind_rows(control_sample, treat_sample)
  
  data= data%>% left_join(pred_y, by = "id") 
  
  #Step 3: Compute mu_0 for those in the treatment group based on the parameters estimated in the linear regression
  change_d = data%>% filter(receive == 1) %>% select(-receive) %>% mutate(receive = 0) 
  change_d$mu_0_treat = predict(reg_only_control, newdata = change_d)
  change_d = change_d %>% select(id, mu_0_treat)
  mu_0_for_controls = control_sample %>% rename(mu_0_controls = pred_y)
  
  #Step 4: Attach mu_0 to matches dataset
  match_data_full = match_data %>% 
    left_join(change_d, by = c("id.1" = "id")) %>% 
    left_join(mu_0_for_controls, by = c("id.0" = "id"))
  
  #Step 5: Compute the bias-corrected estimate of the average treatment effect on the treated
  ##Note: our sample is 1:1, so M = 1 and we automatically sum Y over all matches for each i in the current code
  ##Note 2: pred_y.0 is mu_0 for those in the control group, and predicted_mu_0 is mu_0 for those in the treatment group
  match_data_full = match_data_full %>% mutate(b_i = mu_0_treat - mu_0_controls) 
  
  #Step 6: compute mean bias
  hat_b = mean(match_data_full$b_i)
  
  #Step 7: compute unadjusted ATT and apply adjustment
  unadjusted_att = mean(match_data_full$outcome.1 - match_data_full$outcome.0)
  adjusted_att = unadjusted_att - hat_b
  
  #Step 8: compute variance of the bias-corrected estimate of the average treatment effect on the treated
  ##Step 8.1: compute the number of items each control unit has been used as a match for a treated unit
  k = match_data_full %>% 
    group_by(id.0) %>% 
    summarize(k = n()) %>% 
    ungroup() %>% 
    rename(id = id.0)
  
  #Step 8.2: now compute variance
  aux_var = data%>% 
    left_join(k, by = "id") %>%
    mutate(k = ifelse(is.na(k), 0, k)) %>%
    mutate(m = 1) %>% #because we're using a 1:1 sample
    mutate(term = receive * (outcome - pred_y)^2 + (1 - receive) * (k/m)^2 * (outcome - pred_y)^2) 
  
  var_adjusted_att = (sum(aux_var$term) / sum(aux_var$receive == 1)^2)
  
  #Step 9: compute standard error
  se_adjusted_att = sqrt(var_adjusted_att)
  
  #Step 10: compute 95% confidence interval
  ci_lower = adjusted_att - 1.96 * se_adjusted_att
  ci_upper = adjusted_att + 1.96 * se_adjusted_att
  ci = c(ci_lower, ci_upper)
  
  output = list("unadjusted_att" = unadjusted_att, "adjusted_att" = adjusted_att, "se_adjusted_att" = se_adjusted_att, "ci" = ci)
  return(output)
}
