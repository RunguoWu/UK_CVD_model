
# Functions to be called --------------------------------------------------

### cumulative hazards----
# difference: cumhaz_0 - cumhaz_1
# exponential
cumhaz_exp <- function(dist, lam_b, cf_t, x_t, t, shape = NULL){
  lam <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  cumhaz_0 <- lam * (t - 1)
  cumhaz_1 <- lam * t
  return(cumhaz_0 - cumhaz_1)
}

# weibull
cumhaz_wei <- function(dist, lam_b, cf_t, x_t, t, shape = NULL){
  lam <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  cumhaz_0 <- lam * ((t - 1)^shape)
  cumhaz_1 <- lam * (t^shape)
  return(cumhaz_0 - cumhaz_1)
}

# gompertz
cumhaz_gom <- function(dist, lam_b, cf_t, x_t, t, shape = NULL){
  lam <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  cumhaz_0 <- (lam / shape) * (exp(shape * (t - 1)) - 1)
  cumhaz_1 <- (lam / shape) * (exp(shape * t) - 1)
  return(cumhaz_0 - cumhaz_1)
}

# cox
cumhaz_cox <- function(dist, lam_b, cf_t, x_t, t, shape = NULL){
  # shape is baseline hazard here
  lam <- lam_b * exp(x_t[names(cf_t)] %*% cf_t)
  cumhaz_0 <- log(lam) * shape[t - 1 + 1]
  cumhaz_1 <- log(lam) * shape[t + 1]
  return(cumhaz_0 - cumhaz_1)
}

cumhaz_list <- list(cumhaz_exp = cumhaz_exp, cumhaz_wei = cumhaz_wei,
                    cumhaz_gom = cumhaz_gom, cumhaz_cox = cumhaz_cox)

### probability of an event----
p_event <- function(dist, lam_b, cf_t, x_t, j, shape = NULL) {
  cumhaz_fn <- str_c("cumhaz_", dist)
  cumhaz_dif <- cumhaz_list[[cumhaz_fn]](lam_b = lam_b, cf_t = cf_t, x_t = x_t, t = j, shape = shape)
  return(1 - exp(cumhaz_dif))
}

### simulation by cycles----
gen_next_cycle <- function(prim_flag, j, v_j, 
                           v_b_int, int_names_b,int_names_noage,
                           nonlinage, v_b_int_nonlinage, int_names_b_nonlinage,
                           p_crv, 
                           events_nf_to_predict, events_nf_to_update, events_f, 
                           vars_t_eb, vars_t_e, eventsb,
                           t_since_eb_i, t_since_e_i,
                           alive,
                           dist, lam_b_i, cf_t_i, shape,
                           side_effects_i, had_mvevd, 
                           cvd_index_i, tx_cost, cost_int_str = NULL,
                           adj_cancer=FALSE, cancer_adj_rate_i, 
                           tx_effect_i, 
                           delay_age=0, stop_age=0, younger_45_i=FALSE,
                           adh_vec 
) { 
  
  
  ### time-updated characteristics ------------
  
  CurrAge <- v_j["CurrAge_cent"]
  v_j["cycle"] <- j
  
  # interactions
  # in case there is no-age time-varying variable in interaction
  v_multiply <- c(v_j, v_b_int)
  
  for (v_int in int_names_b){
    
    v_j[str_c("CurrAge_cent_int_", v_int)] <- CurrAge * v_multiply[v_int]
    
  }
  
  if (nonlinage) {
    # TODO: record flexpoint as parameter
    if (CurrAge < 1) {
      v_j["CurrAge_cent_1"] <- CurrAge
      v_j["CurrAge_cent_2"] <- 0
    } else {
      v_j["CurrAge_cent_1"] <- 1
      v_j["CurrAge_cent_2"] <- CurrAge - 1
    }
    for (v_int in int_names_b_nonlinage){
      if (grepl("dm_", v_int)) v_b_int_nonlinage[v_int] <- v_j[v_int] 
      for (gp in 1:2)
        v_j[str_c("CurrAge_cent_", gp, "_int_", v_int)] <- 
          v_j[str_c("CurrAge_cent_", gp)] * v_b_int_nonlinage[v_int]
    }
  }
  
  # event indicators for those events that have happened within 0-3 years
  
  # update time-updated events
  for (e in events_nf_to_update){
    t <- t_since_e_i[[e]] + 1
    # check whether increase changes time bracket
    vars_t_e_temp <- vars_t_e[[e]]
    change_bracket <- (t %in% vars_t_e_temp)
    if (change_bracket) {
      t_right_ind <- min(which(vars_t_e_temp >= t))
      # varname of the previous time bracket
      varname_prev <- if (t_right_ind == 1)
        str_c(e, "_", 0, "_", vars_t_e_temp[t_right_ind]) else 
          str_c(e, "_", vars_t_e_temp[t_right_ind - 1], "_", vars_t_e_temp[t_right_ind])
      # varname of the new time bracket
      # use tolower so that Inf is transformed to inf
      varname_new <- tolower(str_c(e, "_", vars_t_e_temp[t_right_ind], "_", vars_t_e_temp[t_right_ind + 1]))
      # update v_j
      v_j[varname_prev] <- 0
      v_j[varname_new] <- 1
    }
    # update time since event
    t_since_e_i[[e]] <- t
  }
  
  if (prim_flag == "prim" & v_j["two_plus_qol"] == 0){ 
    
    v_j["mi_1_inf_qol_pp"][v_j["mi_1_2"] == 1] <- 1 
    
    v_j["stroke_1_inf_qol_pp"][v_j["stroke_1_2"] == 1] <- 1
  }

  ### predict mi, stroke, crv & cancer -----------
  
  events_nf_rand <- sample(events_nf_to_predict)
  
  # if want CRV straight after MI
  if (!is.null(p_crv) & ("crv" %in% events_nf_rand) & ("mi" %in% events_nf_rand)) {
    # re-arrange if randomly draw a number >p_crv
    if (runif(1) > p_crv) {
      without_crv <- events_nf_rand[events_nf_rand != "crv"]
      without_crv_length <- length(without_crv)
      mi_n <- which(without_crv == "mi")
      events_nf_rand <- c(without_crv[1:mi_n], "crv")
      if (mi_n < without_crv_length)
        events_nf_rand <- c(events_nf_rand, without_crv[(mi_n + 1):without_crv_length])
    }
  }
  
  # to control delay of stop tx
  tx_lam_b <- lam_b_i
  drug_cost <- tx_cost$drug_cost
  
  
  delay <- if (drug_cost != 0 & delay_age != 0 & CurrAge < (delay_age-60)/10) TRUE else FALSE  
  
  stop <- if (drug_cost != 0 & stop_age !=0 &  
              (CurrAge > (stop_age-60)/10 | near(CurrAge, (stop_age-60)/10))) TRUE else FALSE
  
  delay_45_5y <- if (drug_cost != 0 & younger_45_i & j < 6) TRUE else FALSE
  
  # changing adherence 
  adh <- adh_vec[j]
  
  if (delay | stop | delay_45_5y | adh==0) {
    
    for (e in c("mi", "stroke", "crv", "vd")) {
      tx_lam_b[e] <- lam_b_i[e] / tx_effect_i[[e]]
    }
  }
  
  # to control increasing tx 
  if (incr_tx & j >= 6) {
    
    for (e in c("mi", "stroke", "crv", "vd")) {
      tx_lam_b[e] <- tx_lam_b[e] * (exp(v_multiply["ldl_red"])^(-0.0152*(j - 5)))
    }
  }
  
  for (e in events_nf_rand) {
    
    p_rand_e <- runif(1)
    # p_rand_e <- 0
    p_e <- p_event(dist = dist[[e]], lam_b = tx_lam_b[e], 
                   cf_t = cf_t_i[[e]], x_t = v_j, j = j, 
                   shape = shape[[e]])
  
    # treatment side effect on incident diabetes 
    if (e == "dm") {
      p_e <- p_e * side_effects_i$or_dm_low / (1 - p_e + p_e * side_effects_i$or_dm_low)
      p_e <- p_e * side_effects_i$or_dm_high / (1 - p_e + p_e* side_effects_i$or_dm_high)
    }
    
    # incident cancer rate adjustment
    if (adj_cancer == TRUE & e == "cancer_icd" & CurrAge >= 2) {# start at 80
      # use near instead of == due to difference storage format
      p_e <- p_e * cancer_adj_rate_i[near(cancer_adj_rate_i[, "CurrAge_cent"], CurrAge),
                                     "cancer_adj"]
    }
    
    if (p_e > p_rand_e) {
      events_nf_to_predict <- setdiff(events_nf_to_predict, e)
      # update time since event
      t_since_e_i[[e]] <- 0
      # update 0_t variable
      varname_0_1 <- str_c(e, "_0_", vars_t_e[[e]][1])
      # read off name of the 0_t variable
      v_j[varname_0_1] <- 1
      # if incident diabetes happens, recode no diabetes (dmPre) 0
      if (e == "dm") v_j[grepl("dmPre", names(v_j))] <- 0
    }
  }
  
  ### update interactions with events --------
  v_multiply <- c(v_j, v_b_int)
  
  for (v_int in grep("dm_|mi_|stroke_|crv_|cancer", int_names_b, value = T)) {
    
    v_j[str_c("CurrAge_cent_int_", v_int)] <- CurrAge * v_multiply[v_int]
    
  }
  
  # update noage interaction here, because them only exist for VD
  if (length(int_names_noage)>0) {
    
    for (i in 1:length(int_names_noage)) {
      v_j[[str_c(int_names_noage[[i]][1],"_noageint_", int_names_noage[[i]][2])]] <- 
        v_multiply[[int_names_noage[[i]][1]]]*v_multiply[[int_names_noage[[i]][2]]]
    } 
  }
  
  
  ### predict fatal events --------------
  
  events_f_rand <- sample(events_f)
  
  #e <- "nvd"
  
  for (e in events_f_rand) {
    
    p_rand_e <- runif(1)
    p_e <- p_event(dist = dist[[e]], lam_b = tx_lam_b[e], 
                   cf_t = cf_t_i[[e]], x_t = v_j, j = j, 
                   shape = shape[[e]])
    
    # add side effect - NVD death conditioning on rhabdo happen
    if (e == "nvd") p_e <- p_e + side_effects_i$p_nvd_rhabdo * side_effects_i$p_rhabdo
    
    if (p_e > p_rand_e) {
      alive <- 0
      v_j[e] <- 1
      # update d_without_e indicator for events that have not happened
      for (e_pre_death in events_nf_to_predict)
        v_j[str_c("d_without_", e_pre_death)] <- 1
      break
    }
  }
  
  ### hospital cost -----------------------------
  
  # hospital + primary care 
  
  # prepare interaction terms
  for (i in 1:nrow(cost_int_str)) {
    v_j[cost_int_str[i, 1]] <- v_j[cost_int_str[i, 2]] * v_j[cost_int_str[i, 3]]
  }
  
  cf_t_cost_hosp_prob <- cf_t_i[["cost_hosp_prob"]]
  cf_t_cost_hosp_cost <- cf_t_i[["cost_hosp_cost"]]
  cf_t_cost_prim <- cf_t_i[["cost_prim"]]
  
  # hospital
  # part one
  cost_hosp_prob_xb <- lam_b_i["cost_hosp_prob"] + 
    v_j[names(cf_t_cost_hosp_prob)] %*% cf_t_cost_hosp_prob
  
  cost_hosp_prob <- 1/(1 + exp(-cost_hosp_prob_xb))
  
  if (v_j["crv_0_1"] == 1) cost_hosp_prob <- 1 # crv in the same year incur cost 100%
  
  # part two
  cost_hosp_cost <- lam_b_i["cost_hosp_cost"] + 
    v_j[names(cf_t_cost_hosp_cost)] %*% cf_t_cost_hosp_cost
  
  # separate hospital and primary, tx and monitor
  v_j["hosp_cost"] <- cost_hosp_prob * cost_hosp_cost
  
  # adjust adverse effect on cost
  v_j["hosp_cost"] <- v_j["hosp_cost"] + 
    side_effects_i$cost_myop * side_effects_i$p_myop + side_effects_i$cost_rhabdo * side_effects_i$p_rhabdo
  
  # primary care
  cost_prim <- lam_b_i["cost_prim"] + 
    v_j[names(cf_t_cost_prim)] %*% cf_t_cost_prim
  
  # separate hospital and primary, tx and monitor
  v_j["prim_cost"] <- cost_prim
  

    ### sensitivity analyese CVD/incident diabetes costs only---------
  # conversely code all cvd and incident diabetes related costs to be 0
  # use the same name for simplicity
  
  cf_t_cost_hosp_prob_cvd <- cf_t_i[["cost_hosp_prob"]]
  cf_t_cost_hosp_cost_cvd <- cf_t_i[["cost_hosp_cost"]]
  cf_t_cost_prim_cvd <- cf_t_i[["cost_prim"]]

  cf_t_cost_hosp_prob_cvd[!grepl("cancer|nvd", names(cf_t_cost_hosp_prob_cvd))] <- 0
  cf_t_cost_hosp_cost_cvd[!grepl("cancer|nvd", names(cf_t_cost_hosp_cost_cvd))] <- 0
  cf_t_cost_prim_cvd[!grepl("cancer|nvd", names(cf_t_cost_prim_cvd))] <- 0

  if (is.na(t_since_eb_i[["dm"]])) {
    # means do not have baseline diabetes
    # do not count incident diabetes cost
    cf_t_cost_hosp_prob_cvd[grepl("dm", names(cf_t_cost_hosp_prob_cvd))] <- 0
    cf_t_cost_hosp_cost_cvd[grepl("dm", names(cf_t_cost_hosp_cost_cvd))] <- 0
    cf_t_cost_prim_cvd[grepl("dm", names(cf_t_cost_prim_cvd))] <- 0
  }

  # hospital
  # part one
  cost_hosp_prob_xb_cvd <- lam_b_i["cost_hosp_prob"] + 
    v_j[names(cf_t_cost_hosp_prob_cvd)] %*% cf_t_cost_hosp_prob_cvd

  cost_hosp_prob_cvd <- 1/(1 + exp(-cost_hosp_prob_xb_cvd))

  if (v_j["crv_0_1"] == 1) cost_hosp_prob_cvd <- 1 # crv in the same year incur cost 100%

  # part two
  cost_hosp_cost_cvd <- lam_b_i["cost_hosp_cost"] + 
    v_j[names(cf_t_cost_hosp_cost_cvd)] %*% cf_t_cost_hosp_cost_cvd

  # separate hospital and primary, tx and monitor
  hosp_cost_cvd <- cost_hosp_prob_cvd * cost_hosp_cost_cvd

  hosp_cost_cvd <- hosp_cost_cvd +
    side_effects_i$cost_myop * side_effects_i$p_myop + side_effects_i$cost_rhabdo * side_effects_i$p_rhabdo

  # primary care
  cost_prim_cvd <- lam_b_i["cost_prim"] + 
    v_j[names(cf_t_cost_prim_cvd)] %*% cf_t_cost_prim_cvd

  
  ### quality of life prediction -----
  cf_t_qol <- cf_t_i[["qol"]]
  
  if (v_j["two_plus_qol"] == 0 &
      ((v_j["mi_1_inf_qol_pp"]+v_j["stroke_1_inf_qol_pp"] == 2 & cvd_index_i == 0) | # pp
       (v_j["mi_1_2"] + v_j["stroke_1_2"] == 1 & cvd_index_i %in% c(1,4)) | # baseline othCHD only
       (v_j["stroke_1_2"] == 1 & cvd_index_i == 2) | # baseline mi only
       (v_j["mi_1_2"] == 1 & cvd_index_i == 3)# baseline stroke only 
      )){
    
    v_j["two_plus_qol"] <- 1
    v_j["mi_1_inf_qol_pp"] <- v_j["stroke_1_inf_qol_pp"] <- 0
  } 
  
  qol <- lam_b_i["qol"] + v_j[names(cf_t_qol)] %*% cf_t_qol
  
  # add side effect - directly using percentage on qol
  # qol <- qol - (0.017*30/365.25) * side_effects_i$p_myop
  # qol <- qol - (0.5*qol*7.5/365.25 + 0.2*qol*30/365.25) * side_effects_i$p_rhabdo
  
  qol <- qol - (side_effects_i$qol_myop * 30/365.25) * side_effects_i$p_myop
  qol <- qol - (side_effects_i$qol_rhabdo1 *qol*7.5/365.25 + 
                  side_effects_i$qol_rhabdo2 *qol*30/365.25) * side_effects_i$p_rhabdo
  

  ### sensitivity analyses on QoL -------
  # cancer-related QoL decrement back to HSE original coef: -0.127913479
  cf_t_qol_1 <- cf_t_qol
  cf_t_qol_1[grepl("cancer", names(cf_t_qol_1))] <- -0.127913479

  qol_can013 <- lam_b_i["qol"] + v_j[names(cf_t_qol_1)] %*% cf_t_qol_1

  qol_can013 <- qol_can013 - (side_effects_i$qol_myop * 30/365.25) * side_effects_i$p_myop
  qol_can013 <- qol_can013 - (side_effects_i$qol_rhabdo1 *qol_can013*7.5/365.25 +
                  side_effects_i$qol_rhabdo2*qol_can013*30/365.25) * side_effects_i$p_rhabdo

  # two diabetes-related QoL decrement halved
  cf_t_qol_2 <- cf_t_qol
  cf_t_qol_2[grepl("dm_", names(cf_t_qol_2))] <-
    cf_t_qol_2[grepl("dm_", names(cf_t_qol_2))]/2

  qol_dm50 <- lam_b_i["qol"] + v_j[names(cf_t_qol_2)] %*% cf_t_qol_2

  qol_dm50 <- qol_dm50 - (side_effects_i$qol_myop * 30/365.25) * side_effects_i$p_myop
  qol_dm50 <- qol_dm50 - (side_effects_i$qol_rhabdo1 *qol_dm50*7.5/365.25 +
                  side_effects_i$qol_rhabdo2*qol_dm50*30/365.25) * side_effects_i$p_rhabdo

  # Incident CV-events QoL decrement 50%
  cf_t_qol_3 <- cf_t_qol
  cf_t_qol_3[grepl("mi_|stroke_|crv_|two_", names(cf_t_qol_3))] <-
    cf_t_qol_3[grepl("mi_|stroke_|crv_|two_", names(cf_t_qol_3))]/2

  qol_cvdInc50 <- lam_b_i["qol"] + v_j[names(cf_t_qol_3)] %*% cf_t_qol_3

  qol_cvdInc50 <- qol_cvdInc50 - (side_effects_i$qol_myop * 30/365.25) * side_effects_i$p_myop
  qol_cvdInc50 <- qol_cvdInc50 - (side_effects_i$qol_rhabdo1 *qol_cvdInc50*7.5/365.25 +
                  side_effects_i$qol_rhabdo2*qol_cvdInc50*30/365.25) * side_effects_i$p_rhabdo

  # Incident CV-events QoL decrement 150%
  cf_t_qol_4 <- cf_t_qol
  cf_t_qol_4[grepl("mi_|stroke_|crv_|two_", names(cf_t_qol_4))] <-
    cf_t_qol_4[grepl("mi_|stroke_|crv_|two_", names(cf_t_qol_4))] * 1.5

  qol_cvdInc150 <- lam_b_i["qol"] + v_j[names(cf_t_qol_4)] %*% cf_t_qol_4
  
  qol_cvdInc150 <- qol_cvdInc150 - (side_effects_i$qol_myop * 30/365.25) * side_effects_i$p_myop
  qol_cvdInc150 <- qol_cvdInc150 - (side_effects_i$qol_rhabdo1 *qol_cvdInc150*7.5/365.25 +
                  side_effects_i$qol_rhabdo2*qol_cvdInc150*30/365.25) * side_effects_i$p_rhabdo


  
  # life year 
  ly <- 1
  
  ### treatment cost
  # to control delay of stop tx
  v_j["tx_cost"] <- if (delay | stop | delay_45_5y) 0 else drug_cost
  
  # add statin prescribing and monitoring cost into tx_cost
  # full cost regardless of living or death in the year
  # to control delay of stop tx
  
  if (drug_cost != 0) {
    
    if ((delay_age == 0 & j == 1 & !younger_45_i) | 
        (CurrAge > (delay_age-60)/10 & j == 1 & !younger_45_i) | 
        near(CurrAge, (delay_age-60)/10) & !younger_45_i | 
        (younger_45_i & j == 6)) {
      
      # initiation cost for every one: £42.6 #2019/20 price
      # cost in year 1 for all: 12.05 
      # separate hospital and primary, tx and monitor
      # cost inflated to 20/21, 2022-06-06
      # v_j["tx_cost_moni"] <- 56.33 #20/21 price
      v_j["tx_cost_moni"] <- tx_cost$monitor_cost$init_cost

    } else if ((prim_flag == "sec" | had_mvevd ==1) & !stop) {
      
      # separate hospital and primary, tx and monitor
      # cost inflated to 20/21, 2022-06-06
      # v_j["tx_cost_moni"] <- 12.42  #20/21 price
      v_j["tx_cost_moni"] <- tx_cost$monitor_cost$follow_cost
      
      
    } else v_j["tx_cost_moni"] <- 0
    
  } else v_j["tx_cost_moni"] <- 0
  
  ### death situation ------
  if (alive==0) {
    qol <- qol*0.5
    ly <- 0.5
    v_j["tx_cost"] <- 0.5 * v_j["tx_cost"]
    v_j["tx_cost_moni"] <- 0.5 * v_j["tx_cost_moni"]
    
    # sensitivity analyses
    hosp_cost_cvd <- hosp_cost_cvd*0.5
    cost_prim_cvd <- cost_prim_cvd*0.5
    qol_can013 <- qol_can013*0.5
    qol_dm50 <- qol_dm50*0.5
    qol_cvdInc50 <- qol_cvdInc50*0.5
    qol_cvdInc150 <- qol_cvdInc150*0.5
  }
  
  v_j["qol"] <- qol
  
  v_j["ly"] <- ly
  
  # leave the six slots in v_j for discounted values for the six sensitivity analysis values
  v_j["hosp_cost_cvd"] <- hosp_cost_cvd
  v_j["prim_cost_cvd"] <- cost_prim_cvd
  v_j["qol_can013"] <- qol_can013
  v_j["qol_dm50"] <- qol_dm50
  v_j["qol_cvdInc50"] <- qol_cvdInc50
  v_j["qol_cvdInc150"] <- qol_cvdInc150

  
  ### add mvevd -----
  v_j["mvevd"] <- 0
  
  if (had_mvevd == 0){
    
    if (v_j["mi_0_1"]==1 | v_j["stroke_0_1"] ==1 
        | v_j["crv_0_1"] ==1 | v_j["vd"] ==1){
      
      v_j["mvevd"] <- 1
      
      had_mvevd <- 1  
    }
  }
  
  # update time-varying baseline events: cancer, diabetes ---- 	
  
  # as baseline event have been initially defined, we should update it after using it in prediction 
  # update time-updated baseline variables
  # as repeated later with time-updated events
  
  dic <- c("dm"=0, "cancer_bsl"=1) 	
  # baseline diabetes start as 0_5, while baseline cancer start at 1_2 	
  for (eb in eventsb){
    # does the timer need to be updated?
    # ie, was there an event?
    if (!is.na(t_since_eb_i[[eb]])) {
      t <- t_since_eb_i[[eb]] + 1
      # check whether increase changes time bracket
      vars_t_eb_temp <- vars_t_eb[[eb]]
      change_bracket <- (t %in% vars_t_eb_temp)
      if (change_bracket) {
        t_right_ind <- min(which(vars_t_eb_temp >= t))
        # varname of the previous time bracket
        varname_prev <- if (t_right_ind == 1)
          str_c(eb, "_", dic[eb], "_", vars_t_eb_temp[t_right_ind]) else  
            str_c(eb, "_", vars_t_eb_temp[t_right_ind - 1], "_", vars_t_eb_temp[t_right_ind])
        # varname of the new time bracket
        # use tolower so that Inf is transformed to inf
        varname_new <- tolower(str_c(eb, "_", vars_t_eb_temp[t_right_ind], "_", vars_t_eb_temp[t_right_ind + 1]))
        # update v_j
        v_j[varname_prev] <- 0
        v_j[varname_new] <- 1
      }
      # update time since event
      t_since_eb_i[[eb]] <- t
    }
  }
  
  # update age at the end of cycle
  v_j["CurrAge_cent"] <- CurrAge + 1 / 10
  
  ### collate output & return
  ret_list <- list(v_j = v_j, 
                   events_nf_to_predict = events_nf_to_predict,  
                   alive = alive,
                   t_since_eb_i = t_since_eb_i, t_since_e_i = t_since_e_i,
                   had_mvevd = had_mvevd)
  
  
  return(ret_list)
}

###############################################################################
###############################################################################
###############################################################################

# Master function ---------------------------------------------------------

master <- function(rnum = rnum,
                   .input_data_dir, 
                   cf_filename,
                   prim_flag, 
                   mx_b_filename_prefix,
                   mx_t_filename_prefix,
                   events_list,
                   pf_filename,
                   adjust_crv, 
                   calibrated_eqns,
                   nonlinage,
                   tx,
                   dist,
                   stop_expr, 
                   sample_pat,
                   save_by_pat,
                   n_sim,
                   n_cores, 
                   output_dir, 
                   output_filename_prefix,
                   regimen = "none",
                   side_effects_flag = FALSE,
                   side_dm = FALSE,
                   side_mus = FALSE,
                   id_list = "all",
                   adj_cancer = FALSE, 
                   delay_age = 0, 
                   stop_age = 0,
                   delay_45 = FALSE,
                   adh_ctrl = FALSE,  
                   incr_tx = FALSE
                   
){
  
  # set the seed
  set.seed(1234)
  
  ### load data ------
  
  cf_hpc <- readRDS(file.path(.input_data_dir, 
                              str_c(cf_filename, ".rds")))
  
  cf_all <- cf_hpc[[paste0("cf_", rnum-1)]] # cf_0 is the first and deterministic
  
  mx_b <- readRDS(file.path(.input_data_dir, 
                            str_c(mx_b_filename_prefix, "_", prim_flag, ".rds")))
  mx_t <- readRDS(file.path(.input_data_dir, 
                            str_c(mx_t_filename_prefix, "_", prim_flag, ".rds")))
 
  pf <- readRDS(file.path(.input_data_dir,
                          str_c(pf_filename, ".rds")))
  
  p_crv <- if (adjust_crv == FALSE) NULL else 
    pf[["p_crv"]][[prim_flag]]
  
  # events 
  events_nf <- events_list$events_nf
  events_f <- events_list$events_f
  events_all <- c(events_nf, events_f)
  events_to_predict <- pf[["events_to_predict"]][[prim_flag]]
  
  
  # data on time-updated covariates
  vars_t <- pf[["vars_t"]]
  
  
  # duration since baseline cancer/diabetes
  t_since_eb <- pf[["t_since_eb"]][[prim_flag]]
  
  # incident cancer rate adjustment
  cancer_adj_rate <- pf[["cancer_adj_rate"]]
  
  adh_rate <- pf[["adh_rate"]] 
  
  # treatment effects and side effects
  tx_tag <- ""
  if (!identical(tx, FALSE)) {
    tx_tag <- tx
    tx_hpc <- readRDS(file.path(.input_data_dir, str_c(tx, ".rds")))
    # when tx is added, baseline LDL is replaced by pre-treated LDL 
    mx_b[,"NEWB_LDL_CL_cent"] <- mx_b[,"LDL_nostatin_cent"]
  }
  tx <- tx_hpc[[paste0("tx_", rnum-1)]]
  
  ### prepare coefficients --------------
  
  # identify interaction terms
  int_names_b <- sapply(strsplit(grep("_int_", colnames(mx_t), value = TRUE), "_int_"), "[", 2)
  int_names_b_nonlinage <- if (nonlinage)
    sapply(strsplit(grep("1_int_", colnames(mx_t), value = TRUE), "_int_"), "[", 2) else
      NULL
  
  # TODO
  int_names_noage <- strsplit(grep("_noageint_", colnames(mx_t), value = TRUE), "_noageint_")
  if (length(int_names_noage)==0) int_names_noage <- NULL
  
  # avoid error when no interaction terms
  if (length(int_names_b_nonlinage)==0) int_names_b_nonlinage <- NULL
  
  # extract correct equations
  if (calibrated_eqns) {
    # TODO: this seems to only extract calibrated equations for secondary?
    # need to check
    cf_temp <- list( 
      mi = cf_all[["mi"]][["calibrated"]][[prim_flag]][[dist[["mi"]]]],
      stroke = cf_all[["stroke"]][["calibrated"]][[prim_flag]][[dist[["stroke"]]]],
      crv = cf_all[["crv"]][["calibrated"]][[prim_flag]][[dist[["crv"]]]],
      cancer_icd = cf_all[["cancer"]][["calibrated"]][[prim_flag]][[dist[["cancer_icd"]]]], 
      vd = cf_all[["vd"]][["calibrated"]][[prim_flag]][[dist[["vd"]]]],
      nvd = cf_all[["nvd"]][["calibrated"]][[prim_flag]][[dist[["nvd"]]]], 
      dm = cf_all[["diabetes"]][[prim_flag]][[dist[["dm"]]]]
    ) 
  } else {
    cf_temp <- list(
      mi = cf_all[["mi"]][[prim_flag]][[dist[["mi"]]]],
      stroke = cf_all[["stroke"]][[prim_flag]][[dist[["stroke"]]]],
      crv = cf_all[["crv"]][[prim_flag]][[dist[["crv"]]]],
      cancer_icd = cf_all[["cancer"]][[prim_flag]][[dist[["cancer_icd"]]]], 
      vd = cf_all[["vd"]][[prim_flag]][[dist[["vd"]]]],
      nvd = cf_all[["nvd"]][[prim_flag]][[dist[["nvd"]]]],
      dm = cf_all[["diabetes"]][[prim_flag]][[dist[["dm"]]]]
    )
  }
  
  # cost equation 

  # cf_temp[["cost"]] <- cf_all[["cost"]][[prim_flag]]
  # hospital + primary care
  cf_temp[["cost_hosp"]] <- cf_all[["cost_hosp"]][[prim_flag]]
  cf_temp[["cost_prim"]] <- cf_all[["cost_prim"]][[prim_flag]]
  
  
  # QoL 
  # cf_temp[["qol"]] <- cf_all[["qol"]][[prim_flag]]
  cf_temp[["qol"]] <- cf_all[["qol"]]
  
  # baseline time-updated variables
  vars_t_eb <- vars_t[["b"]]
  eventsb <- names(vars_t_eb)
  # within-simulation updated variables
  vars_t_e <- vars_t[["e"]]
  
  # lam_b, cf_t and shape for all patients
  lam_b <- list()
  cf_t <- list()
  shape <- list()
  for (e in events_all) {
    cf_b <- cf_temp[[e]][["cf_b"]]
    
    if(e == "dm" | !calibrated_eqns) intercept_calibrated <- 0 else 
        intercept_calibrated <- cf_temp[[e]][["intercept_calibrated"]]
    
    lam_b[[e]] <- exp(mx_b[, names(cf_b)] %*% cf_b)* exp(intercept_calibrated)
    
    cf_t[[e]] <- cf_temp[[e]][["cf_t"]]
    shape[[e]] <- cf_temp[[e]][["shape"]]
  }

  # prepare tx effect ------
  
  # treat effect
  noeffect <- rep(1, nrow(mx_b))
  tx_effect <- list("mi"=noeffect, "stroke"=noeffect , "crv"=noeffect, "vd"=noeffect)
  
  if (!identical(tx, FALSE)) {
    
    # treatment regimen
    # get the absolute LDL reduction
    ldl_red <- (mx_b[, "NEWB_LDL_CL_cent"] + 3.6) * tx$statin_ldl_red[regimen]
    mx_b <- cbind(mx_b, ldl_red)
    
    for (e in c("mi", "stroke", "crv", "vd")){
      # generate a variable "ldl_red" of LDL reduction in mmol/L
      tx_effect[[e]] <- exp(mx_b[, "ldl_red"] * tx$tx_effect[[e]])
      lam_b[[e]] <- lam_b[[e]] * tx_effect[[e]]
    }
    
    # high-intensity statin or not
    statin_high <- tx$statin_ldl_red[regimen]>=0.45
    
    # treatment cost
    # drug cost
    tx_cost <- list()
    tx_cost$drug_cost <- tx$drug_annual_cost[regimen]
    tx_cost$monitor_cost <- tx$monitor_cost
    
  } else {
    ldl_red <- 0
    mx_b <- cbind(mx_b, ldl_red)
  }
  
  # treatment side effects ----
  
  side_effect_default <- list(or_dm_low = 1,
                              or_dm_high = 1,
                              p_myop = 0,
                              p_rhabdo = 0,
                              p_nvd_rhabdo = 0, 
                              cost_myop = 0, 
                              cost_rhabdo = 0,
                              qol_myop = 0,
                              qol_rhabdo1 = 0, 
                              qol_rhabdo2 = 0
                              )
  
  side_effects <- side_effect_default
  
  if (!identical(tx, FALSE)) {
    
    if (regimen=="none") side_effects_flag <- FALSE
    
    if (side_effects_flag) side_effects <- tx$side_effects
    
    if (side_effects_flag & !statin_high) side_effects$or_dm_high <- 1
    
    if (!side_dm) {
      side_effects$or_dm_low <- 1
      side_effects$or_dm_high <- 1
    }
    
    if (!side_mus) {
      side_effects$p_myop <- 0
      side_effects$p_rhabdo <- 0
      side_effects$p_nvd_rhabdo <- 0
    }
  }
  
  ### Hospital + primary care costs-----
  # hosptial  
  mx_b1 <- mx_b
  mx_b1[, "lnbcreann"] <- (mx_b1[, "lnbcreann"] - 4.4)/0.2
  
  cf_b_hosp_prob <- cf_temp[["cost_hosp"]][["p2_prob"]][["cf_b"]]
  lam_b[["cost_hosp_prob"]] <- mx_b1[, names(cf_b_hosp_prob)] %*% cf_b_hosp_prob
  cf_t[["cost_hosp_prob"]] <- cf_temp[["cost_hosp"]][["p2_prob"]][["cf_t"]]
  
  cf_b_hosp_cost <- cf_temp[["cost_hosp"]][["p2_cost"]][["cf_b"]]
  lam_b[["cost_hosp_cost"]] <- mx_b1[, names(cf_b_hosp_cost)] %*% cf_b_hosp_cost
  cf_t[["cost_hosp_cost"]] <- cf_temp[["cost_hosp"]][["p2_cost"]][["cf_t"]]
  
  # primary care
  cf_b_prim <- cf_temp[["cost_prim"]][["cf_b"]]
  lam_b[["cost_prim"]] <- mx_b[, names(cf_b_prim)] %*% cf_b_prim
  cf_t[["cost_prim"]] <- cf_temp[["cost_prim"]][["cf_t"]]
  
  # interation names in cost equations, cf_t
  cost_int <- unique(names(c(cf_t[["cost_hosp_prob"]], cf_t[["cost_hosp_cost"]], 
                             cf_t[["cost_prim"]])))
  cost_int <- cost_int[grepl("_int_", cost_int)]
  
  cost_int_str <- matrix(nrow = length(cost_int), ncol = 3)
  
  for (int in 1:length(cost_int)) {
    
    name_str <- unlist(strsplit(grep("_int_", cost_int[int], value = TRUE), "_int_"))
    
    cost_int_str[int, ] <- c(cost_int[int], name_str) 
  }
  
  ### Quality of life-----
  # cf_t[["qol"]] <- cf_temp[["qol"]][["cf_t"]]
  cf_b <- cf_temp[["qol"]][["cf_b"]]
  lam_b[["qol"]] <- mx_b[, names(cf_b)] %*% cf_b
  if (prim_flag == "prim") 
    cf_t[["qol"]] <- c(cf_temp[["qol"]][["cf_t"]], cf_temp[["qol"]][["no_cvd"]])
  
  
  ### template output matrix-----
  # independent of each simulation
  out_names <- c(colnames(mx_t)) 
  N <- length(out_names)
  
  vars_to_summarise_e <- str_c(names(vars_t_e), "_0_", unlist(lapply(vars_t_e, "[", 1)))
  # death variables
  vars_to_summarise_d <- grep("vd|without", out_names, value = TRUE)
  
  # addon variables
  var_addon <- c("hosp_cost", "prim_cost", "tx_cost", "tx_cost_moni", "qol", "ly",
                 "hosp_cost_cvd", "prim_cost_cvd", "qol_can013", 
                 "qol_dm50", "qol_cvdInc50", "qol_cvdInc150", "mvevd")
  
  # combine
  vars_to_summarise <- c(vars_to_summarise_e, vars_to_summarise_d, var_addon)
  
  # based on N
  v_output_length <- length(var_addon) + 1 # for nsim
  
  v_output_null <- matrix(nrow = 0, ncol = N + v_output_length)


  # Start loops -------------------------------------------------------------

  ### loop across patients-----
  
  if (is.numeric(sample_pat)) ids <- sort(sample(1:nrow(mx_b), sample_pat)) else
    if (id_list !="all")
      ids <- readRDS(file.path(.input_data_dir, str_c(id_list, ".rds")))[,1] else
        ids <- 1:nrow(mx_b) 

  set.seed(2021)
  
  n_patient <- nrow(mx_b)
  
  mat_rand <- matrix(sample(1:(n_patient*n_sim), n_sim * n_patient, replace = T), nrow = n_patient)
  
  #   # initiate parallel
  cl <- makeSOCKcluster(n_cores)
  clusterExport(cl, c("cumhaz_list", "p_event", "gen_next_cycle", 
                      "adj_cancer", "delay_45", "incr_tx", 
                      "stop_expr", "adh_ctrl","delay_age", "stop_age",
                      "save_by_pat", "output_dir", "output_filename_prefix"))

  clusterEvalQ(cl, library(tidyverse))
  clusterEvalQ(cl, library(data.table))
  registerDoSNOW(cl)
  
  retval <- parLapply(cl = cl, ids, function(i){
 
    # initialise output dataframe
    df_output_i <- v_output_null
    
    ### events to model
    events_to_predict_i <- events_to_predict[i, -1]
    events_nf_i <- names(events_to_predict_i)[which(events_to_predict_i == 1)]
    
    events_nf_to_update <- c() 
    
    # initiate baseline values
    # independent for each simulation
    v_b <- c(mx_b[i, ])
    v_b_int <- v_b[intersect(colnames(mx_b), c(int_names_b, unlist(int_names_noage)))]
    v_b_int_nonlinage <- v_b[int_names_b_nonlinage]
    
    # increasing tx 
    # add ldl_red in v_b_int for use in the loop
    v_b_int <- c(v_b_int, v_b["ldl_red"])
    
    # incident cancer rate adjustment
    if (v_b["male"]==0) cancer_adj_rate_i <- cancer_adj_rate$female else
      if (v_b["male"]==1) cancer_adj_rate_i <- cancer_adj_rate$male
    
    # move to outside simulation
    tx_effect_i <- list()
    for (e in c("mi", "stroke", "crv", "vd")) {
      tx_effect_i[[e]] <- tx_effect[[e]][i]
    }

    # create cf_t_i to avoid changing cf_t permanently 
    cf_t_i <- cf_t
    
    cvd_index_i <- 0
    
    if (prim_flag == "sec"){
      
      cf_t_i[["qol"]] <- cf_temp[["qol"]][["cf_t"]]
      
      if (v_b["othchd_only"]==1) {
        cvd_index_i <- 1
        cf_t_i[["qol"]] <- c(cf_temp[["qol"]][["cf_t"]], cf_temp[["qol"]][["othchd_only"]])
      }
      if (v_b["h_mi_raw_only"]==1) {
        cvd_index_i <- 2
        cf_t_i[["qol"]] <- c(cf_temp[["qol"]][["cf_t"]], cf_temp[["qol"]][["mi_only"]])
      }
      if (v_b["cvd_only"]==1) {
        cvd_index_i <- 3
        cf_t_i[["qol"]] <- c(cf_temp[["qol"]][["cf_t"]], cf_temp[["qol"]][["stroke_only"]])
      }
      if (v_b["h_pad_raw_only"]==1) {
        cvd_index_i <- 4
        cf_t_i[["qol"]] <- c(cf_temp[["qol"]][["cf_t"]], cf_temp[["qol"]][["pad_only"]])
      }
    }  
    
    younger_45_i <- FALSE
    
    if (delay_45) {
      
      younger_45_i <- if (mx_t[i, "CurrAge_cent"] < -1.5) TRUE else FALSE
      
    } 
    
    ### Monte Carlo simulation -----------------
    
    df_output_i <- lapply(1:n_sim, function(n){
      
      set.seed(mat_rand[i, n])
      
      # baseline value of time-update characteristics
      v_j <- mx_t[i, ]
      
      # stopping value
      J <- max(eval(stop_expr) - 1, 1)
      
      # default output matrix 
      mx_output_null <- matrix(nrow = J, ncol = N + v_output_length)
      colnames(mx_output_null) <- c(out_names, var_addon, "nsim") 
      
      # define / reset baseline values
      events_nf_to_predict <- events_nf_i
      alive <- 1
      # To check
      had_mvevd <- 0
      
      # side effect 
      side_effects_i <- side_effects

      lam_b_i <- sapply(lam_b, "[[", i)
      
      # time since events
      # NA indicates no event
      # baseline events
      t_since_eb_i <- list()
      for (eb in eventsb)
        t_since_eb_i[[eb]] <- as.numeric(t_since_eb[[eb]][which(t_since_eb[[eb]][, "ids"] == i), "t"])
      # simulated events
      t_since_e_i <- list()
      for (e in events_nf_i) 
        t_since_e_i[[e]] <- NA
      
      # adherence control
      if (adh_ctrl) { 
        
        rand_vec <- runif(J)
        
        adh_vec <- ifelse (adh_rate[1:J] > rand_vec, 1, 0)
        
      } else adh_vec <- rep(1, J)
      
      ### loop across years --------------------
      
      # output for each simulation
      mx_output <- mx_output_null
      
      #j <- 1
      for (j in 1:J) {
        
        alpha <- gen_next_cycle(prim_flag = prim_flag, 
                                j = j, 
                                v_j = v_j, 
                                v_b_int = v_b_int, 
                                int_names_b = int_names_b,
                                int_names_noage = int_names_noage,
                                nonlinage = nonlinage,
                                v_b_int_nonlinage = v_b_int_nonlinage, 
                                int_names_b_nonlinage = int_names_b_nonlinage,
                                p_crv = p_crv,
                                events_nf_to_predict = events_nf_to_predict, 
                                events_nf_to_update = events_nf_to_update, 
                                events_f = events_f,
                                eventsb = eventsb,
                                vars_t_eb = vars_t_eb, vars_t_e = vars_t_e,
                                t_since_eb_i = t_since_eb_i, t_since_e_i = t_since_e_i,
                                alive = alive,
                                dist = dist, 
                                lam_b_i = lam_b_i, 
                                cf_t_i = cf_t_i, 
                                shape = shape,
                                side_effects_i = side_effects_i,
                                had_mvevd = had_mvevd, 
                                cvd_index_i = cvd_index_i, 
                                tx_cost = tx_cost,
                                cost_int_str = cost_int_str,
                                adj_cancer = adj_cancer,
                                cancer_adj_rate_i = cancer_adj_rate_i, 
                                tx_effect_i = tx_effect_i, 
                                delay_age = delay_age,
                                stop_age = stop_age,
                                younger_45_i = younger_45_i,
                                adh_vec = adh_vec 
        )
        
        # save recursive output
        v_j <- alpha$v_j
        events_nf_to_predict <- alpha$events_nf_to_predict
        events_nf_to_update <- setdiff(events_nf_i, events_nf_to_predict)
        alive <- alpha$alive
        
        had_mvevd <- alpha$had_mvevd 
        
        # values of the time-updated covariates
        t_since_eb_i <- alpha$t_since_eb_i
        t_since_e_i <- alpha$t_since_e_i
        
        # record in the output matrix
        # automatically match mx_output and v_j
        v_j <- c(v_j, "nsim"=n)
        
        mx_output[j, ] <- v_j[colnames(mx_output)]
        
        # check whether the patient still alive
        if (alive == 0) {
          break
        }
      }
      
      return(mx_output) 
    })
    
    df_output_i <- map_df(df_output_i, ~as.data.frame(.x))
    
    # cleanup & return
    if (save_by_pat)
      saveRDS(df_output_i, compress = T, 
              file = file.path(output_dir, 
                               str_c(output_filename_prefix, "_pat_", i, "_", prim_flag, ".rds")))
    
    df_output_i <- df_output_i %>% as.data.frame() %>% group_by(nsim) %>%
      mutate(dm_0_10 = pmax(dm_0_10 - lag(dm_0_10, default = 0), 0)) 
    # use pmax to avoid minus when jumping from dm_0_10 to dm_10_inf
    
    # summarise
    df_combined <- data.table(df_output_i)[!is.na(id_new),
                                           lapply(.SD, function(x) sum(x) / n_sim), 
                                           by = .(id_new, cycle),
                                           .SDcols = vars_to_summarise]
    
    return(df_combined)
    
  }) 
  
  stopCluster(cl)
  
  # combine and cleanup
  df_output <- do.call("rbind", retval)
  
  # rename event columns
  for (e in events_nf) {
    varname_0_1 <- str_c(e, "_0_", vars_t_e[[e]][1])
    colnames(df_output)[which(colnames(df_output) == varname_0_1)] <- e
  }
  
  df_output <- as.data.frame(df_output)
  df_output <- df_output[ , -grep("_without_", colnames(df_output))]
  
  # save the output
  # tx_tag <- if (!identical(tx, FALSE)) "tx" else ""
  
  regimen_tag <- str_replace_all(regimen, fixed(" "), "") # remove space
  
  side_effct_tag <- if (side_effects_flag) "sideon" else "sideoff"
  
  if (side_effects_flag & !side_dm) side_effct_tag <-"sideondmoff"
  
  if (side_effects_flag & !side_mus) side_effct_tag <-"sideonmusoff"
  
  if (!side_dm & !side_mus) side_effct_tag <- "sideoff"
  
  adj_cancer_tag <- if (adj_cancer) "AdjCan_" else ""
  
  delay_tag <- if (delay_age > 0) paste0("D", delay_age, "_") else "" 
  
  stop_tag <- if (stop_age > 0) paste0("S", stop_age, "_") else ""
  
  delay_45_tag <- if (delay_45) "Delay45_" else ""
  
  sample_tag <- if (length(id_list)> 0 ) paste0(id_list, "_") else ""
  
  adh_ctrl_tag <- if (adh_ctrl) "Dyn_adh_" else ""
  
  incr_tx_tag <- if (incr_tx) "Incr_tx_" else ""
  
  # length_tag <- if (is.expression(stop_expr)) "lifetime" else paste0(stop_expr,"yr")
  
  saveRDS(df_output, compress = T, 
          file = file.path(output_dir, 
                           str_c(output_filename_prefix,"_",
                                 prim_flag, "_",
                                 tx_tag, "_",
                                 regimen_tag, "_", 
                                 side_effct_tag, "_",
                                 adj_cancer_tag, 
                                 delay_tag, 
                                 stop_tag, 
                                 delay_45_tag, 
                                 adh_ctrl_tag,
                                 incr_tx_tag,
                                 sample_tag, 
                                 n_sim, "_",
                                 rnum,
                                 ".rds")))
  
  return(df_output)
}

