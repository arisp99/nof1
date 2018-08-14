#' Read json data in as an R object
#'
#' @param data input data. see sample input.json in the github repo
#' @param metadata metadata. see sample input.json in the github repo
#' @export

read_input_data <- function(data, metadata){

  if(metadata$user_age < 14){
    Outcome <- data$parent_response
  } else{
    if(sum(is.na(data$parent_response)) < sum(is.na(data$child_response))){
      Outcome <- data$parent_response
    }else{
      Outcome <- data$child_response
    }
  }

  Treatment <- data[,"treatment"]
  Treatment[Treatment == "strict"] = "A"
  Treatment[Treatment == "liberalized"] = "B"

  length_each <- sapply(Outcome[,"daily_stool_consistency"], length)
  Treat <- rep(Treatment, times = length_each)

  length_each2 <- sapply(Outcome[,"promis_gi_symptoms"]$`t-score`, length)
  Treatment_weekly <- rep(Treatment, time = length_each2)

  stool_consistency <- unlist(Outcome$daily_stool_consistency)
  stool_consistency <- ifelse(1 < stool_consistency & stool_consistency < 7, 0, 1)

  stool_frequency <- unlist(Outcome$daily_stool_frequency)
  pain_interference <- as.vector(unlist(Outcome$promis_pain_interference))
  gi_symptoms <- as.vector(unlist(Outcome$promis_gi_symptoms))

  list(Treatment = Treat, Treatment_weekly = Treatment_weekly, stool_consistency = stool_consistency, stool_frequency = stool_frequency, pain_interference = pain_interference,  gi_symptoms = gi_symptoms)
}

washout <- function(read_data){

  with(read_data,{

    change_point <- cumsum(rle(Treatment)$lengths)
    change_point <- change_point[-length(change_point)]

    delete_obs_daily <- NULL
    for(i in 1:length(change_point)){
      delete_obs_daily <- c(delete_obs_daily, (change_point[i]+1):(change_point[i]+7))
    }
    delete_obs_daily


    change_point2 <- cumsum(rle(Treatment_weekly)$lengths)
    change_point2 <- change_point2[-length(change_point2)]
    delete_obs_weekly <- NULL
    for(i in 1:length(change_point2)){
      delete_obs_weekly <- c(delete_obs_weekly, (change_point2[i]+1))
    }
    delete_obs_weekly

    stool_consistency[delete_obs_daily] <- NA
    stool_frequency[delete_obs_daily] <- NA
    pain_interference[delete_obs_weekly] <- NA
    gi_symptoms[delete_obs_weekly] <- NA

    list(Treatment = Treatment, Treatment_weekly = Treatment_weekly, stool_consistency = stool_consistency, stool_frequency = stool_frequency, pain_interference = pain_interference, gi_symptoms = gi_symptoms)
  })
}


change <- function(x){
  x = ifelse(x==0,1,x)
  x = ifelse(x==100,99,x)
  return(x)
}

check_enough_data <- function(Treatment, x){
  length(table(Treatment[!is.na(x)])) == 3
}

check_success <- function(x){
  ifelse(is.list(x), TRUE, x)
}

link_function <- function(x, response){
  answer <-
    if(response == "poisson"){
      exp(x)
    } else if(response == "binomial"){
      inv_logit(x)
    } else if(response == "normal"){
      x
    }
}

inv_logit <- function(a){
  1/(1+exp(-a))
}

find_raw_mean <- function(Y, Treat, response){

  raw_mean <- c(mean(Y[Treat == "baseline"], na.rm = TRUE), mean(Y[Treat == "A"], na.rm = TRUE), mean(Y[Treat == "B"], na.rm = TRUE))
  raw_mean[is.nan(raw_mean)] <- NA
  raw_mean
}

round_number <- function(raw_mean, response){

  if(response == "poisson" || response == "normal"){
    round(raw_mean,1)
  } else if(response == "binomial"){
    round(raw_mean*100)
  }
}

transform_using_link <- function(coef, response){

  coef_alpha <- coef_beta_A <- coef_beta_B <- NA

  if("alpha" %in% colnames(coef)){
    coef_alpha <- coef[,"alpha", drop = F]
  }

  if("beta_A" %in% colnames(coef)){
    coef_beta_A <- coef[,"beta_A", drop = F]
  }

  if("beta_B" %in% colnames(coef)){
    coef_beta_B <- coef[,"beta_B",drop = F]
  }

  base <- link_function(coef_alpha, response)
  scd <- link_function(coef_alpha + coef_beta_A, response)
  mscd <- link_function(coef_alpha + coef_beta_B, response)

  return(list(base = base, scd = scd, mscd = mscd))
}


find_mean_difference <- function(coef, response, raw_mean){

  trans <- transform_using_link(coef, response)

  mean_difference <- with(trans, {
    c(base_vs_scd = mean(scd - base), base_vs_mscd = mean(mscd - base), mscd_vs_scd = mean(scd - mscd))
  })

  rounded <- round_number(mean_difference, response)

  if(response == "binomial"){
    rounded[rounded==0 & !is.na(rounded)] <- 1
    rounded[rounded==100  & !is.na(rounded)] <- 99
  }
  return(rounded)
}

calculate_p_threshold <- function(coef, response){

  upper <-
    if(response == "poisson"){
      1.1
    } else if(response == "binomial"){
      1.1
    } else if(response == "normal"){
      2.9
    }

  lower <-
    if(response == "poisson"){
      0.9
    } else if(response == "binomial"){
      0.9
    } else if(response == "normal"){
      -2.9
    }

  trans <- transform_using_link(coef, response)

  with(trans, {

  if("beta_A" %in% colnames(coef)){
    if(response == "normal"){
      base_vs_scd <- list(greater_than_threshold = round(mean(scd - base > upper, na.rm = TRUE)*100),
                          lower_than_threshold = round(mean(scd - base < lower, na.rm = TRUE)*100))
    } else if(response %in% c("poisson", "binomial")){
      base_vs_scd <- list(greater_than_threshold = round(mean(scd/base > upper, na.rm = TRUE)*100),
                          lower_than_threshold = round(mean(scd/base < lower, na.rm = TRUE)*100))
    }
    base_vs_scd <- rapply(base_vs_scd, change, how = "replace")
  } else{
    base_vs_scd <- list(greater_than_threshold = NA, lower_than_threshold = NA)
  }

  if("beta_B" %in% colnames(coef)){
    if(response == "normal"){
      base_vs_mscd <- list(greater_than_threshold = round(mean(mscd - base > upper, na.rm = TRUE)*100),
                           lower_than_threshold = round(mean(mscd - base < lower, na.rm = TRUE)*100))
    } else if(response %in% c("poisson", "binomial")){
      base_vs_mscd <- list(greater_than_threshold = round(mean(mscd/base > upper, na.rm = TRUE)*100),
                          lower_than_threshold = round(mean(mscd/base < lower, na.rm = TRUE)*100))
    }
    base_vs_mscd <- rapply(base_vs_mscd, change, how = "replace")
  } else{
    base_vs_mscd <- list(greater_than_threshold = NA, lower_than_threshold = NA)
  }

  if("beta_A" %in% colnames(coef) & "beta_B" %in% colnames(coef)){

    if(response == "normal"){
      mscd_vs_scd <- list(greater_than_threshold = round(mean(scd - mscd > upper, na.rm = TRUE)*100),
                           lower_than_threshold = round(mean(scd - mscd < lower, na.rm = TRUE)*100))
    } else if(response %in% c("poisson", "binomial")){
      mscd_vs_scd <- list(greater_than_threshold = round(mean(scd/mscd > upper, na.rm = TRUE)*100),
                          lower_than_threshold = round(mean(scd/mscd < lower, na.rm = TRUE)*100))
    }
    mscd_vs_scd <- rapply(mscd_vs_scd, change, how = "replace")
  } else{
    mscd_vs_scd <- list(greater_than_threshold = NA, lower_than_threshold = NA)
  }

  return(list(base_vs_scd = base_vs_scd, base_vs_mscd = base_vs_mscd, mscd_vs_scd = mscd_vs_scd))
  })
}

#' Summarizes the result from the model into json format
#'
#' @export

summarize_nof1 <- function(nof1, result){

  with(c(nof1, result),{

    samples <- do.call(rbind, samples)
    raw_mean <- find_raw_mean(Y, Treat, response)
    rounded_raw_mean <- round_number(raw_mean, response)

    coef <- samples[,colnames(samples) %in% c("alpha", "beta_A", "beta_B")]
    diff <- find_mean_difference(coef, response, raw_mean)

    raw_mean <- list(base = rounded_raw_mean[1], scd = rounded_raw_mean[2], mscd = rounded_raw_mean[3])
    mean_difference <- list(base_vs_scd = diff[1], base_vs_mscd = diff[2], mscd_vs_scd = diff[3])
    three_group_comparison <- list(raw_mean = raw_mean, mean_difference = mean_difference)

    gauge_graph <- calculate_p_threshold(coef, response)
    return(list(three_group_comparison = three_group_comparison, gauge_graph = gauge_graph))
  })
}



#' Wrapper function that runs the n-of-1 model
#'
#' @param json.file input json data
#' @export

wrap <- function(data, metadata){

  read_data <- tryCatch({
    read_dummy <- read_input_data(data, metadata)
    if(length(rle(read_dummy$Treatment)$lengths) > 1) read_dummy <- washout(read_dummy)
    read_dummy
  }, error = function(error){
    return(paste("input read error: ", error))
  })

  print(read_data)

  stool_frequency <- tryCatch({
    data_freq <- list(Treat = read_data$Treatment, Y = read_data$stool_frequency)
    nof1_freq <- with(data_freq, {
      nof1.data(Y, Treat, response = "poisson")
    })
    result_freq <- nof1.run(nof1_freq)
    summarize_nof1(nof1_freq, result_freq)
  }, error = function(error){
    return(paste("stool_frequency run error: ", error))
  })

  stool_consistency <- tryCatch({
    data_cons <- list(Treat = read_data$Treatment, Y = read_data$stool_consistency)
    nof1_cons <- with(data_cons, {
      nof1.data(Y, Treat, response = "binomial")
    })
    result_cons <- nof1.run(nof1_cons)
    summarize_nof1(nof1_cons, result_cons)
  }, error = function(error){
    return(paste("stool_consistency run error: ", error))
  })

  pain_interference <- tryCatch({
    data_pain <- list(Treat = read_data$Treatment_weekly, Y = read_data$pain_interference)
    nof1_pain <- with(data_pain, {
      nof1.data(Y, Treat, response = "normal")
    })
    result_pain <- nof1.run(nof1_pain)
    summarize_nof1(nof1_pain, result_pain)
  }, error = function(error){
    return(paste("pain_interference run error: ", error))
  })

  gi_symptoms <- tryCatch({
    data_gi <- list(Treat = read_data$Treatment_weekly, Y = read_data$gi_symptoms)
    nof1_gi <- with(data_gi, {
      nof1.data(Y, Treat, response = "normal")
    })
    result_gi <- nof1.run(nof1_gi)
    summarize_nof1(nof1_gi, result_gi)
  }, error = function(error){
    return(paste("gi_symptoms run error: ", error))
  })

  metadata <- list(successful_input_reading = check_success(read_data),
                   successful_run_stool_frequency = check_success(stool_frequency),
                   successful_run_stool_consistency = check_success(stool_consistency),
                   successful_run_pain_interference = check_success(pain_interference),
                   successful_run_gi_symptoms = check_success(gi_symptoms),
                   enough_stool_consistency = check_enough_data(read_data$Treatment, read_data$stool_consistency),
                   enough_stool_frequency = check_enough_data(read_data$Treatment, read_data$stool_frequency),
                   enough_pain_interference = check_enough_data(read_data$Treatment_weekly, read_data$pain_interference),
                   enough_gi_symptoms = check_enough_data(read_data$Treatment_weekly, read_data$gi_symptoms),
                   user_id = metadata$user_id,
                   timestamp_trialist_completed = Sys.time(),
                   trialist_version_id = 2,
                   trialist_version_date = "12/04/2017",
                   trialist_version_note = "")

  summary_graph <- tryCatch({
    base_vs_scd <- find_summary_graph(metadata, "base_vs_scd", stool_frequency, stool_consistency, pain_interference, gi_symptoms)
    base_vs_mscd <- find_summary_graph(metadata, "base_vs_mscd", stool_frequency, stool_consistency, pain_interference, gi_symptoms)
    mscd_vs_scd <- find_summary_graph(metadata, "mscd_vs_scd", stool_frequency, stool_consistency, pain_interference, gi_symptoms)
    list(base_vs_scd = base_vs_scd, base_vs_mscd = base_vs_mscd, mscd_vs_scd = mscd_vs_scd)
  }, error = function(error){
    return(paste("error in summary step:", error))
  })

  metadata2 <- list(successful_summary_graph = check_success(summary_graph))
  metadata <- c(metadata2, metadata)

  final <- list(metadata = metadata, stool_frequency = stool_frequency, stool_consistency = stool_consistency,
                pain_interference = pain_interference, gi_symptoms = gi_symptoms, summary_graph = summary_graph)

  return(final)
}

find_summary_graph <- function(metadata, treatment_comparison,stool_frequency, stool_consistency, pain_interference, gi_symptoms){

  summary <- list()

  if(metadata$successful_run_stool_frequency == TRUE){
    a <- stool_frequency[["gauge_graph"]][[treatment_comparison]][["greater_than_threshold"]]
    b <- stool_frequency[["gauge_graph"]][[treatment_comparison]][["lower_than_threshold"]]
    summary$stool_frequency <- round(ifelse(a >= b, -a, b))
  }

  if(metadata$successful_run_stool_consistency == TRUE){
    c <- stool_consistency[["gauge_graph"]][[treatment_comparison]][["greater_than_threshold"]]
    d <- stool_consistency[["gauge_graph"]][[treatment_comparison]][["lower_than_threshold"]]
    summary$stool_consistency <- round(ifelse(c >= d, -c, d))
  }

  if(metadata$successful_run_pain_interference == TRUE){
    e <- pain_interference[["gauge_graph"]][[treatment_comparison]][["greater_than_threshold"]]
    f <- pain_interference[["gauge_graph"]][[treatment_comparison]][["lower_than_threshold"]]
    summary$pain_interference <- round(ifelse(e >= f, -e, f))
  }

  if(metadata$successful_run_gi_symptoms == TRUE){
    g <- gi_symptoms[["gauge_graph"]][[treatment_comparison]][["greater_than_threshold"]]
    h <- gi_symptoms[["gauge_graph"]][[treatment_comparison]][["lower_than_threshold"]]
    summary$gi_symptoms <- round(ifelse(g >= h, -g, h))
  }
  return(summary)
}

