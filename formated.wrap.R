#Input File Format:
#The file has two parts, the data section and the metadata section. The data section will contain all the observations taken during the trial. For each observation, we will have a list of the treatment the patient was on at the time followed by a list of all the data points collected. These two lists should be the same length. The metadata section will then contain the  user_id, the trigger, the design of the trial, whether or not a washout  period will be used (if not in metadata, default is TRUE), the alpha value for the confidence interval, the start and end  date (YYYY-MM-DD), followed by the response type of each of the observations. Either binomial, poisson, or normal. The response types must be in the same order as the observations in the data section and must be the last pieces of information in the metadata section.

#To read in a data set, as a .json file use the fromJSON command.
#ex: json.file <- fromJSON("NAME-OF-FILE.json")

#Some packages needed for file to run
library(nof1)
library(jsonlite)

#This function essentially reads in the data and changes the data that was entered binomially to contain only 0 and 1.
formated_read_input <- function(data, metadata){

  outcome_names <- names(data)
  output <- data
  #chaning binomials to 0, 1
  for(i in 1:length(outcome_names)){
  	if(metadata[[length(metadata)-length(outcome_names)+i]] == "binomial"){
  		fixed <- unlist(data[,i]$result)
  		if(all(fixed %in% c(0, 1, NA))){
  			output[[outcome_names[i]]]$result = list(fixed)
  		}
  		else{
  		   max <- max(fixed, na.rm=TRUE)
  		   fixed <- ifelse(max > fixed & (!is.nan(fixed)), 0, 1)
  		   output[[outcome_names[i]]]$result = list(fixed)
  		}
  	}
  }
  output
}

#Washout function. in some cases when we switch from A to B, for example, the first couple data points in B could be corrupted because the effects of A are still there. Thus we ingore the first couple data points (set them to NA).
washout <- function(read_data, metadata){

  #calculating the duration of the study
  finish <- as.Date(metadata$trial_end_date, format="%Y-%m-%d")
  start <- as.Date(metadata$trial_start_date, format="%Y-%m-%d")
  date_diff <-as.numeric(finish-start+1)

  #calculating the number of days and weeks
  num_days <- date_diff
  num_weeks <- date_diff/7

  #loop over the number of different observations we have
  for(i in 1:length(change_point)){
    vec_treat <- unlist(read_data[[i]]$treatment)
   	change_point <- cumsum(rle(vec_treat)$lengths)
   	change_point <- change_point[-length(change_point)]

    #we only change weekly results and daily results, more than weekly does
    #not make much sense as enough time would have passed so that the previous 	#treatment would not make a difference anymore
    #NB: at some point, may want to account for different frequencies, ie.
    #multiple time a day
    delete_obs_daily <- NULL
    delete_obs_weekly <- NULL

  	if(length(vec_treat) == num_days){
   		for(j in 1:length(change_point)){
      		delete_obs_daily <- c(delete_obs_daily,
      		(change_point[j]+1):(change_point[j]+7))
   		}
    	delete_obs_daily
    }
    else if(length(vec_treat) == num_weeks){
    	for(j in 1:length(change_point)){
     		delete_obs_weekly <- c(delete_obs_weekly,
     		(change_point[j]+1))
    	}
    	delete_obs_weekly
    }
    read_data[[i]]$result[[1]][delete_obs_daily] <- NA
    read_data[[i]]$result[[1]][delete_obs_weekly] <- NA
  }
  read_data
}

#Finds the raw mean of the input vector. Returns raw mean for baseline and all other treatments.
find_raw_mean <- function(Y, Treat){

  raw_mean <- c(mean(Y[Treat == "baseline"], na.rm = TRUE))
  for(i in 1:(length(unique(unlist(json.file$data[,1]$treat)))-1)){
    raw_mean <- c(raw_mean, mean(Y[Treat == LETTERS[i]], na.rm = TRUE))
    }
  raw_mean[is.nan(raw_mean)] <- NA
  raw_mean
}

#Finds the raw median of the input vector. Returns raw median for baseline and all other treatments.
find_raw_median <- function(Y, Treat){

  raw_median <- c(median(Y[Treat == "baseline"], na.rm = TRUE))
  for(i in 1:(length(unique(unlist(json.file$data[,1]$treat)))-1)){
    raw_median <- c(raw_median, median(Y[Treat == LETTERS[i]], na.rm = TRUE))
    }
  raw_median[is.nan(raw_median)] <- NA
  raw_median
}

#Checking that for each outcome given, we have the correct number of treatments
check_nof_treatments <- function(treatment, data, nof_treat){
  length(table(treatment[!is.na(data)])) == nof_treat
}

#i dont think this is needed at all. i dont rlly see the point of it tbh...
# check_success <- function(x){
  # ifelse(is.list(x), TRUE, x)
# }

#Rounds the raw mean. Unchanged if poisson or normal, multiply by 100 if binomial
round_number <- function(raw_mean, response){

  if(response == "poisson" || response == "normal"){
    round(raw_mean,1)
  } else if(response == "binomial"){
    round(raw_mean*100)
  }
}

#This was used for calculate p-threshold, which is used for the graphs...
#Other stuff used for graphs not implemented: find_mean_difference, calculate_p_threshold, find_summary_graph
change <- function(x){
  x = ifelse(x==0,1,x)
  x = ifelse(x==100,99,x)
  return(x)
}

#Our model file will not present our posterior distribution in an explicit form, we may need to exponentiate or use the inverse logit function.
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

#Used for link_function. Defines inv_logit.
inv_logit <- function(a){
  1/(1+exp(-a))
}

#Function to present a summary of our results. When we take our draws we have coefficients: alpha, beta_A, beta_B, etc. Baseline=alpha, treatmentA=alpha+beta_A, treatmentB=alpha+beta_B, etc. Function returns: input mean and median, mean and median for the treatments, mean and median for the coefs, P(treatment>0), P(coef>0), CI(alpha=0.5) for treatment and coef.
summarize_nof1 <- function(nof1, result, nof_treat, alpha){

  with(c(nof1, result),{

    samples <- do.call(rbind, samples)

    #creating our list of treatment names. base, treat_A, treat_B, etc.
    treat_n = list("base")
    for(i in 1:(nof_treat-1)){
    	treat_n <- c(treat_n, paste("treat", LETTERS[i], sep="_"))
    }

    #mean and median for input vectors
    input_mean = c(mean(Y[Treat == "baseline"], na.rm = TRUE))
    input_median = c(median(Y[Treat == "baseline"], na.rm = TRUE))
    for(i in 1:(nof_treat-1)){
    	input_mean <-
    	c(input_mean, mean(Y[Treat == LETTERS[i]], na.rm = TRUE))
    	input_median <-
    	c(input_median, median(Y[Treat == LETTERS[i]], na.rm = TRUE))
    	}

    input_mean[is.nan(input_mean)] <- NA
    input_median[is.nan(input_median)] <- NA

	#rounding mean and median
    input_mean <- sapply(input_mean, round_number, response)
    input_median <- sapply(input_median, round_number, response)

    names(input_mean) <- treat_n
    names(input_median) <- treat_n

    #creating treatment vectors
    treatment = list(samples[,1])
    for(i in 2:nof_treat){
    	treatment <-
    	c(treatment, list(link_function(samples[,1] + samples[,i], response)))
    }
    names(treatment) <- treat_n

    #mean and median for treatment lists
    treat_mean <- list()
    treat_median <- list()
    treat_mean <- sapply(treatment, mean)
    treat_median <- sapply(treatment, median)
    names(treat_mean) <- treat_n
    names(treat_median) <- treat_n

    #mean and median for coef draws
    coef_mean <- list()
    coef_median <- list()
    for(i in 1:nof_treat){
    	coef_mean <- c(coef_mean, mean(samples[,i]))
    	coef_median <- c(coef_median, median(samples[,i]))
    }
    names(coef_mean) <- treat_n
    names(coef_median) <- treat_n
    coef_mean <- unlist(coef_mean)
    coef_median <- unlist(coef_median)

    #probability that the treatment draw is greater than zero
    treat_greater_zero <- list()
    for(i in 1:nof_treat){
    	treat_greater_zero <-
    	c(treat_greater_zero, sum(treatment[[i]]>0)/length(treatment[[i]]))
    }
    names(treat_greater_zero) <- treat_n
    treat_greater_zero <- unlist(treat_greater_zero)

    #probability that the coef draw is greater than zero
    coef_greater_zero <- list()
    for(i in 1:nof_treat){
    	coef_greater_zero <-
    	c(coef_greater_zero, sum(samples[,i]>0)/length(samples[,i]))
    }
    names(coef_greater_zero) <- treat_n
    coef_greater_zero <- unlist(coef_greater_zero)

    #confidence interval for the treatment draw.
    treat_confidence <- list()
    for(i in 1:nof_treat){
    	treat_confidence <-
    	rbind(treat_confidence, quantile(treatment[[i]], c(alpha, 1-alpha)))
    }
    rownames(treat_confidence) <- treat_n

    #confidence interval for the coef draw.
    coef_confidence <- list()
    for(i in 1:nof_treat){
    	coef_confidence <-
    	rbind(coef_confidence, quantile(samples[,i], c(alpha, 1-alpha)))
    }
    rownames(coef_confidence) <- treat_n

    final <- list(
    input_mean = input_mean,
    input_median = input_median,
    treat_mean = treat_mean,
    treat_median = treat_median,
    coef_mean = coef_mean,
    coef_median = coef_median,
    treat_greater_zero = treat_greater_zero,
    coef_greater_zero = coef_greater_zero,
    treat_confidence = treat_confidence,
    coef_confidence = coef_confidence
    )

    return(final)
  })
}

#Wrapper function to run the nof1 model
gen_wrap <- function(data, metadata){

  #getting our data read and formated
  read_data <- tryCatch({
    read_dummy <- formated_read_input(data, metadata)
    if(metadata$washout == "FALSE" || is.null(metadata$washout)){
    	read_dummy
    }
    else{
    	read_dummy <- washout(read_dummy, metadata)
    	read_dummy
    }
  }, error = function(error){
    return(paste("input read error: ", error))
  })

  print(read_data)

  #initializing some values
  names <- names(data)
  nof_responses <- length(data)
  read_len <- length(read_data)
  nof_treat <- length(unique(unlist(data[,1]$treat)))

  #running our algorithm
  returns = list()
  for(i in 1:nof_responses){
    returns[[i]] <- wrap_helper(read_data[,i], names[i], metadata[length(metadata)-nof_responses+i], nof_treat, as.numeric(metadata["confidence"])/2)
  }
  names(returns) <- as.list(names)

  #checking that the number of treatments is correct for each observation
  correct_treat <- list()
  for(i in 1:nof_responses){
  	correct_treat <-
  	c(correct_treat,
  	check_nof_treatments(
  	read_data[,i]$treatment,
  	read_data[,i]$result,
  	nof_treat))
  }
  names(correct_treat) <- as.list(names)

  #listing some info about algorithm run
  system_info <- list(enough_data = correct_treat,
                      user_id = metadata$user_id,
                      trigger = metadata$trigger,
                      design = metadata$design,
                      timestap_completion = Sys.time(),
                      version_id = 1,
                      version_date = "8/10/2018"
                      )

  final <- list(system_info = system_info, model_results = returns)

  return(final)
}

#Helper function for our wrap function. Creates the neccesary objects needed to run the nof1 algorithim, runs the algorithim, and then summarizes the result.
wrap_helper <-
function(specific_data, outcome_name, response_type, nof_treat, alpha){
  summary <- tryCatch({
  	data_out <-
  	list(Treat = specific_data$treatment[[1]], Y = specific_data$result[[1]])
    nof1_out <- with(data_out, {
      nof1.data(Y, Treat, response = response_type)
    })
    result_out <- nof1.run(nof1_out)
    summarize_nof1(nof1_out, result_out, nof_treat, alpha)
  }, error = function(error){
    return(paste(outcome_name, "run error: ", error))
  })
}
