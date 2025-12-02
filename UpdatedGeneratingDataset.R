generate_vec <- function(shape1, rate1, shape2, rate2, tmax = 500) {
  y1 <- numeric()   # empty numeric vector
  y2 <- numeric()
  
  while (length(y1) == 0 || max(y1) < tmax) {
    t_i <- rgamma(1, shape1, scale = 1 / rate1)
    y1next <- ifelse(length(y1) == 0, t_i, y1[length(y1)] + t_i)
    if (y1next > tmax) {
      break
    }
    y1 <- c(y1, y1next)
  }
  
  while (length(y2) == 0 || max(y2) < tmax) {
    t_i <- rgamma(1, shape2, scale = 1 / rate2)
    y2next <- ifelse(length(y2) == 0, t_i, y2[length(y2)] + t_i)
    if (y2next > tmax) {
      break
    }
    y2 <- c(y2, y2next)
  }
  
  n <- length(y1)
  m <- length(y2)
  #Create Superposition vector and Category vector
  sup_vec <- rep(NA, n+m)
  cat_vec <- rep(NA, n+m)

  #Make counters to iterate through the two processes
  counter1 <- 1
  counter2 <- 1
  #Add a inf to the end of each process to get the whole superposition vec
  y1 <- append(y1, Inf)
  y2 <- append(y2, Inf)
  
  for (i in 1:(n+m)) {
    if(y1[counter1] <= y2[counter2]){
      sup_vec[i] <- y1[counter1]
      cat_vec[i] <- 1
      counter1 <- counter1 + 1
    }
    else{
      sup_vec[i] <- y2[counter2]
      cat_vec[i] <- 2
      counter2 <- counter2 + 1
    }
  }
  #Remove all obs before a certain random time
  random_time <- runif(1, 0, tmax/10)
  idx <- sup_vec > random_time
  sup_vec <- sup_vec[idx]
  cat_vec <- cat_vec[idx]
  sup_vec <- sup_vec - random_time
  tmax <- tmax - random_time
  #return the superposition vector, category vector and the new t_max
  return(list(sup_vec = sup_vec, cat_vec = cat_vec, tmax = tmax))
}