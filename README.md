# The source codes and datasets of Iuchi and Hamada, 2020

[Figures 2 and 3](/Figures2and3/)  
[Figure 4](/Figure4/)  
[Figures 5, 6 and 7](/Figures5,6and7.R/)  
[Figure 8](/Figure8.R)

# Usage
```
set.seed(8)
#functions and a parameter to calculate tau
trial <- 1000
calc.null.dist <- function(dummy, d){
  d <- sample(d, length(d))
  seq1 <- d[1:(length(d) / 2)]
  seq2 <- d[(length(d) / 2 + 1):length(d)]
  
  tau <- 0
  for(i in 1:length(seq1)) for(j in 1:(length(seq1))) if(i < j)
    tau <- tau + (sign(seq1[j] - seq1[i]) * sign(seq2[j] - seq2[i]))
  tau <- tau / (length(seq1) * (length(seq1) - 1) / 2)
  return(tau)
}
calc.tau <- function(d){
  control.seq <- as.numeric(d)[(length(d) / 2 + 1):length(d)]
  case.seq <- as.numeric(d)[1:(length(d) / 2)]
  
  tau <- 0
  for(i in 1:length(control.seq)) for(j in 1:(length(control.seq))) if(i < j)
    tau <- tau + (sign(control.seq[j] - control.seq[i]) * sign(case.seq[j] - case.seq[i]))
  tau <- tau / (length(control.seq) * (length(control.seq) - 1) / 2)

  dist <- sapply(1:trial, calc.null.dist, d)
  return(c(tau = tau, p = sum(dist < tau) / trial))
}

#load the data
load("data.RData")

#calculate tau and p-value
result <- apply(data, 1, calc.tau)
```
