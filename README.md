# The source codes and datasets of Iuchi and Hamada, 2020

[Figures 2 and 3](/Figures2and3/)  
[Figure 4](/Figure4/)  
[Figures 5, 6 and 7](/Figures5,6and7.R/)  
[Figure 8](/Figure8.R)

# Usage
`set.seed(8)<br>
#functions and a parameter to calculate tau<br>
trial <- 1000<br>
calc.null.dist <- function(dummy, d){<br>
  d <- sample(d, length(d))<br>
  seq1 <- d[1:(length(d) / 2)]<br>
  seq2 <- d[(length(d) / 2 + 1):length(d)]<br>
  
  tau <- 0<br>
  for(i in 1:length(seq1)) for(j in 1:(length(seq1))) if(i < j)<br>
    tau <- tau + (sign(seq1[j] - seq1[i]) * sign(seq2[j] - seq2[i]))<br>
  tau <- tau / (length(seq1) * (length(seq1) - 1) / 2)<br>
  return(tau)<br>
}<br>
calc.tau <- function(d){<br>
  control.seq <- as.numeric(d)[(length(d) / 2 + 1):length(d)]<br>
  case.seq <- as.numeric(d)[1:(length(d) / 2)]<br>
  
  tau <- 0<br>
  for(i in 1:length(control.seq)) for(j in 1:(length(control.seq))) if(i < j)<br>
    tau <- tau + (sign(control.seq[j] - control.seq[i]) * sign(case.seq[j] - case.seq[i]))<br>
  tau <- tau / (length(control.seq) * (length(control.seq) - 1) / 2)<br>

  dist <- sapply(1:trial, calc.null.dist, d)<br>
  return(c(tau = tau, p = sum(dist < tau) / trial))<br>
}<br>

#load the data<br>
load("data.RData")<br>

#calculate tau and p-value<br>
result <- apply(data, 1, calc.tau)<br>`
