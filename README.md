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
```
The example dataset (data.RData) is [here](/data.RData).
```
#load the data
load("data.RData")
data
##      control_TP0 control_TP1 control_TP2 control_TP3 control_TP4 control_TP5
##gene1    88.33333         106    73.00000    57.00000    33.00000          30
##gene2    36.66667          37    71.66667    81.33333    99.66667          79
##      control_TP6 control_TP7 case_TP0 case_TP1 case_TP2 case_TP3 case_TP4
##gene1    18.33333    17.66667 154.6667  81.0000 117.3333       65 60.33333
##gene2   113.00000    83.66667 206.0000 188.3333 126.3333      129 58.33333
##      case_TP5 case_TP6  case_TP7
##gene1       24       21 17.666667
##gene2       23        2  6.666667
```
![plot](https://user-images.githubusercontent.com/2022405/99774217-d96b5480-2b50-11eb-8243-ee6f146dd57e.png)
```
#calculate tau and p-value
result <- apply(data, 1, calc.tau)
result
##        gene1      gene2
##tau 0.8571429 -0.7142857
##p   0.9980000  0.0050000
```
