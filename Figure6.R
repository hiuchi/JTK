#This is the source code to reproduce Figure 6.

require(tidyverse)
require(poweRlaw)
require(ImpulseDE2)
require(maSigPro)
require(splineTimeR)
require(plotROC)
require(tictoc)
theme_set(theme_bw(base_size = 20))
set.seed(8)

#parameters
gene.num <- 100
replicate.num <- 3
alpha <- 2.5
min <- 30
tp <- 8
t <- tp
phi <- 0.05
trial <- 1000

#functions for simulation
f1 <- function(){
  a <- rplcon(1, min, alpha)
  b <- rplcon(1, min, alpha)
  mu <- a * log(seq + 1) + b
  mu[mu < 0] <- 0
  d <- NULL
  for(i in 1:t){
    p <- 1 / (1 + mu[i] * phi)
    d <- c(d, rnbinom(n = replicate.num, size = r, prob = p))
  }
  return(d)
}
f2 <- function(){
  a <- rplcon(1, min, alpha)
  b <- rplcon(1, min, alpha)
  mu <- a * (10 / (seq + 1)) + b
  mu[mu < 0] <- 0
  d <- NULL
  for(i in 1:t){
    p <- 1 / (1 + mu[i] * phi)
    d <- c(d, rnbinom(n = replicate.num, size = r, prob = p))
  }
  return(d)
}
f3 <- function(){
  a <- rplcon(1, min, alpha)
  b <- rplcon(1, min, alpha)
  mu <- a * cos(seq / 2) + 10 + b
  mu[mu < 0] <- 0
  d <- NULL
  for(i in 1:t){
    p <- 1 / (1 + mu[i] * phi)
    d <- c(d, rnbinom(n = replicate.num, size = r, prob = p))
  }
  return(d)
}
f4 <- function(){
  a <- rplcon(1, min, alpha)
  b <- rplcon(1, min, alpha)
  mu <- a * sin(seq / 2) + 10 + b
  mu[mu < 0] <- 0
  d <- NULL
  for(i in 1:t){
    p <- 1 / (1 + mu[i] * phi)
    d <- c(d, rnbinom(n = replicate.num, size = r, prob = p))
  }
  return(d)
}
sim.data <- function(){
  data <- matrix(data = NA, nrow = gene.num * 2, ncol = tp * replicate.num)
  for(i in 1:(gene.num / 4)){
    data[i,] <- f1()
    data[i + gene.num,] <- f1()
    data[gene.num / 4 + i,] <- f2()
    data[gene.num / 4 + i + gene.num,] <- f2()
    data[gene.num / 4 * 2 + i,] <- f3()
    data[gene.num / 4 * 2 + i + gene.num,] <- f3()
    data[gene.num / 4 * 3 + i,] <- f4()
    data[gene.num / 4 * 3 + i + gene.num,] <- f4()
  }
  return(data)
}
generate.teg <- function(data){
  first <- c(1:(gene.num / 4), (1:(gene.num / 4) + (gene.num / 4) * 4))
  second <- first + gene.num / 4
  third <- second + gene.num / 4
  fourth <- third + gene.num / 4
  teg.num <- nrow(data) / 2
  
  sequence <- NULL
  for(i in 1:teg.num){
    if(i <= teg.num / 4) sequence <- c(sequence, sample(setdiff(1:(2 * gene.num), first), 1))
    else if(sum(i == second) == 1) sequence <- c(sequence, sample(setdiff(1:(2 * gene.num), second), 1))
    else if(sum(i == third) == 1) sequence <- c(sequence, sample(setdiff(1:(2 * gene.num), third), 1))
    else if(sum(i == fourth) == 1) sequence <- c(sequence, sample(setdiff(1:(2 * gene.num), fourth), 1))
  }
  sequence <- c(sequence, (gene.num + 1):(gene.num * 2))
  data <- data[sequence,]
  return(data)
}

#functions for JTK
calc.JTK <- function(d){
  control.seq <- as.numeric(d)[(length(d) / 2 + 1):length(d)]
  case.seq <- as.numeric(d)[1:(length(d) / 2)]
  
  tau <- 0
  for(i in 1:length(control.seq)) for(j in 1:(length(control.seq))) if(i < j)
    tau <- tau + (sign(control.seq[j] - control.seq[i]) * sign(case.seq[j] - case.seq[i]))
  tau <- tau / (length(control.seq) * (length(control.seq) - 1) / 2)
  
  dist <- sapply(1:trial, calc.p, d)
  return(sum(dist < tau) / trial)
}
calc.p <- function(dummy, d){
  d <- sample(d, length(d))
  seq1 <- d[1:(length(d) / 2)]
  seq2 <- d[(length(d) / 2 + 1):length(d)]
  
  tau <- 0
  for(i in 1:length(seq1)) for(j in 1:(length(seq1))) if(i < j)
    tau <- tau + (sign(seq1[j] - seq1[i]) * sign(seq2[j] - seq2[i]))
  tau <- tau / (length(seq1) * (length(seq1) - 1) / 2)
  return(tau)
}

r <- phi^-1
t <- tp
seq <- seq(0, t - 1, length = tp)
data <- cbind(sim.data(), generate.teg(sim.data()))
data <- replace_na(data, 0)
rownames(data) <- paste0("gene", 1:(gene.num * 2))
colnames(data) <- c(paste0("control_TP", rep(seq, each = replicate.num), "_", rep(1:replicate.num)),
                    paste0("case_TP", rep(seq, each = replicate.num), "_", rep(1:replicate.num)))
data.df <- as_tibble(data)
data.df$gene <- 1:(gene.num * 2)
data.df <- data.df %>% pivot_longer(colnames(data), names_to = "label", values_to = "counts")
label <- str_split(data.df$label, "_", simplify = TRUE)
data.df$sample <- label[,1]
data.df$TP <- str_sub(label[,2], start = 3) %>% as.numeric
data.df$replicate <- label[,3]
data.df$sample <- factor(data.df$sample, levels = c("control", "case"))
data.smr <- data.df %>% group_by(gene, sample, TP) %>% summarise(mean = mean(counts), sd = sd(counts))
data.smr$label <- paste0(data.smr$sample, "_TP", data.smr$TP)
data.mat <- data.smr %>% pivot_wider(id_cols = -c(sample, TP, sd), names_from = "label", values_from = "mean") %>% as.matrix
rownames(data.mat) <- data.mat[,1]
data.mat <- data.mat[,-1]

result <- tibble(method = NA, param = NA, reps = NA, time = NA)
#result <- tibble(method = NA, param = NA, time = NA)
for(rep in 1:3){
  #JTK
  tic()
  apply(data.mat, 1, calc.JTK)
  tmp <- toc()
  result <- result %>% add_row(method = "JTK",
                               param = NA,
                               reps = rep,
                               time = tmp$toc - tmp$tic)
  
  for(df in c(3, 5, 7)){
    #maSigPro
    timepoints <- paste(rep((seq - 1) * 3, each = replicate.num), "h", sep = "")
    header <- paste(timepoints, rep(1:replicate.num), sep = "-")
    
    # create design matrix
    ctrl <- rep(c(1,0), each = length(header))
    mat <- cbind(Time = as.numeric(sub("h","", timepoints)),
                 Replicate = rep(1:(as.integer(tp) * 2), each = as.integer(replicate.num)),
                 Control = ctrl,
                 Treatment = as.numeric(ctrl == 0))
    rownames(mat) <- colnames(data)
    
    # run differential expression analysis
    tic()
    NBp <- p.vector(data,design = make.design.matrix(mat, degree = df), counts = TRUE, Q = 1)
    NBt <- T.fit(NBp, step.method = "backward")
    tmp <- toc()
    result <- result %>% add_row(method = "maSigPro",
                                 param = df,
                                 reps = rep,
                                 time = tmp$toc - tmp$tic)
    
    #splineTC
    design <- data.frame(row.names = colnames(data),
                         "SampleName" = colnames(data),
                         "Time" = rep(rep(seq, each = replicate.num), 2), 
                         "Treatment" = rep(c("control", "case"), each = tp * replicate.num),
                         "Replicate" = rep(1:replicate.num, tp * 2))
    phenoData <- new("AnnotatedDataFrame", data = design)
    d <- ExpressionSet(assayData = as.matrix(data), phenoData = phenoData)
    tic()
    diffExprs <- splineDiffExprs(eSetObject = d, df = df, reference = c("control", "case")[1], intercept = FALSE)
    tmp <- toc()
    result <- result %>% add_row(method = "splineTC",
                                 param = df,
                                 reps = rep,
                                 time = tmp$toc - tmp$tic)
  }
  
  #ImpulseDE2
  # specify experimental design
  design <- data.frame("Sample" = colnames(data),
                       "Condition" = rep(c("control", "case"), each = t * replicate.num),
                       "Time" = rep(rep(1:t, each = replicate.num), 2),
                       "Batch" = rep("B_NULL", ncol(data)), 
                       row.names = colnames(data))
  # DEG analysis
  tic()
  runImpulseDE2(matCountData = data,
                dfAnnotation = design,
                boolCaseCtrl = TRUE,
                scaNProc = 1,
                scaQThres = 1,
                boolIdentifyTransients = TRUE)
  tmp <- toc()
  result <- result %>% add_row(method = "ImpulseDE2",
                               param = NA,
                               reps = rep,
                               time = tmp$toc - tmp$tic)
}

result <- result %>% filter(!is.na(time))
result$method <- factor(result$method, levels = c("JTK", "maSigPro", "splineTC", "ImpulseDE2"))
result$param <- as.factor(result$param)
save(result, file = "result.RData")

result$label <- paste(result$method, result$param, sep = "_")
df <- result %>% group_by(label) %>% summarise(elapsed_time = mean(time), sd = sd(time))
df$label <- c("ImpulseDE2", "JTK", "maSigPro (dgree=3)", "maSigPro (dgree=5)", "maSigPro (dgree=7)",
              "splineTC (df=3)", "splineTC (df=5)", "splineTC (df=7)")
df$label <- fct_relevel(df$label, "JTK", "maSigPro (dgree=3)", "maSigPro (dgree=5)", "maSigPro (dgree=7)",
                        "splineTC (df=3)", "splineTC (df=5)", "splineTC (df=7)", "ImpulseDE2")

theme_set(theme_bw(base_size = 20))
g <- ggplot(df, aes(x = label, y = elapsed_time))
g <- g + geom_point() + xlab(NULL) + ylab("Elapsed time (sec, log)")
g <- g + geom_errorbar(aes(ymin = elapsed_time - sd, ymax = elapsed_time + sd))
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1))
g <- g + scale_y_log10()
g
ggsave(g, file = "time.eps", dpi = 300)

