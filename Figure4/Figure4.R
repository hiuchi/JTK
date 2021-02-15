#This is the source code to reproduce Figure 4.

require(tidyverse)
require(poweRlaw)
require(ImpulseDE2)
require(maSigPro)
require(splineTimeR)
require(limma)
require(limorhyde)
require(data.table)
require(foreach)
require(plotROC)
theme_set(theme_bw(base_size = 20))
set.seed(8)

#parameters
gene.num <- 200
replicate.num <- 3
alpha <- 2.5
min <- 30
phi <- 0.05
r <- phi^-1
tp.range <- c(4, 8, 16, 32)
qvalRhyCutoff <- 0.15
qvalDrCutoff <- 0.05

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
  data <- matrix(data = NA, nrow = gene.num, ncol = tp * replicate.num)
  for(i in 1:(gene.num / 4)){
    data[i,] <- f1()
    data[gene.num / 4 + i,] <- f2()
    data[gene.num / 4 * 2 + i,] <- f3()
    data[gene.num / 4 * 3 + i,] <- f4()
  }
  return(data)
}

#function to calculate tau
calc.tau <- function(d){
  control.seq <- as.numeric(d)[(length(d) / 2 + 1):length(d)]
  case.seq <- as.numeric(d)[1:(length(d) / 2)]
  
  tau <- 0
  for(i in 1:length(control.seq)) for(j in 1:(length(control.seq))) if(i < j)
    tau <- tau + (sign(control.seq[j] - control.seq[i]) * sign(case.seq[j] - case.seq[i]))
  tau <- tau / (length(control.seq) * (length(control.seq) - 1) / 2)
  return(tau)
}

#Figure 4
result <- tibble(method = NA, TP = NA, i = NA, j = NA, AUC = NA)
for(tp in tp.range){
  t <- tp
  seq <- seq(0, t - 1, length = tp)
  d <- sim.data()
  d <- replace_na(d, 0)
  rownames(d) <- paste0("gene", 1:gene.num)
  colnames(d) <- paste0("TP", rep(seq, each = replicate.num), "_", rep(1:replicate.num))
  
  for(i in 1:4){
    for(j in 2:4){
      if(i < j){
        data <- cbind(d[((i - 1) * gene.num / 4 + 1):(i * gene.num / 4), ],
                      d[((j - 1) * gene.num / 4 + 1):(j * gene.num / 4), ])
        random.sampling <- sample(c(((i - 1) * gene.num / 4 + 1):(i * gene.num / 4)))
        tmp <- cbind(d[((i - 1) * gene.num / 4 + 1):(i * gene.num / 4), ],
                     d[random.sampling,])
        data <- rbind(data, tmp)
        colnames(data) <- c(paste("control", colnames(data)[1:(ncol(data) / 2)], sep = "_"),
                            paste("case", colnames(data)[(ncol(data) / 2 + 1):ncol(data)], sep = "_"))
        rownames(data) <- paste("gene", 1:nrow(data), sep = "_")
        data.df <- as_tibble(data)
        data.df$gene <- 1:(gene.num / 2)
        data.df <- data.df %>% pivot_longer(colnames(data), names_to = "label", values_to = "counts")
        label <- str_split(data.df$label, "_", simplify = TRUE)
        data.df$sample <- label[,1]
        data.df$TP <- str_sub(label[,2], start = 3) %>% as.numeric
        data.df$replicate <- label[,3]
        data.smr <- data.df %>% group_by(gene, sample, TP) %>% summarise(mean = mean(counts), sd = sd(counts))
        data.smr$label <- paste(data.smr$sample, "TP", data.smr$TP, sep = "_")
        data.mat <- data.smr %>% pivot_wider(id_cols = -c(sample, TP, sd), names_from = "label", values_from = "mean") %>% as.matrix
        rownames(data.mat) <- data.mat[,1]
        data.mat <- data.mat[,-1]
        
        #JTK
        jtk.df <- apply(data.mat, 1, calc.tau)
        df <- data.frame(D = rep(c(1, 0), each = gene.num / 4), M = 1 - as.numeric(jtk.df))
        g <- ggplot(df, aes(d = D, m = M))
        g <- g + geom_roc()
        auc <- calc_auc(g)
        result <- result %>% add_row(method = "JTK",
                                     TP = tp,
                                     i = i,
                                     j = j,
                                     AUC = auc$AUC)
        
        #maSigPro
        timepoints <- paste(rep(seq * 3, each = replicate.num), "h", sep = "")
        header <- paste(timepoints, rep(1:replicate.num), sep = "-")
        
        # create design matrix
        ctrl <- rep(c(1,0), each = length(header))
        mat <- cbind(Time = as.numeric(sub("h","", timepoints)),
                     Replicate = rep(1:(as.integer(tp) * 2), each = as.integer(replicate.num)),
                     Control = ctrl,
                     Treatment = as.numeric(ctrl == 0))
        rownames(mat) <- colnames(data)
        
        # run differential expression analysis
        NBp <- p.vector(data,design = make.design.matrix(mat, degree = 3), counts = TRUE, Q = 1)
        NBt <- T.fit(NBp, step.method = "backward")
        masigpro.df <- NBt$sol[,c(1, 4)]
        df <- data.frame(D = rep(c(1, 0), each = gene.num / 4), M = 1 - as.numeric(masigpro.df$p.valor_TreatmentvsControl))
        g <- ggplot(df, aes(d = D, m = M))
        g <- g + geom_roc()
        auc <- calc_auc(g)
        result <- result %>% add_row(method = "maSigPro",
                                     TP = tp,
                                     i = i,
                                     j = j,
                                     AUC = auc$AUC)
        
        #splineTC
        design <- data.frame(row.names = colnames(data),
                             "SampleName" = colnames(data),
                             "Time" = rep(rep(seq, each = replicate.num), 2), 
                             "Treatment" = rep(c("control", "case"), each = tp * replicate.num),
                             "Replicate" = rep(1:replicate.num, tp * 2))
        phenoData <- new("AnnotatedDataFrame", data = design)
        set <- ExpressionSet(assayData = as.matrix(data), phenoData = phenoData)
        diffExprs <- splineDiffExprs(eSetObject = set, df = 3, reference = c("control", "case")[1], intercept = FALSE)
        splineTC.df <- tibble(gene = rownames(diffExprs), value = diffExprs$P.Value)
        splineTC.df$gene.num <- as.numeric(str_split(splineTC.df$gene, pattern = "_", simplify = TRUE)[,2])
        splineTC.df <- splineTC.df %>% arrange(gene.num)
        df <- data.frame(D = rep(c(1, 0), each = gene.num / 4), M = 1 - as.numeric(splineTC.df$value))
        g <- ggplot(df, aes(d = D, m = M))
        g <- g + geom_roc()
        auc <- calc_auc(g)
        result <- result %>% add_row(method = "splineTC",
                                     TP = tp,
                                     i = i,
                                     j = j,
                                     AUC = auc$AUC)
        
        #LimoRhyde
        sm <- data.frame(title = colnames(data),
                         time = as.numeric(str_sub(str_split(colnames(data), pattern = "_", simplify = TRUE)[,2], 3)),
                         cond = str_split(colnames(data), pattern = "_", simplify = TRUE)[,1])
        sm <- cbind(sm, limorhyde(sm$time, 'time_')) %>% as.data.table
        
        #Identify rhythmic genes
        rhyLimma <- foreach(condNow = unique(sm$cond), .combine = rbind) %do% {
          design <- model.matrix(~ time_cos + time_sin, data = sm[cond == condNow])
          fit <- lmFit(data[, sm$cond == condNow], design)
          fit <- eBayes(fit, trend = TRUE)
          rhyNow <- data.table(topTable(fit, coef = 2:3, number = Inf), keep.rownames = TRUE)
          setnames(rhyNow, 'rn', 'geneId')
          rhyNow[, cond := condNow]
        }
        rhyLimmaSummary <- rhyLimma[, .(P.Value = min(P.Value)), by = geneId]
        rhyLimmaSummary[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
        
        #Identify differentially rhythmic genes
        design <- model.matrix(~ cond * (time_cos + time_sin), data = sm)
        fit <- lmFit(data, design)
        fit <- eBayes(fit, trend = TRUE)
        drLimma <- data.table(topTable(fit, coef = 5:6, number = Inf), keep.rownames = TRUE)
        setnames(drLimma, 'rn', 'geneId')
        drLimma <- drLimma[geneId %in% rhyLimmaSummary[adj.P.Val <= qvalRhyCutoff]$geneId]
        drLimma[, adj.P.Val := p.adjust(P.Value, method = 'BH')]
        drLimma <- drLimma[drLimma$adj.P.Val <= qvalDrCutoff]
        LimoRhyde.df <- tibble(gene = drLimma$geneId, value = drLimma$P.Value)
        
        #Identify differentially expressed genes
        design <- model.matrix(~ cond + time_cos + time_sin, data = sm)
        fit <- lmFit(data, design)
        fit <- eBayes(fit, trend = TRUE)
        deLimma <- data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
        setnames(deLimma, 'rn', 'geneId')
        deLimma <- deLimma[!(geneId %in% drLimma[adj.P.Val <= qvalDrCutoff]$geneId)]
        
        LimoRhyde.df <- LimoRhyde.df %>% add_row(gene = deLimma$geneId, value = deLimma$P.Value)
        LimoRhyde.df$gene.num <- as.numeric(str_split(LimoRhyde.df$gene, pattern = "_", simplify = TRUE)[,2])
        LimoRhyde.df <- LimoRhyde.df %>% arrange(gene.num)
        
        df <- data.frame(D = rep(c(1, 0), each = gene.num / 4), M = 1 - as.numeric(LimoRhyde.df$value))
        g <- ggplot(df, aes(d = D, m = M))
        g <- g + geom_roc()
        auc <- calc_auc(g)
        result <- result %>% add_row(method = "LimoRhyde",
                                     TP = tp,
                                     i = i,
                                     j = j,
                                     AUC = auc$AUC)
        
        #ImpulseDE2
        # specify experimental design
        design <- data.frame("Sample" = colnames(data),
                             "Condition" = rep(c("control", "case"), each = t * replicate.num),
                             "Time" = rep(rep(1:t, each = replicate.num), 2),
                             "Batch" = rep("B_NULL", ncol(data)), 
                             row.names = colnames(data))
        # DEG analysis
        impulse.df <- runImpulseDE2(matCountData = data,
                                    dfAnnotation = design,
                                    boolCaseCtrl = TRUE,
                                    scaNProc = 4,
                                    scaQThres = 1,
                                    boolIdentifyTransients = TRUE)
        df <- data.frame(D = rep(c(1, 0), each = gene.num / 4), M = 1 - as.numeric(impulse.df$dfImpulseDE2Results$p))
        g <- ggplot(df, aes(d = D, m = M))
        g <- g + geom_roc()
        auc <- calc_auc(g)
        result <- result %>% add_row(method = "ImpulseDE2",
                                     TP = tp,
                                     i = i,
                                     j = j,
                                     AUC = auc$AUC)
      }
    }
  }
}
result <- result %>% filter(!is.na(AUC))
result$label <- paste0("Function", result$i, " vs Function", result$j)
result$method <- factor(result$method, levels = c("JTK", "maSigPro", "splineTC", "ImpulseDE2", "LimoRhyde"))
save(result, file = "result.RData")

for(tp in tp.range){
  d <- result %>% dplyr::filter(TP == tp)
  g <- ggplot(d, aes(x = method, y = AUC, group = label))
  g <- g + geom_bar(stat = "identity")
  g <- g + facet_grid(j ~ i) + ggtitle(paste0(tp, " time points"))
  g <- g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(g, file = paste0("./plots/fig4_TP=", tp, ".eps"), dpi = 300)
}
