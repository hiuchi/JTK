#This is the scripts to reproduce Figure 4A and 4B.

require(tidyverse)
require(snowfall)
require(maSigPro)
require(splineTimeR)
require(ImpulseDE2)
require(stringr)
require(plotROC)
theme_set(theme_bw(base_size = 20))
set.seed(8)

#parameters
trial <- 1000
file <- "./data/30M_3rep_4TP_contr.txt"

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
    d <- sample(d, length(d)) %>% as.numeric
    seq1 <- d[1:(length(d) / 2)]
    seq2 <- d[(length(d) / 2 + 1):length(d)]
    
    tau <- 0
    for(i in 1:length(seq1)) for(j in 1:(length(seq1))) if(i < j)
        tau <- tau + (sign(seq1[j] - seq1[i]) * sign(seq2[j] - seq2[i]))
    tau <- tau / (length(seq1) * (length(seq1) - 1) / 2)
    return(tau)
}

# parse parameters
par <- strsplit(basename(file), "_")[[1]]
repl <- as.numeric(strsplit(par[2], "rep")[[1]])
TP <- as.numeric(strsplit(par[3], "TP")[[1]])
time <- ((1:TP) - 1) * 3
header <- paste(rep(((1:TP) - 1) * 3, each = repl), "h-", rep(1:repl), sep = "")

# read data
contr <- read.table(file, header = T, stringsAsFactors = F)
treat <- read.table(sub("contr", "treat", file), header = T, stringsAsFactors = F)

# merge data into 1 data.framea
d <- cbind(contr, treat)
type <- c("control", "case")
colnames(d) <- paste(rep(type, each = length(header)), rep(header, 2), sep = "_")
d <- as.matrix(d[rowSums(d) > 0,]) # otherwise crash

#JTK
data <- as.data.frame(d)
rm(d)
data$gene <- rownames(data)
data <- data %>% pivot_longer(cols = 1:(ncol(data) - 1), names_to = "sample", values_to = "expression")
design <- str_split(data$sample, "_", simplify = TRUE)
data$cc <- design[,1]
data$TP <- str_split(design[,2], "-", simplify = TRUE)[,1]
data$rep <- str_split(design[,2], "-", simplify = TRUE)[,2]
data$sample <- str_split(data$sample, "-", simplify = TRUE)[,1]
data <- data %>% group_by(gene, sample, cc, TP) %>% summarise(expression = mean(expression))
data$TP <- factor(data$TP, levels = c("0h", "3h", "6h", "9h", "12h", "15h", "18h", "21h", "24h", "27h", "30h", "33h"))
data.mat <- data %>% pivot_wider(id_cols = -c(TP, cc), names_from = "sample", values_from = "expression") %>% as.matrix
rownames(data.mat) <- data.mat[,1]
data.mat <- data.mat[,-1]

#JTK
jtk.df <- apply(data.mat, 1, calc.JTK)
JTK <- tibble(gene = str_sub(names(jtk.df), start = 6) %>% as.numeric,
              method = "JTK",
              value = as.numeric(jtk.df))

#splineTC
inData <- cbind(contr, treat)
type <- c("contr", "treat")
colnames(inData) <- paste(rep(type, each = length(header)), rep(header, 2), sep = "_")
inData <- as.matrix(inData[rowSums(inData)>0,]) # otherwise crash

design <- data.frame(row.names = colnames(inData),
                     "SampleName" = colnames(inData),
                     "Time" = rep(rep(time, each = repl), 2), 
                     "Treatment" = rep(type, each = TP * repl),
                     "Replicate" = rep(1:repl, TP * 2))
phenoData <- new("AnnotatedDataFrame", data = design)
data <- ExpressionSet(assayData = as.matrix(inData), phenoData = phenoData)
diffExprs <- splineDiffExprs(eSetObject = data, df = 3, reference = type[1], intercept = TRUE)
splineTC.df <- diffExprs[,c("P.Value", "adj.P.Val")]
splineTC <- tibble(gene = str_sub(rownames(splineTC.df), start = 6) %>% as.numeric,
                   method = "splineTC",
                   value = splineTC.df$P.Value)

#ImpulseDE2
# merge data into 1 data.framea
rownames(treat) <- rownames(contr)
inData <- cbind(contr, treat)
type <- c("control", "case")
colnames(inData) <- paste(rep(type, each = length(header)), rep(header, 2), sep = "_")
inData <- as.matrix(inData[rowSums(inData)>0,]) # otherwise crash

# specify experimental design
design <- data.frame("Sample" = colnames(inData),
                     "Condition" = rep(type,each = TP * repl),
                     "Time" = rep(rep(time,each = repl), 2),
                     "Batch" = rep("B_NULL", ncol(inData)), 
                     row.names = colnames(inData))
# DEG analysis
impulse.df <- runImpulseDE2(matCountData = inData,
                            dfAnnotation = design,
                            boolCaseCtrl = TRUE,
                            scaNProc = 4,
                            scaQThres = 1,
                            boolIdentifyTransients = TRUE)

ImpulseDE2 <- tibble(gene = str_sub(impulse.df$dfImpulseDE2Results$Gene, start = 6) %>% as.numeric,
                     method = "ImpulseDE2",
                     value = impulse.df$dfImpulseDE2Results$p)

#maSigPro
# merge data frames
data <- merge(contr, treat, by = "row.names")
rownames(data) <- data$Row.names
data$Row.names <- NULL
colnames(data) <- paste(header, rep(c("contr","treat"), each = as.numeric(TP) * as.numeric(repl)), sep = "_")

# create design matrix
ctrl <- rep(c(1,0), each = length(header))
timepoints <- paste(rep(((1:TP)-1)*3,each=repl),"h",sep="")
mat <- cbind(Time = as.numeric(sub("h","", timepoints)), Replicate = rep(1:(as.integer(TP)*2), each = as.integer(repl)), Control = ctrl, Treatment = as.numeric(ctrl == 0))
rownames(mat) <- colnames(data)

# run differential expression analysis
NBp <- p.vector(data, design = make.design.matrix(mat, degree = 3), counts = TRUE, Q = 1)
NBt <- T.fit(NBp, step.method = "backward")

get <- get.siggenes(NBt, vars="groups")
diff <- get$sig.genes$TreatmentvsControl$sig.pvalues[,c(1,2)]

masigpro.df <- NBt$sol[,c(1, 4)]
maSigPro <- tibble(gene = str_sub(rownames(masigpro.df), start = 6) %>% as.numeric,
                   method = "maSigPro",
                   value = masigpro.df$`p-value`)

result <- rbind(JTK, maSigPro, splineTC, ImpulseDE2)

degs <- read.table("./DEG_IDs.txt")
ids <- read.table("./SIM_IDs.txt")
find.deg <- function(gene.name) return(as.character(ids$V3[which(ids$V2 == gene.name)]))
deg.ids <- sapply(as.character(degs$V1), find.deg)

result$DEG <- TRUE
for(i in 1:nrow(result)) if(result$gene[i] > 1200) result$DEG[i] <- FALSE
result$method <- factor(result$method, levels = c("JTK", "maSigPro", "splineTC", "ImpulseDE2"))
save(result, file = "result.RData")

result.roc <- tibble(D = result$DEG, method = result$method, value = 1 - result$value)
g <- ggplot(result.roc, aes(d = D, m = value, colour = method))
g <- g + geom_roc(n.cuts = FALSE)
g <- g + xlab("False positive rate") + ylab("True positive rate")
g <- g + theme(legend.title = element_text(size=15),legend.text = element_text(size=15)) 
g
ggsave(g, file = ".Figure4B.eps")

#plots
theme_set(theme_bw(base_size = 15))
# parse parameters
par <- strsplit(basename(file), "_")[[1]]
repl <- as.numeric(strsplit(par[2], "rep")[[1]])
TP <- as.numeric(strsplit(par[3], "TP")[[1]])
time <- ((1:TP) - 1) * 3
header <- paste(rep(((1:TP) - 1) * 3, each = repl), "h-", rep(1:repl), sep = "")

# read data
contr <- read.table(file, header = T, stringsAsFactors = F)
treat <- read.table(sub("contr", "treat", file), header = T, stringsAsFactors = F)
rownames(treat) <- rownames(contr)

# merge data into 1 data.framea
d <- cbind(contr, treat)
type <- c("control", "case")
colnames(d) <- paste(rep(type, each = length(header)), rep(header, 2), sep = "_")
d <- as.matrix(d[rowSums(d) > 0,]) # otherwise crash

data <- as.data.frame(d)
data$gene <- rownames(data) %>% str_sub(start = 6) %>% as.numeric
data <- data %>% pivot_longer(cols = 1:(ncol(data) - 1), names_to = "sample", values_to = "expression")
design <- str_split(data$sample, "_", simplify = TRUE)
data$cc <- design[,1]
data$TP <- str_split(design[,2], "-", simplify = TRUE)[,1] %>% str_sub(end = -2) %>% as.numeric
data$rep <- str_split(design[,2], "-", simplify = TRUE)[,2]
data$sample <- str_split(data$sample, "-", simplify = TRUE)[,1]
data <- data %>% group_by(gene, cc, TP) %>% summarise(mean = mean(expression), sd = sd(expression))
data$cc <- factor(data$cc, levels = c("control", "case"))

jtk <- result %>% filter(method == "JTK")
jtk <- jtk %>% arrange(gene)
masigpro <- result %>% filter(method == "maSigPro")
masigpro$value <- p.adjust(masigpro$value, method = "BH")
masigpro <- masigpro %>% arrange(gene)
splinetc <- result %>% filter(method == "splineTC")
splinetc$value <- p.adjust(splinetc$value, method = "BH")
splinetc <- splinetc %>% arrange(gene)
impulsede2 <- result %>% filter(method == "ImpulseDE2")
impulsede2$value <- p.adjust(impulsede2$value, method = "BH")

genes <- jtk %>% filter(gene <= 1200, value >= 0.8)
genes <- as.vector(genes$gene)
for(i in genes){
    jtk.value <- jtk %>% filter(gene == i)
    masigpro.value <- masigpro %>% filter(gene == i)
    splinetc.value <- splinetc %>% filter(gene == i)
    impulsede2.value <- impulsede2 %>% filter(gene == i)
    title <- paste0("JTK's q=", jtk.value$value %>% signif(2),
                    ", maSigPro's q=", masigpro.value$value %>% signif(2),
                    ",\nsplineTC's q=",  splinetc.value$value %>% signif(2),
                    ", ImpulseDE2's q=", impulsede2.value$value %>% signif(2))
    d <- data %>% filter(gene == i)
    g <- ggplot(d, aes(x = TP, y = mean, colour = cc))
    g <- g + geom_point() + ggtitle(title) + xlab("Time (h)") + ylab("Counts")
    g <- g + geom_line(data = d, aes(x = TP, y = mean)) + guides(colour=FALSE) 
    g <- g + ylim(0, NA)+ geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5)
    g <- g + annotate("text", label = title, x = 2.5, y = 10) 
    ggsave(g, file = paste0("./plots/gene", i, ".pdf"))
}

selected.genes <- c(682, 684, 783, 792, 1171, 1174)
d$cc <- factor(d$cc, levels = c("control", "case"))
for(i in selected.genes){
    jtk.value <- jtk %>% filter(gene == i)
    masigpro.value <- masigpro %>% filter(gene == i)
    splinetc.value <- splinetc %>% filter(gene == i)
    impulsede2.value <- impulsede2 %>% filter(gene == i)
    title <- paste0("JTK's q=", jtk.value$value %>% signif(2),
                    "\nmaSigPro's q=", masigpro.value$value %>% signif(2),
                    "\nsplineTC's q=",  splinetc.value$value %>% signif(2),
                    "\nImpulseDE2's q=", impulsede2.value$value %>% signif(2))
    d <- data %>% filter(gene == i)
    g <- ggplot(d, aes(x = TP, y = mean, colour = cc))
    g <- g + geom_point() + xlab("Time (h)") + ylab("Counts")
    g <- g + geom_line(data = d, aes(x = TP, y = mean))# + guides(colour=FALSE) 
    g <- g + ylim(0, NA)+ geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5)
    g <- g + annotate("text", label = title, x = 7.5, y = 400) 
    ggsave(g, file = paste0("./examples/gene", i, ".pdf"))
}

jtk.value <- jtk %>% filter(gene %in% selected.genes)
masigpro.value <- masigpro %>% filter(gene %in% selected.genes)
splinetc.value <- splinetc %>% filter(gene %in% selected.genes)
impulsede2.value <- impulsede2 %>% filter(gene %in% selected.genes)
title <- paste0("JTK's q=", jtk.value$value %>% signif(2),
                "\nmaSigPro's q=", masigpro.value$value %>% signif(2),
                "\nsplineTC's q=",  splinetc.value$value %>% signif(2),
                "\nImpulseDE2's q=", impulsede2.value$value %>% signif(2))
d <- data %>% filter(gene %in% selected.genes)
d$description <- c(rep("Up early slow", each = 16),
                   rep("Down early slow", each = 16),
                   rep("Mixed slow", each = 16))
d$gene <- rep(c(1, 2), each = 8)

text.df <- tibble(gene = rep(c(1, 2), 3),
                  description = rep(c("Up early slow", "Down early slow", "Mixed slow"), each = 2),
                  text = title, cc = NA)
text.df$description <- factor(text.df$description,
                              levels = c("Up early slow", "Down early slow", "Mixed slow"))
text.df$cc <- factor(text.df$cc, levels = c("control", "case"))
g <- ggplot(d, aes(x = TP, y = mean, colour = cc))
g <- g + geom_point() + xlab("Time (h)") + ylab("Counts")
g <- g + geom_line(data = d, aes(x = TP, y = mean))# + guides(colour=FALSE) 
g <- g + ylim(0, NA)+ geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5)
g <- g + geom_text(data = text.df, mapping = aes(x = 5, y = 600, label = text))
g <- g + facet_grid(gene ~ description)
g
ggsave(g, file = "Figure4A.eps", dpi = 300)
