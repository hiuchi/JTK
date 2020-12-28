#This is the source code to reproduce Figures 5, 6 and 7.

require(tidyverse)
require(biomaRt)
require(snowfall)
require(maSigPro)
require(splineTimeR)
require(limorhyde)
require(data.table)
require(ImpulseDE2)
require(gplots)
theme_set(theme_bw(base_size = 20))
set.seed(8)

#loading data and formatting
#mouse
mm <- read_tsv("./files/E-MTAB-6798.sdrf.txt")
mm <- mm %>% dplyr::select(c("Comment[ENA_RUN]", "Characteristics[organism part]", "Characteristics[developmental stage]", "Characteristics[age]"))
mm$`Characteristics[developmental stage]`[mm$`Characteristics[developmental stage]` == "embryo"] <- "e"
mm$`Characteristics[developmental stage]`[mm$`Characteristics[developmental stage]` == "postnatal"] <- "P"
mm$label <- paste0(mm$`Characteristics[organism part]`, "_", mm$`Characteristics[developmental stage]`, mm$`Characteristics[age]`)
mm$replicates <- NA
for(i in unique(mm$label)){
  count <- sum(mm$label == i)
  for(j in 1:count){
    mm$replicates[mm$label == i][j] <- j
  }
}
mm$label <- paste(mm$label, mm$replicates, sep = "_")

files <- list.files("./TPM/mouse", full.names = TRUE)
tmp <- read_tsv(files[1])
mouse <- tibble(gene_id = str_split(tmp$gene_id, pattern = "_", simplify = TRUE)[,1])
for(i in files){
  print(i)
  tmp <- read_tsv(i) %>% dplyr::select(TPM)
  title <- str_split(i, pattern = "/", simplify = TRUE)[,4]
  title <- str_split(title, pattern = "[.]", simplify = TRUE)[,1]
  title <- mm$label[mm$`Comment[ENA_RUN]` == title]
  colnames(tmp) <- paste("mouse", title, sep = "_")
  mouse <- bind_cols(mouse, tmp)
}

#rat
rn <- read_tsv("./files/E-MTAB-6811.sdrf.txt")
rn <- rn %>% dplyr::select(c("Comment[ENA_RUN]", "Characteristics[organism part]", "Characteristics[developmental stage]", "Characteristics[age]"))
rn$`Characteristics[developmental stage]`[rn$`Characteristics[developmental stage]` == "embryo"] <- "e"
rn$`Characteristics[developmental stage]`[rn$`Characteristics[developmental stage]` == "postnatal"] <- "P"
rn$label <- paste0(rn$`Characteristics[organism part]`, "_", rn$`Characteristics[developmental stage]`, rn$`Characteristics[age]`)
rn$replicates <- NA
for(i in unique(rn$label)){
  count <- sum(rn$label == i)
  for(j in 1:count){
    rn$replicates[rn$label == i][j] <- j
  }
}
rn$label <- paste(rn$label, rn$replicates, sep = "_")

files <- list.files("./TPM/rat", full.names = TRUE)
tmp <- read_tsv(files[1])
rat <- tibble(gene_id = str_split(tmp$gene_id, pattern = "_", simplify = TRUE)[,1])
for(i in files){
  print(i)
  tmp <- read_tsv(i) %>% dplyr::select(TPM)
  title <- str_split(i, pattern = "/", simplify = TRUE)[,4]
  title <- str_split(title, pattern = "[.]", simplify = TRUE)[,1]
  title <- rn$label[rn$`Comment[ENA_RUN]` == title]
  colnames(tmp) <- paste("rat", title, sep = "_")
  rat <- bind_cols(rat, tmp)
}

genes <- getLDS(attributes = "ensembl_gene_id",
                filters = "ensembl_gene_id", values = rat$gene_id,
                mart = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl"), 
                attributesL = "ensembl_gene_id",
                martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), uniqueRows = T)
genes <- genes %>% distinct(Gene.stable.ID, .keep_all = TRUE)
genes <- genes %>% distinct(Gene.stable.ID.1, .keep_all = TRUE)

for(i in 1:length(rat$gene_id)){
  gene <- rat$gene_id[i]
  if(sum(genes[,1] == gene) == 0){
    rat$gene_id[i] <- NA
  }else{
    rat$gene_id[i] <- genes$Gene.stable.ID.1[genes[,1] == gene]
  }
}
rat <- drop_na(rat)

tpm <- inner_join(mouse, rat, by = "gene_id")
tpm <- tpm %>% gather(2:ncol(tpm), key = label, value = expression)
title <- str_split(tpm$label, pattern = "_", simplify = TRUE)
tpm$organism <- title[,1]
tpm$organ <- title[,2]
tpm$stage <- title[,3]
tpm$replicates <- title[,4]

tpm$stage.num <- NA
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e10.5"] <- 0
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e11.5"] <- 1
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e12.5"] <- 2
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e13.5"] <- 3
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e14.5"] <- 4
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e15.5"] <- 5
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e16.5"] <- 6
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "e17.5"] <- 7
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "P0"] <- 8
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "P3"] <- 9
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "P14"] <- 10
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "P28"] <- 11
tpm$stage.num[tpm$organism == "mouse" & tpm$stage == "P63"] <- 12

tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e11"] <- 0
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e12"] <- 1
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e14"] <- 2
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e15"] <- 3
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e16"] <- 4
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e18"] <- 5
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e19"] <- 6
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "e20"] <- 7
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "P0"] <- 8
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "P3"] <- 9
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "P14"] <- 10
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "P42"] <- 11
tpm$stage.num[tpm$organism == "rat" & tpm$stage == "P112"] <- 12
tpm <- na.omit(tpm)
tpm.smr <- tpm %>% group_by(gene_id, organism, organ, stage.num) %>% summarise(mean = mean(expression), sd = sd(expression))
save(tpm.smr, file = "tpm.smr.RData")
rm(tpm, tmp, title)

#preparation of counts data
files <- list.files("./counts/mouse", full.names = TRUE)
tmp <- read_tsv(files[1], col_names = FALSE)
mouse <- tibble(gene_id = tmp$X1)
for(i in files){
  print(i)
  tmp <- read_tsv(i, col_names = FALSE)[,2]
  title <- str_split(i, pattern = "/", simplify = TRUE)[,4]
  title <- str_split(title, pattern = "[.]", simplify = TRUE)[,1]
  title <- mm$label[mm$`Comment[ENA_RUN]` == title]
  colnames(tmp) <- paste("mouse", title, sep = "_")
  mouse <- bind_cols(mouse, tmp)
}

#rat
files <- list.files("./counts/rat", full.names = TRUE)
tmp <- read_tsv(files[1], col_names = FALSE)
rat <- tibble(gene_id = tmp$X1)
for(i in files){
  print(i)
  tmp <- read_tsv(i, col_names = FALSE)[,2]
  title <- str_split(i, pattern = "/", simplify = TRUE)[,4]
  title <- str_split(title, pattern = "[.]", simplify = TRUE)[,1]
  title <- rn$label[rn$`Comment[ENA_RUN]` == title]
  colnames(tmp) <- paste("rat", title, sep = "_")
  rat <- bind_cols(rat, tmp)
}

for(i in 1:length(rat$gene_id)){
  gene <- rat$gene_id[i]
  if(sum(genes[,1] == gene) == 0){
    rat$gene_id[i] <- NA
  }else{
    rat$gene_id[i] <- genes$Gene.stable.ID.1[genes[,1] == gene]
  }
}
rat <- drop_na(rat)

counts <- inner_join(mouse, rat, by = "gene_id")
counts <- counts %>% gather(2:ncol(counts), key = label, value = expression)
title <- str_split(counts$label, pattern = "_", simplify = TRUE)
counts$organism <- title[,1]
counts$organ <- title[,2]
counts$stage <- title[,3]
counts$replicates <- title[,4]

counts$stage.num <- NA
counts$stage.num[counts$organism == "mouse" & counts$stage == "e10.5"] <- 0
counts$stage.num[counts$organism == "mouse" & counts$stage == "e11.5"] <- 1
counts$stage.num[counts$organism == "mouse" & counts$stage == "e12.5"] <- 2
counts$stage.num[counts$organism == "mouse" & counts$stage == "e13.5"] <- 3
counts$stage.num[counts$organism == "mouse" & counts$stage == "e14.5"] <- 4
counts$stage.num[counts$organism == "mouse" & counts$stage == "e15.5"] <- 5
counts$stage.num[counts$organism == "mouse" & counts$stage == "e16.5"] <- 6
counts$stage.num[counts$organism == "mouse" & counts$stage == "e17.5"] <- 7
counts$stage.num[counts$organism == "mouse" & counts$stage == "P0"] <- 8
counts$stage.num[counts$organism == "mouse" & counts$stage == "P3"] <- 9
counts$stage.num[counts$organism == "mouse" & counts$stage == "P14"] <- 10
counts$stage.num[counts$organism == "mouse" & counts$stage == "P28"] <- 11
counts$stage.num[counts$organism == "mouse" & counts$stage == "P63"] <- 12

counts$stage.num[counts$organism == "rat" & counts$stage == "e11"] <- 0
counts$stage.num[counts$organism == "rat" & counts$stage == "e12"] <- 1
counts$stage.num[counts$organism == "rat" & counts$stage == "e14"] <- 2
counts$stage.num[counts$organism == "rat" & counts$stage == "e15"] <- 3
counts$stage.num[counts$organism == "rat" & counts$stage == "e16"] <- 4
counts$stage.num[counts$organism == "rat" & counts$stage == "e18"] <- 5
counts$stage.num[counts$organism == "rat" & counts$stage == "e19"] <- 6
counts$stage.num[counts$organism == "rat" & counts$stage == "e20"] <- 7
counts$stage.num[counts$organism == "rat" & counts$stage == "P0"] <- 8
counts$stage.num[counts$organism == "rat" & counts$stage == "P3"] <- 9
counts$stage.num[counts$organism == "rat" & counts$stage == "P14"] <- 10
counts$stage.num[counts$organism == "rat" & counts$stage == "P42"] <- 11
counts$stage.num[counts$organism == "rat" & counts$stage == "P112"] <- 12
counts <- na.omit(counts)
rm(tmp, title)
save(counts, file = "counts.RData")

#function and variables for JTK
trial <- 1000
calc.null.dist <- function(dummy, d){
  seq1 <- d[1:(length(d) / 2)][sample(1:(length(d) / 2), replace = FALSE)] %>% as.numeric
  seq2 <- d[(length(d) / 2 + 1):length(d)][sample(1:(length(d) / 2), replace = FALSE)] %>% as.numeric
  
  tau <- 0
  for(i in 1:length(seq1)) for(j in 1:(length(seq1))) if(i < j)
    tau <- tau + (sign(seq1[j] - seq1[i]) * sign(seq2[j] - seq2[i]))
  tau <- tau / (length(seq1) * (length(seq1) - 1) / 2)
  return(tau)
}
calc.p <- function(d){
  rn.seq <- as.numeric(d)[(length(d) / 2 + 1):length(d)]
  mm.seq <- as.numeric(d)[1:(length(d) / 2)]
  
  tau <- 0
  for(i in 1:length(rn.seq)) for(j in 1:(length(rn.seq))) if(i < j)
    tau <- tau + (sign(rn.seq[j] - rn.seq[i]) * sign(mm.seq[j] - mm.seq[i]))
  tau <- tau / (length(rn.seq) * (length(rn.seq) - 1) / 2)
  
  dist <- sapply(1:trial, calc.null.dist, d)
  return(sum(dist < tau) / trial)
}
calc.JTK <- function(gene.name){
  d <- mat[rownames(mat) == gene.name,]
  return(calc.p(d))
}

result <- tibble(organ = NA, gene_id = NA, method = NA, p = NA, q = NA)
for(org in c("forebrain", "heart", "hindbrain", "kidney", "liver", "ovary", "testis")){
  print(org)
  if(org %in% c("forebrain", "hindbrain")){
    d <- tpm.smr %>% filter(organ %in% c(org, "brain"))
    d <- d %>% filter(!(organism == "rat" & organ == "brain" & stage.num %in% c(2, 3, 4)))
  }else{
    d <- tpm.smr %>% filter(organ == org)
  }
  m <- filter(d, organism == "mouse")
  r <- filter(d, organism == "rat")
  seq <- c(unique(r$stage.num), unique(m$stage.num))[duplicated(c(unique(r$stage.num), unique(m$stage.num)))]
  d <- d %>% filter(stage.num %in% seq)
  d$label <- paste(d$organism, d$stage.num, sep = "_")
  d <- d %>% pivot_wider(id_cols = -c(organism, organ, stage.num, sd),
                         names_from = label, values_from = mean)
  mat <- matrix(NA, nrow = nrow(d), ncol = ncol(d) - 1)
  mat <- as.matrix(d[,-1])
  colnames(mat) <- colnames(d)[-1]
  rownames(mat) <- d$gene_id
  mat <- mat[rowMeans(mat) >= 1,]
  gene.list <- tibble(gene_id = rownames(mat), value = mat[,1])
  
  #JTK
  sfInit(parallel = TRUE, cpus = 4)
  sfLibrary(base)
  sfLibrary(tidyverse)
  sfExportAll()
  jtk <- sfSapply(rownames(mat), calc.JTK)
  sfRemoveAll()
  sfStop()
  result <- result %>% add_row(organ = org,
                               gene_id= rownames(mat),
                               method = "JTK",
                               p = as.numeric(jtk),
                               q = p.adjust(jtk, method = "BH"))
  for(param in c(3, 5)){
    if(param < (ncol(mat) / 2)){
      #maSigPro
      tmp <- str_split(colnames(mat), "_", simplify = TRUE)
      timepoints <- as.numeric(tmp[,2])
      
      # create design matrix
      ctrl <- rep(c(1,0), each = length(timepoints) / 2)
      design <- cbind(Time = timepoints, Replicate = 1, Control = ctrl,
                      Treatment = as.numeric(ctrl == 0))
      rownames(design) <- colnames(mat)
      
      # run differential expression analysis
      NBp <- p.vector(mat,
                      design = make.design.matrix(design, degree = param),
                      counts = FALSE,
                      Q = 0.05)
      NBt <- T.fit(NBp, alfa = 1, step.method = "backward")
      result <- result %>% add_row(organ = org,
                                   gene_id= rownames(NBt$sol),
                                   method = paste0("maSigPro_degree=", param),
                                   p = NBt$sol$p.valor_TreatmentvsControl,
                                   q = p.adjust(NBt$sol$p.valor_TreatmentvsControl, method = "BH"))
      
      #splineTC
      design <- data.frame(row.names = colnames(mat),
                           "SampleName" = colnames(mat),
                           "Time" = timepoints,
                           "Treatment" = rep(c("rat", "mouse"), each = ncol(mat) / 2),
                           "Replicate" = rep(1, ncol(mat)))
      phenoData <- new("AnnotatedDataFrame", data = design)
      d <- ExpressionSet(assayData = mat, phenoData = phenoData)
      diffExprs <- splineDiffExprs(eSetObject = d, df = param, reference = c("rat", "mouse")[1], intercept = TRUE)
      result <- result %>% add_row(organ = org,
                                   gene_id= rownames(diffExprs),
                                   method = paste0("splineTC_df=", param),
                                   p = diffExprs$P.Value,
                                   q = p.adjust(diffExprs$P.Value, method = "BH"))
    }
  }
  
  #LimoRhyde
  sm <- data.frame(title = colnames(mat),
                   time = as.numeric(str_split(colnames(mat), pattern = "_", simplify = TRUE)[,2]),
                   cond = str_split(colnames(mat), pattern = "_", simplify = TRUE)[,1])
  sm <- cbind(sm, limorhyde(sm$time, 'time_'))
  
  design <- model.matrix(~ cond + time_cos + time_sin, data = sm)
  fit <- lmFit(mat, design)
  fit <- eBayes(fit, trend = TRUE)
  deLimma <- data.table(topTable(fit, coef = 2, number = Inf), keep.rownames = TRUE)
  result <- result %>% add_row(organ = org,
                               gene_id = deLimma$rn,
                               method = "LimoRhyde",
                               p = deLimma$P.Value,
                               q = p.adjust(deLimma$P.Value, method = "BH"))
  
  #ImpulseDE2
  if(org %in% c("forebrain", "hindbrain")){
    d <- counts %>% filter(organ %in% c(org, "brain"))
    d <- d %>% filter(!(organism == "rat" & organ == "brain" & stage.num %in% c(2, 3, 4)))
  }else{
    d <- counts %>% filter(organ == org)
  }
  r <- filter(d, organism == "rat")
  m <- filter(d, organism == "mouse")
  seq <- c(unique(r$stage.num), unique(m$stage.num))[duplicated(c(unique(r$stage.num), unique(m$stage.num)))]
  d <- d %>% filter(stage.num %in% seq)
  d$label <- paste(d$organism, d$stage.num, d$replicates, sep = "_")
  d <- d %>% pivot_wider(id_cols = -c(organism, organ, stage, replicates, stage.num),
                         names_from = label, values_from = expression)
  test <- inner_join(d, gene.list, by = "gene_id")
  d <- test[,1:(ncol(test) - 1)]
  mat <- matrix(NA, nrow = nrow(d), ncol = ncol(d) - 1)
  mat <- as.matrix(d[,-1])
  colnames(mat) <- colnames(d)[-1]
  rownames(mat) <- d$gene_id
  
  # specify experimental design
  labels <- str_split(colnames(mat), "_", simplify = TRUE)
  labels[labels == "mouse"] <- "control"
  labels[labels == "rat"] <- "case"
  design <- data.frame("Sample" = colnames(mat),
                       "Condition" = labels[,1],
                       "Time" = labels[,2] %>% as.numeric,
                       "Batch" = rep("B_NULL", ncol(mat)),
                       row.names = colnames(mat))
  # DEG analysis
  impulse.df <- NULL
  impulse.df <- runImpulseDE2(matCountData = mat,
                              dfAnnotation = design,
                              boolCaseCtrl = TRUE,
                              scaNProc = 4,
                              scaQThres = 1,
                              boolIdentifyTransients = TRUE)
  result <- result %>% add_row(organ = org,
                               gene_id = str_sub(impulse.df$dfImpulseDE2Results$Gene),
                               method = "ImpulseDE2",
                               p = impulse.df$dfImpulseDE2Results$p,
                               q = p.adjust(impulse.df$dfImpulseDE2Results$p, method = "BH"))
}
result <- result %>% filter(!is.na(q))
save(result, file = "result.RData")

#draw venn diagram
dir.create("./plots/venndiagram/")
for(org in c("forebrain", "heart", "hindbrain", "kidney", "liver", "ovary", "testis")){
  data <- result %>% filter(organ == org)
  jtk <- data %>% filter(method == "JTK")
  jtk <- jtk$gene_id[jtk$q < 0.05]
  masigpro <- data %>% filter(method == "maSigPro_degree=5")
  masigpro <- masigpro$gene_id[masigpro$q < 0.05]
  splinetc <- data %>% filter(method == "splineTC_df=5")
  splinetc <- splinetc$gene_id[splinetc$q < 0.05]
  limoryde <- data %>% filter(method == "LimoRhyde")
  limoryde <- limoryde$gene_id[limoryde$q < 0.05]
  impulsede2 <- data %>% filter(method == "ImpulseDE2")
  impulsede2 <- impulsede2$gene_id[impulsede2$q < 0.05]
  
  d <- list(JTK = jtk, maSigPro = masigpro, splineTC = splinetc, LimoRhyde = limoryde, ImpulseDE2 = impulsede2)
  postscript(paste0("./plots/venndiagram/", org, ".eps"), title = org, width = 8, height = 8)
  venn(d)
  dev.off()
}

#draw individual plots
dir.create("./plots/JTK_unique/")
dir.create("./plots/pc/")
dir.create("./plots/except_JTK/")
for(org in c("forebrain", "heart", "hindbrain", "kidney", "liver", "ovary", "testis")){
  data <- result %>% filter(organ == org)
  jtk <- data %>% filter(method == "JTK")
  jtk <- jtk$gene_id[jtk$q < 0.05]
  masigpro <- data %>% filter(method == "maSigPro_degree=5")
  masigpro <- masigpro$gene_id[masigpro$q < 0.05]
  splinetc <- data %>% filter(method == "splineTC_df=5")
  splinetc <- splinetc$gene_id[splinetc$q < 0.05]
  limorhyde <- data %>% filter(method == "LimoRhyde")
  limorhyde <- limorhyde$gene_id[limorhyde$q < 0.05]
  impulsede2 <- data %>% filter(method == "ImpulseDE2")
  impulsede2 <- impulsede2$gene_id[impulsede2$q < 0.05]
  
  jtk.unique <- setdiff(jtk, c(masigpro, splinetc, limorhyde))
  dir.create(paste0("./plots/JTK_unique/", org))
  for(i in jtk.unique){
    d <- tpm.smr %>% filter(organ == org, gene_id == i)
    g <- ggplot(d, aes(x = stage.num, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (TPM)")
    g <- g + guides(colour = FALSE) + xlim(-0.5, 12.5) + ylim(0, NA)
    ggsave(g, file = paste0("./plots/JTK_unique/", org, "/", i, ".eps"))
  }
  
  pc <- intersect(jtk, c(masigpro, splinetc, limorhyde, impulsede2))
  dir.create(paste0("./plots/pc/", org))
  for(i in pc){
    d <- tpm.smr %>% filter(organ == org, gene_id == i)
    g <- ggplot(d, aes(x = stage.num, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (TPM)")
    g <- g + guides(colour=FALSE) + xlim(-0.5, 12.5) + ylim(0, NA)
    ggsave(g, file = paste0("./plots/pc/", org, "/", i, ".eps"))
  }
  
  jtk <- data %>% filter(method == "JTK")
  jtk <- jtk$gene_id[jtk$q > 0.05]
  d <- intersect(jtk, c(masigpro, splinetc, limorhyde, impulsede2))
  d <- sample(d, 3)
  dir.create(paste0("./plots/except_JTK/", org))
  for(i in d){
    d <- tpm.smr %>% filter(organ == org, gene_id == i)
    g <- ggplot(d, aes(x = stage.num, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (TPM)")
    g <- g + guides(colour=FALSE) + xlim(-0.5, 12.5) + ylim(0, NA)
    ggsave(g, file = paste0("./plots/except_JTK/", org, "/", i, ".eps"))
  }
}
