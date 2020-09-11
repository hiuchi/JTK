require(tidyverse)
require(edgeR)
require(biomaRt)
library(GenomicFeatures)
require(snowfall)
require(maSigPro)
require(splineTimeR)
require(ImpulseDE2)
require(ggVennDiagram)
theme_set(theme_bw(base_size = 20))

#read files
setwd("~/Research/200905_Cardoso/")
mouse.list <- list.files("./results", full.names = TRUE, pattern = "Mouse")
size.files <- file.info(mouse.list)$size
mouse.list <- mouse.list[size.files != 0]
mouse <- read_table2(mouse.list[1], col_names = FALSE)[,1]
mouse <- mouse[1:(nrow(mouse) - 5),]
colnames(mouse) <- "gene"
for(i in mouse.list){
  mouse <- bind_cols(mouse, read_table2(i, col_names = FALSE)[1:nrow(mouse),2])
  name <- str_split(i, "/", simplify = TRUE)[3] %>% str_split("\\.", simplify = TRUE)
  if(length(name) == 8){
    colnames(mouse)[ncol(mouse)] <- name[,3:5] %>% str_flatten(collapse = "_")
  }else{
    colnames(mouse)[ncol(mouse)] <- paste(name[,3:5] %>% str_flatten(collapse = "_"), name[,6], sep = ".")
  }
}
gene.name <- mouse$gene
mouse <- mouse[,-1]
mouse <- as.matrix(mouse)
rownames(mouse) <- gene.name

#calculate rpkm
setwd("~/Research/200911_Cardoso/")
length <- makeTxDbFromGFF('Mus_musculus.GRCm38.101.gtf', format = 'gtf')
length <- exonsBy(length, by = 'gene')
length <- lapply(length, function(x){sum(width(reduce(x)))})
mouse.length <- NULL
for(i in gene.name) mouse.length <- c(mouse.length, length[i])
mouse <- rpkm(mouse, as.numeric(mouse.length))
mouse <- as_tibble(mouse)
mouse$gene <- gene.name

setwd("~/Research/200905_Cardoso/")
rat.list <- list.files("./results", full.names = TRUE, pattern = "Rat")
size.files <- file.info(rat.list)$size
rat.list <- rat.list[size.files != 0]
rat <- read_table2(rat.list[1], col_names = FALSE)[,1]
rat <- rat[1:(nrow(rat) - 5),]
colnames(rat) <- "gene"
for(i in rat.list){
  rat <- bind_cols(rat, read_table2(i, col_names = FALSE)[1:nrow(rat),2])
  name <- str_split(i, "/", simplify = TRUE)[3] %>% str_split("\\.", simplify = TRUE)
  name <- name[,3:5] %>% str_flatten(collapse = "_")
  colnames(rat)[ncol(rat)] <- name
}
gene.name <- rat$gene
rat <- rat[,-1]
rat <- as.matrix(rat)
rownames(rat) <- gene.name

setwd("~/Research/200911_Cardoso/")
length <- makeTxDbFromGFF('Rattus_norvegicus.Rnor_6.0.101.gtf', format = 'gtf')
length <- exonsBy(length, by = 'gene')
length <- lapply(length, function(x){sum(width(reduce(x)))})
rat.length <- NULL
for(i in gene.name) rat.length <- c(rat.length, length[i])
rat <- rpkm(rat, as.numeric(rat.length))
rat <- as_tibble(rat)
rat$gene <- gene.name

genes <- getLDS(attributes = "ensembl_gene_id",
                filters = "ensembl_gene_id", values = mouse$gene,
                mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"), 
                attributesL = "ensembl_gene_id",
                martL = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl"), uniqueRows = T)
genes <- genes %>% distinct(Gene.stable.ID, .keep_all = TRUE)
genes <- genes %>% distinct(Gene.stable.ID.1, .keep_all = TRUE)

convert <- NULL
for(i in rat$gene){
  gene.name <- genes[,1][genes[,2] == i][1]
  if(length(gene.name) == 0) gene.name <- NA
  convert <- c(convert, gene.name)
}
rat$gene <- convert
rat <- na.omit(rat)

rpkm <- inner_join(mouse, rat, by = "gene")
rpkm <- rpkm[, c(316, 1:315, 317:ncol(rpkm))]
rpkm <- rpkm[rowSums(rpkm[,-1]) != 0,]

column.names <- str_split(colnames(rpkm), pattern = "[.]", n = Inf, simplify = TRUE)[,1:2]
for(i in 1:nrow(column.names)) if(column.names[i, 2] != "") column.names[i, 1] <- paste(column.names[i, 1], column.names[i, 2], sep = ".")
for(i in column.names){
  if(sum(column.names == i) > 1){
    replicates <- sum(column.names == i)
    column.names[column.names == i] <- paste(column.names[column.names == i], 1:replicates, sep = "_")
  }else{
    column.names[column.names == i] <- paste(column.names[column.names == i], 1, sep = "_")
  }
}
column.names[1] <- "gene"
colnames(rpkm) <- column.names
rpkm <- rpkm %>% pivot_longer(-gene, names_to = "label", values_to = "rpkm")
labels <- rpkm$label %>% str_split(pattern = "_", simplify = TRUE)
rpkm$organism <- labels[,1]
rpkm$organ <- labels[,2]
rpkm$stage <- labels[,3]
rpkm$replicates <- labels[,4]

rpkm$sequence <- NA
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 10.5] <- 0
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 11.5] <- 1
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 12.5] <- 2
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 13.5] <- 3
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 14.5] <- 4
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 15.5] <- 5
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 16.5] <- 6
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == 17.5] <- 7
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == "0dpb"] <- 8
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == "3dpb"] <- 9
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == "2wpb"] <- 10
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == "4wpb"] <- 11
rpkm$sequence[rpkm$organism == "Mouse" & rpkm$stage == "9wpb"] <- 12

rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "11"] <- 0
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "12"] <- 1
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "14"] <- 2
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "15"] <- 3
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "16"] <- 4
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "18"] <- 5
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "19"] <- 6
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "20"] <- 7
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "0dpb"] <- 8
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "3dpb"] <- 9
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "2wpb"] <- 10
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "6wpb"] <- 11
rpkm$sequence[rpkm$organism == "Rat" & rpkm$stage == "16wpb"] <- 12
rpkm <- na.omit(rpkm)
rpkm <- rpkm %>% filter(organ != 'WholeBrain')
rpkm.smr <- rpkm %>% group_by(gene, organism, organ, sequence) %>% summarise(mean = mean(rpkm), sd = sd(rpkm))

save(rpkm, file ="rpkm.RData")
save(rpkm.smr, file ="rpkm.smr.RData")

#ImpulseDE2
setwd("~/Research/200905_Cardoso/")
mouse.counts.list <- list.files("./results", full.names = TRUE, pattern = "Mouse")
size.files <- file.info(mouse.counts.list)$size
mouse.counts.list <- mouse.counts.list[size.files != 0]
mouse.counts <- read_table2(mouse.counts.list[1], col_names = FALSE)[,1]
mouse.counts <- mouse.counts[1:(nrow(mouse.counts) - 5),]
colnames(mouse.counts) <- "gene"
for(i in mouse.counts.list){
  mouse.counts <- bind_cols(mouse.counts, read_table2(i, col_names = FALSE)[1:nrow(mouse.counts),2])
  name <- str_split(i, "/", simplify = TRUE)[3] %>% str_split("\\.", simplify = TRUE)
  if(length(name) == 8){
    colnames(mouse.counts)[ncol(mouse.counts)] <- name[,3:5] %>% str_flatten(collapse = "_")
  }else{
    colnames(mouse.counts)[ncol(mouse.counts)] <- paste(name[,3:5] %>% str_flatten(collapse = "_"), name[,6], sep = ".")
  }
}
gene.name <- mouse.counts$gene
mouse.counts <- mouse.counts[,-1]
mouse.counts <- as_tibble(mouse.counts)
mouse.counts$gene <- gene.name

rat.counts.list <- list.files("./results", full.names = TRUE, pattern = "Rat")
size.files <- file.info(rat.counts.list)$size
rat.counts.list <- rat.counts.list[size.files != 0]
rat.counts <- read_table2(rat.counts.list[1], col_names = FALSE)[,1]
rat.counts <- rat.counts[1:(nrow(rat.counts) - 5),]
colnames(rat.counts) <- "gene"
for(i in rat.counts.list){
  rat.counts <- bind_cols(rat.counts, read_table2(i, col_names = FALSE)[1:nrow(rat.counts),2])
  name <- str_split(i, "/", simplify = TRUE)[3] %>% str_split("\\.", simplify = TRUE)
  name <- name[,3:5] %>% str_flatten(collapse = "_")
  colnames(rat.counts)[ncol(rat.counts)] <- name
}
gene.name <- rat.counts$gene
rat.counts <- rat.counts[,-1]
rat.counts <- as_tibble(rat.counts)
rat.counts$gene <- gene.name
rat.counts$gene <- convert
rat.counts <- na.omit(rat.counts)

counts <- inner_join(mouse.counts, rat.counts, by = "gene")
counts <- counts[, c(316, 1:315, 317:ncol(counts))]
counts <- counts[rowSums(counts[,-1]) != 0,]

column.names <- str_split(colnames(counts), pattern = "[.]", n = Inf, simplify = TRUE)[,1:2]
for(i in 1:nrow(column.names)) if(column.names[i, 2] != "") column.names[i, 1] <- paste(column.names[i, 1], column.names[i, 2], sep = ".")
for(i in column.names){
  if(sum(column.names == i) > 1){
    replicates <- sum(column.names == i)
    column.names[column.names == i] <- paste(column.names[column.names == i], 1:replicates, sep = "_")
  }else{
    column.names[column.names == i] <- paste(column.names[column.names == i], 1, sep = "_")
  }
}
column.names[1] <- "gene"
colnames(counts) <- column.names
counts <- counts %>% pivot_longer(-gene, names_to = "label", values_to = "counts")
labels <- counts$label %>% str_split(pattern = "_", simplify = TRUE)
counts$organism <- labels[,1]
counts$organ <- labels[,2]
counts$stage <- labels[,3]
counts$replicates <- labels[,4]

counts$sequence <- NA
counts$sequence[counts$organism == "Mouse" & counts$stage == 10.5] <- 0
counts$sequence[counts$organism == "Mouse" & counts$stage == 11.5] <- 1
counts$sequence[counts$organism == "Mouse" & counts$stage == 12.5] <- 2
counts$sequence[counts$organism == "Mouse" & counts$stage == 13.5] <- 3
counts$sequence[counts$organism == "Mouse" & counts$stage == 14.5] <- 4
counts$sequence[counts$organism == "Mouse" & counts$stage == 15.5] <- 5
counts$sequence[counts$organism == "Mouse" & counts$stage == 16.5] <- 6
counts$sequence[counts$organism == "Mouse" & counts$stage == 17.5] <- 7
counts$sequence[counts$organism == "Mouse" & counts$stage == "0dpb"] <- 8
counts$sequence[counts$organism == "Mouse" & counts$stage == "3dpb"] <- 9
counts$sequence[counts$organism == "Mouse" & counts$stage == "2wpb"] <- 10
counts$sequence[counts$organism == "Mouse" & counts$stage == "4wpb"] <- 11
counts$sequence[counts$organism == "Mouse" & counts$stage == "9wpb"] <- 12

counts$sequence[counts$organism == "Rat" & counts$stage == "11"] <- 0
counts$sequence[counts$organism == "Rat" & counts$stage == "12"] <- 1
counts$sequence[counts$organism == "Rat" & counts$stage == "14"] <- 2
counts$sequence[counts$organism == "Rat" & counts$stage == "15"] <- 3
counts$sequence[counts$organism == "Rat" & counts$stage == "16"] <- 4
counts$sequence[counts$organism == "Rat" & counts$stage == "18"] <- 5
counts$sequence[counts$organism == "Rat" & counts$stage == "19"] <- 6
counts$sequence[counts$organism == "Rat" & counts$stage == "20"] <- 7
counts$sequence[counts$organism == "Rat" & counts$stage == "0dpb"] <- 8
counts$sequence[counts$organism == "Rat" & counts$stage == "3dpb"] <- 9
counts$sequence[counts$organism == "Rat" & counts$stage == "2wpb"] <- 10
counts$sequence[counts$organism == "Rat" & counts$stage == "6wpb"] <- 11
counts$sequence[counts$organism == "Rat" & counts$stage == "16wpb"] <- 12
counts <- na.omit(counts)
counts <- counts %>% filter(organ != 'WholeBrain')
setwd("~/Research/200911_Cardoso/")
save(counts, file ="counts.RData")

#function and variables
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

set.seed(8)
cpus <- 4
result <- tibble(organ = NA, gene = NA, method = NA, p = NA, q = NA)
minimum.tp <- 10
for(org in rpkm.smr$organ %>% unique){
  print(org)
  d <- rpkm.smr %>% filter(organ == org)
  m <- filter(d, organism == "Mouse")
  r <- filter(d, organism == "Rat")
  d <- d %>% filter(sequence %in% intersect(unique(r$sequence), unique(m$sequence)))
  d$label <- paste(d$organism, d$sequence, sep = "_")
  d <- d %>% pivot_wider(id_cols = -c(organism, organ, sequence, sd),
                         names_from = label, values_from = mean)
  mat <- matrix(NA, nrow = nrow(d), ncol = ncol(d) - 1)
  mat <- as.matrix(d[,-1])
  colnames(mat) <- colnames(d)[-1]
  rownames(mat) <- d$gene
  mat.seq <- NULL
  for(i in 1:nrow(mat)){
    seq1 <- mat[i, 1:(ncol(mat) / 2)]
    seq2 <- mat[i, (ncol(mat) / 2 + 1):ncol(mat)]
    if(sum(seq1 != 0) >= minimum.tp & sum(seq2 != 0) >= minimum.tp){
      mat.seq <- c(mat.seq, TRUE)
    }else{
      mat.seq <- c(mat.seq, FALSE)
    }
  }
  mat <- mat[mat.seq,]
  
  #JTK
  sfInit(parallel = TRUE, cpus = cpus)
  sfLibrary(base)
  sfLibrary(tidyverse)
  sfExportAll()
  jtk <- sfSapply(rownames(mat), calc.JTK)
  sfRemoveAll()
  sfStop()
  result <- result %>% add_row(organ = org,
                               gene = rownames(mat),
                               method = "JTK",
                               p = as.numeric(jtk),
                               q = p.adjust(jtk, method = "BH"))
  
  for(parameter in c(3, 5, 7)){
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
                    design = make.design.matrix(design, degree = parameter),
                    counts = FALSE,
                    Q = 0.05)
    NBt <- T.fit(NBp, alfa = 1, step.method = "backward")
    result <- result %>% add_row(organ = org,
                                 gene = rownames(NBt$sol),
                                 method = paste0("maSigPro(degree=", parameter, ")"),
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
    diffExprs <- splineDiffExprs(eSetObject = d, df = parameter, reference = c("rat", "mouse")[1], intercept = TRUE)
    result <- result %>% add_row(organ = org,
                                 gene = rownames(diffExprs),
                                 method = paste0("splineTC(df=", parameter, ")"),
                                 p = diffExprs$P.Value,
                                 q = p.adjust(diffExprs$P.Value, method = "BH"))
  }
  
  #ImpuseDE2
  d <- counts %>% filter(organ == org)
  m <- filter(d, organism == "Mouse")
  r <- filter(d, organism == "Rat")
  d <- d %>% filter(sequence %in% intersect(unique(m$sequence), unique(r$sequence)))
  d$label <- paste(d$organism, d$sequence, d$replicates, sep = "_")
  d <- d %>% pivot_wider(id_cols = -c(organism, organ, stage, replicates, sequence),
                         names_from = label, values_from = counts)
  mat <- matrix(NA, nrow = nrow(d), ncol = ncol(d) - 1)
  mat <- as.matrix(d[,-1])
  colnames(mat) <- colnames(d)[-1]
  rownames(mat) <- d$gene
  mat <- mat[mat.seq,]
  
  # specify experimental design
  labels <- str_split(colnames(mat), "_", simplify = TRUE)
  labels[labels == "Mouse"] <- "control"
  labels[labels == "Rat"] <- "case"
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
                              scaNProc = cpus,
                              scaQThres = 1,
                              boolIdentifyTransients = TRUE)
  result <- result %>% add_row(organ = org,
                               gene = impulse.df$dfImpulseDE2Results$Gene,
                               method = "ImpulseDE2",
                               p = impulse.df$dfImpulseDE2Results$p,
                               q = p.adjust(impulse.df$dfImpulseDE2Results$p, method = "BH"))
}
result <- result %>% filter(!is.na(q))
save(result, file = "result.RData")

for(org in result$organ %>% unique){
  data <- result %>% filter(organ == org)
  jtk <- data %>% filter(method == "JTK")
  jtk <- jtk$gene[jtk$q < 0.01]
  masigpro <- data %>% filter(method == "maSigPro(degree=3)")
  masigpro <- masigpro$gene[masigpro$q < 0.01]
  splinetc <- data %>% filter(method == "splineTC(df=3)")
  splinetc <- splinetc$gene[splinetc$q < 0.01]
  impulsede2 <- data %>% filter(method == "ImpulseDE2")
  impulsede2 <- impulsede2$gene[impulsede2$q < 0.05]
  
  d <- list(JTK = jtk, maSigPro = masigpro, splineTC = splinetc, ImpulseDE2 = impulsede2)
  g <- ggVennDiagram(d, label = "count") + scale_fill_gradient(low = "white", high = "white", guide = FALSE)
  g <- g + ggtitle(paste0(org, " (", nrow(data %>% filter(method == "JTK")), " genes)"))
  g <- g + xlim(0, NA)
  ggsave(g, file = paste0("./plots/venndiagram/", org, ".eps"))
  
  jtk.unique <- setdiff(jtk, c(masigpro, splinetc))
  dir.create(paste0("./plots/", org))
  
  for(i in jtk.unique){
    d <- rpkm.smr %>% filter(organ == org, gene == i)
    g <- ggplot(d, aes(x = sequence, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (rpkm)")
    g <- g + guides(colour=FALSE)
    ggsave(g, file = paste0("./plots/", org, "/", i, ".eps"))
  }
  
  pc <- intersect(jtk, c(masigpro, splinetc, impulsede2))
  for(i in pc){
    d <- rpkm.smr %>% filter(organ == org, gene == i)
    g <- ggplot(d, aes(x = sequence, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (rpkm)")
    g <- g + guides(colour=FALSE)
    ggsave(g, file = paste0("./plots/pc/", i, ".eps"))
  }
}

for(org in result$organ %>% unique){
  df <- tibble(method = NA, threshold = NA, value = NA)
  for(threshold in c(0.05, 0.005, 0.0005)){
    data <- result %>% filter(organ == org)
    jtk <- data %>% filter(method == "JTK")
    total <- nrow(jtk)
    jtk <- jtk$gene[jtk$q < threshold]
    masigpro <- data %>% filter(method == "maSigPro(degree=3)")
    masigpro <- masigpro$gene[masigpro$q < threshold]
    masigpro.5 <- data %>% filter(method == "maSigPro(degree=5)")
    masigpro.5 <- masigpro.5$gene[masigpro.5$q < threshold]
    masigpro.7 <- data %>% filter(method == "maSigPro(degree=7)")
    masigpro.7 <- masigpro.7$gene[masigpro.7$q < threshold]
    splinetc <- data %>% filter(method == "splineTC(df=3)")
    splinetc <- splinetc$gene[splinetc$q < threshold]
    splinetc.5 <- data %>% filter(method == "splineTC(df=5)")
    splinetc.5 <- splinetc.5$gene[splinetc.5$q < threshold]
    splinetc.7 <- data %>% filter(method == "splineTC(df=7)")
    splinetc.7 <- splinetc.7$gene[splinetc.7$q < threshold]
    impulsede2 <- data %>% filter(method == "ImpulseDE2")
    impulsede2 <- impulsede2$gene[impulsede2$q < threshold]
    
    df <- df %>% add_row(method = c("JTK", "maSigPro(degree=3)", "maSigPro(degree=5)", "maSigPro(degree=7)",
                                    "splineTC(df=3)", "splineTC(df=5)", "splineTC(df=7)", "ImpulseDE2"),
                         threshold = threshold,
                         value = c(length(jtk) / total * 100, length(masigpro) / total * 100, length(masigpro.5) / total * 100, length(masigpro.7) / total * 100,
                                   length(splinetc) / total * 100, length(splinetc.5) / total * 100, length(splinetc.7) / total * 100, length(impulsede2) / total * 100))
  }
  df <- na.omit(df)
  df$method <- factor(df$method, levels = c("JTK", "maSigPro(degree=3)", "maSigPro(degree=5)", "maSigPro(degree=7)",
                                            "splineTC(df=3)", "splineTC(df=5)", "splineTC(df=7)", "ImpulseDE2"))
  df$threshold <- factor(df$threshold, levels = c("0.05", "0.005", "5e-04"))
  
  g <- ggplot(df, aes(x = threshold, y = value, fill = method))
  g <- g + geom_bar(stat = "identity", position = "dodge") + theme(legend.position="none")
  g <- g + ylab("Persentage of \nsignificant genes (%)") + ggtitle(org)
  g
  ggsave(g, file = paste0("./plots/significant/", org, ".eps"), dpi = 300)
}
