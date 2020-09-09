#This is the source code to reproduce Figures 5 and 6.

require(tidyverse)
require(edgeR)
require(biomaRt)
require(snowfall)
require(maSigPro)
require(splineTimeR)
require(ImpulseDE2)
require(ggVennDiagram)
theme_set(theme_bw(base_size = 20))

cpus <- 4

#read files
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
mouse <- cpm(mouse)
mouse <- as_tibble(mouse)
mouse$gene <- gene.name

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
rat <- cpm(rat)
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

cpm <- inner_join(mouse, rat, by = "gene")
cpm <- cpm[, c(316, 1:315, 317:ncol(cpm))]
cpm <- cpm[rowSums(cpm[,-1]) != 0,]

column.names <- str_split(colnames(cpm), pattern = "[.]", n = Inf, simplify = TRUE)[,1:2]
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
colnames(cpm) <- column.names
cpm <- cpm %>% pivot_longer(-gene, names_to = "label", values_to = "CPM")
labels <- cpm$label %>% str_split(pattern = "_", simplify = TRUE)
cpm$organism <- labels[,1]
cpm$organ <- labels[,2]
cpm$stage <- labels[,3]
cpm$replicates <- labels[,4]

cpm$sequence <- NA
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 10.5] <- 0
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 11.5] <- 1
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 12.5] <- 2
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 13.5] <- 3
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 14.5] <- 4
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 15.5] <- 5
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 16.5] <- 6
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == 17.5] <- 7
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == "0dpb"] <- 8
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == "3dpb"] <- 9
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == "2wpb"] <- 10
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == "4wpb"] <- 11
cpm$sequence[cpm$organism == "Mouse" & cpm$stage == "9wpb"] <- 12

cpm$sequence[cpm$organism == "Rat" & cpm$stage == "11"] <- 0
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "12"] <- 1
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "14"] <- 2
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "15"] <- 3
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "16"] <- 4
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "18"] <- 5
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "19"] <- 6
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "20"] <- 7
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "0dpb"] <- 8
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "3dpb"] <- 9
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "2wpb"] <- 10
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "6wpb"] <- 11
cpm$sequence[cpm$organism == "Rat" & cpm$stage == "16wpb"] <- 12
cpm <- na.omit(cpm)
cpm <- cpm %>% filter(organ != 'WholeBrain')
cpm.smr <- cpm %>% group_by(gene, organism, organ, sequence) %>% summarise(mean = mean(CPM), sd = sd(CPM))

save(cpm, file ="cpm.RData")
save(cpm.smr, file ="cpm.smr.RData")

#ImpulseDE2
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
mouse <- as_tibble(mouse)
mouse$gene <- gene.name

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
rat <- as_tibble(rat)
rat$gene <- gene.name
rat$gene <- convert
rat <- na.omit(rat)

counts <- inner_join(mouse, rat, by = "gene")
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
counts <- counts %>% arrange(organism, sequence)
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
result <- tibble(organ = NA, gene = NA, method = NA, p = NA, q = NA)
minimum.tp <- 10
for(org in cpm$organ %>% unique){
  print(org)
  d <- cpm.smr %>% filter(organ == org)
  m <- filter(d, organism == "Mouse")
  r <- filter(d, organism == "Rat")
  seq <- c(unique(r$sequence), unique(m$sequence))[duplicated(c(unique(r$sequence), unique(m$sequence)))]
  d <- d %>% filter(sequence %in% seq)
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
  r <- filter(d, organism == "Rat")
  m <- filter(d, organism == "Mouse")
  seq <- c(unique(r$sequence), unique(m$sequence))[duplicated(c(unique(r$sequence), unique(m$sequence)))]
  d <- d %>% filter(sequence %in% seq)
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
  ggsave(g, file = paste0("./plots/venndiagram/", org, ".eps"))
  
  jtk.unique <- setdiff(jtk, c(masigpro, splinetc))
  dir.create(paste0("./plots/", org))
  
  for(i in jtk.unique){
    d <- cpm.smr %>% filter(organ == org, gene == i)
    g <- ggplot(d, aes(x = sequence, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (CPM)")
    ggsave(g, file = paste0("./plots/", org, "/", i, ".eps"))
  }
  
  pc <- intersect(jtk, c(masigpro, splinetc, impulsede2))
  for(i in pc){
    d <- cpm.smr %>% filter(organ == org, gene == i)
    g <- ggplot(d, aes(x = sequence, y = mean, colour = organism))
    g <- g + geom_point() + geom_line() + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))
    g <- g + ggtitle(paste(org, i, sep = ", ")) + xlab("Stage") + ylab("Expression (CPM)")
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
