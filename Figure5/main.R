require(tidyverse)
require(edgeR)
require(biomaRt)
require(snowfall)
require(maSigPro)
require(splineTimeR)
require(ImpulseDE2)
theme_set(theme_bw(base_size = 20))

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

data <- inner_join(mouse, rat, by = "gene")
data <- data[, c(231, 1:230, 232:ncol(data))]
data <- data[rowSums(data[,-1]) != 0,]

column.names <- str_split(colnames(data), pattern = "[.]", n = Inf, simplify = TRUE)[,1:2]
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
colnames(data) <- column.names
data <- data %>% pivot_longer(-gene, names_to = "label", values_to = "CPM")
labels <- data$label %>% str_split(pattern = "_", simplify = TRUE)
data$organism <- labels[,1]
data$organ <- labels[,2]
data$stage <- labels[,3]
data$replicates <- labels[,4]

data$sequence <- NA
data$sequence[data$organism == "Mouse" & data$stage == 10.5] <- 1
data$sequence[data$organism == "Mouse" & data$stage == 11.5] <- 2
data$sequence[data$organism == "Mouse" & data$stage == 12.5] <- 3
data$sequence[data$organism == "Mouse" & data$stage == 13.5] <- 4
data$sequence[data$organism == "Mouse" & data$stage == 14.5] <- 5
data$sequence[data$organism == "Mouse" & data$stage == 15.5] <- 6
data$sequence[data$organism == "Mouse" & data$stage == 16.5] <- 7
data$sequence[data$organism == "Mouse" & data$stage == 17.5] <- 8
data$sequence[data$organism == "Mouse" & data$stage == "0dpb"] <- 9
data$sequence[data$organism == "Mouse" & data$stage == "3dpb"] <- 10
data$sequence[data$organism == "Mouse" & data$stage == "2wpb"] <- 11
data$sequence[data$organism == "Mouse" & data$stage == "4wpb"] <- 12
data$sequence[data$organism == "Mouse" & data$stage == "9wpb"] <- 13

data$sequence[data$organism == "Rat" & data$stage == "11"] <- 1
data$sequence[data$organism == "Rat" & data$stage == "12"] <- 2
data$sequence[data$organism == "Rat" & data$stage == "14"] <- 3
data$sequence[data$organism == "Rat" & data$stage == "15"] <- 4
data$sequence[data$organism == "Rat" & data$stage == "16"] <- 5
data$sequence[data$organism == "Rat" & data$stage == "18"] <- 6
data$sequence[data$organism == "Rat" & data$stage == "19"] <- 7
data$sequence[data$organism == "Rat" & data$stage == "20"] <- 8
data$sequence[data$organism == "Rat" & data$stage == "0dpb"] <- 9
data$sequence[data$organism == "Rat" & data$stage == "3dpb"] <- 10
data$sequence[data$organism == "Rat" & data$stage == "2wpb"] <- 11
data$sequence[data$organism == "Rat" & data$stage == "6wpb"] <- 12
data$sequence[data$organism == "Rat" & data$stage == "16wpb"] <- 13
data <- na.omit(data)
data <- data %>% filter(organ != 'WholeBrain')
data.smr <- data %>% group_by(gene, organism, organ, sequence) %>% summarise(mean = mean(CPM))

save(data, file ="data.RData")
save(data.smr, file ="data.smr.RData")

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
for(org in data$organ %>% unique){
  d <- data.smr %>% filter(organ == org & !is.na(mean))
  r <- filter(d, organism == "Rat")
  m <- filter(d, organism == "Mouse")
  seq <- c(unique(r$sequence), unique(m$sequence))[duplicated(c(unique(r$sequence), unique(m$sequence)))]
  d <- d %>% filter(sequence %in% seq)
  d$label <- paste(d$organism, d$sequence, sep = "_")
  d <- d %>% pivot_wider(id_cols = -c(organism, organ, sequence),
                         names_from = label, values_from = mean)
  mat <- matrix(NA, nrow = nrow(d), ncol = ncol(d) - 1)
  mat <- as.matrix(d[,-1])
  colnames(mat) <- colnames(d)[-1]
  rownames(mat) <- d$gene
  mat.seq <- NULL
  for(i in 1:nrow(mat)){
    seq1 <- mat[i, 1:(ncol(mat) / 2)]
    seq2 <- mat[i, (ncol(mat) / 2 + 1):ncol(mat)]
    if(sum(seq1 == 0) > 7 | sum(seq2 == 0) > 7){
      mat.seq <- c(mat.seq, FALSE)
    }else{
      mat.seq <- c(mat.seq, TRUE)
    }
  }
  mat <- mat[mat.seq,]
  
  if(ncol(mat) != 0){
    #JTK
    sfInit(parallel = TRUE, cpus = 4)
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
                                     gene = rownames(NBt$sol),
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
                                     gene = rownames(diffExprs),
                                     method = paste0("splineTC_df=", param),
                                     p = diffExprs$P.Value,
                                     q = p.adjust(diffExprs$P.Value, method = "BH"))
      }
    }
  }
}
result <- result %>% filter(!is.na(q))
save(result, file = "result.RData")

#ImpulseDE2
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

data <- inner_join(mouse, rat, by = "gene")
data <- data[, c(231, 1:230, 232:ncol(data))]
data <- data[rowSums(data[,-1]) != 0,]

column.names <- str_split(colnames(data), pattern = "[.]", n = Inf, simplify = TRUE)[,1:2]
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
colnames(data) <- column.names
data <- data %>% pivot_longer(-gene, names_to = "label", values_to = "counts")
labels <- data$label %>% str_split(pattern = "_", simplify = TRUE)
data$organism <- labels[,1]
data$organ <- labels[,2]
data$stage <- labels[,3]
data$replicates <- labels[,4]

data$sequence <- NA
data$sequence[data$organism == "Mouse" & data$stage == 10.5] <- 1
data$sequence[data$organism == "Mouse" & data$stage == 11.5] <- 2
data$sequence[data$organism == "Mouse" & data$stage == 12.5] <- 3
data$sequence[data$organism == "Mouse" & data$stage == 13.5] <- 4
data$sequence[data$organism == "Mouse" & data$stage == 14.5] <- 5
data$sequence[data$organism == "Mouse" & data$stage == 15.5] <- 6
data$sequence[data$organism == "Mouse" & data$stage == 16.5] <- 7
data$sequence[data$organism == "Mouse" & data$stage == 17.5] <- 8
data$sequence[data$organism == "Mouse" & data$stage == "0dpb"] <- 9
data$sequence[data$organism == "Mouse" & data$stage == "3dpb"] <- 10
data$sequence[data$organism == "Mouse" & data$stage == "2wpb"] <- 11
data$sequence[data$organism == "Mouse" & data$stage == "4wpb"] <- 12
data$sequence[data$organism == "Mouse" & data$stage == "9wpb"] <- 13

data$sequence[data$organism == "Rat" & data$stage == "11"] <- 1
data$sequence[data$organism == "Rat" & data$stage == "12"] <- 2
data$sequence[data$organism == "Rat" & data$stage == "14"] <- 3
data$sequence[data$organism == "Rat" & data$stage == "15"] <- 4
data$sequence[data$organism == "Rat" & data$stage == "16"] <- 5
data$sequence[data$organism == "Rat" & data$stage == "18"] <- 6
data$sequence[data$organism == "Rat" & data$stage == "19"] <- 7
data$sequence[data$organism == "Rat" & data$stage == "20"] <- 8
data$sequence[data$organism == "Rat" & data$stage == "0dpb"] <- 9
data$sequence[data$organism == "Rat" & data$stage == "3dpb"] <- 10
data$sequence[data$organism == "Rat" & data$stage == "2wpb"] <- 11
data$sequence[data$organism == "Rat" & data$stage == "6wpb"] <- 12
data$sequence[data$organism == "Rat" & data$stage == "16wpb"] <- 13
data <- na.omit(data)
data <- data %>% filter(organ != 'WholeBrain')
data <- data %>% arrange(organism, sequence)
save(data, file ="counts.RData")

impulsede2 <- tibble(organ = NA, gene = NA, method = NA, p = NA, q = NA)
for(org in data$organ %>% unique){
  d <- data %>% filter(organ == org & !is.na(counts))
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
  mat.seq <- NULL
  for(i in 1:nrow(mat)){
    seq1 <- mat[i, 1:(ncol(mat) / 2)]
    seq2 <- mat[i, (ncol(mat) / 2 + 1):ncol(mat)]
    if(sum(seq1 == 0) > 7 | sum(seq2 == 0) > 7){
      mat.seq <- c(mat.seq, FALSE)
    }else{
      mat.seq <- c(mat.seq, TRUE)
    }
  }
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
                              scaNProc = 4,
                              scaQThres = 1,
                              boolIdentifyTransients = TRUE)
  result <- result %>% add_row(organ = org,
                               gene = str_sub(impulse.df$dfImpulseDE2Results$Gene, start = 5),
                               method = "ImpulseDE2",
                               p = impulse.df$dfImpulseDE2Results$p,
                               q = p.adjust(impulse.df$dfImpulseDE2Results$p, method = "BH"))
}
result <- result %>% filter(!is.na(q))
save(result, file = "result.RData")
