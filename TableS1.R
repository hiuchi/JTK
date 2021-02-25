#This is the source code to reproduce Table S1.

require(tidyverse)
require(gprofiler2)

load("./result.RData")
load("./tpm.smr.RData")

table <- NULL
for(org in unique(result$organ)){
  data <- result %>% filter(organ == org, method == "LimoRhyde", q < 0.05)
  
  #backgraound
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
  bg <- rownames(mat)
  
  gostres <- gost(query = data$gene_id,
                  organism = "mmusculus", ordered_query = FALSE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "fdr", 
                  domain_scope = "annotated", custom_bg = bg, 
                  numeric_ns = "", sources = "GO", as_short_link = FALSE)
  if(!is.null(gostres)){
    gostres <- as.data.frame(gostres$result) %>% mutate(organ = org) %>% select(organ, source, term_id, term_name, p_value)
    gostres$p_value <- sprintf('%.2e', gostres$p_value)
    gostres$source <- str_split(gostres$source, pattern = ":", simplify = TRUE)[,2]
    colnames(gostres) <- c("Organ", "Category", "Term ID", "Term name", "P value")
    table <- rbind(table, gostres)
  }
}
write.csv(table, quote = FALSE, row.names = FALSE, file = "enrich.csv")
