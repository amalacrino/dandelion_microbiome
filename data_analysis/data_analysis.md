Below we report the functions we used for the major steps in our analyses, together with an example of how we used them.

### Calculate diversity index

```r
calculate.div <- function(psobject, compartment){
  ps.div <- prune_samples(sample_data(x)$sample_type == y, x)
  otus <- as.data.frame(otu_table(ps.div))
  diversity <- estimate_richness(ps.div, split = TRUE)
  diversity <- cbind(sample_data(ps.div), diversity)
  return(diversity)
}

div.16S <- calculate.div(ps.16S, y = "rhizosphere")
```

### Bayesian model for diversity

```r
bayes.div <- function(df, metric){
  fit.rst <- rstanarm::stan_glm(
                    get(metric) ~ plant_genotype * treatment,
                    data = df, cores=4, chains=4, seed = 100)
  return(fit.rst)
}

get.p <- function(samples){
  if(mean(samples) > 0){
    p <- sum(samples <= 0)/length(samples)
  }else{
    p <- sum(samples >= 0)/length(samples)
  }
  return(max(p*2, (1/4000)*2))
}

bayes.test.posteriors <- function(x){
  aaa <- as.data.frame(x)
  aaa <- aaa[,c(1:6)]
  colnames(aaa)[1] <- c("Intercept")
  
  colnames(aaa) <- aaa %>% rename_all(~stringr::str_replace(.,"^b_","")) %>% colnames(.) 
  
  aaa$RNAi_herb <- aaa$Intercept
  aaa$RNAi_control <- aaa$Intercept + aaa$treatmentnoherb 
  aaa$RNAi_wounding <- aaa$Intercept + aaa$treatmentwounding 
  aaa$Wt_herb <- aaa$Intercept + aaa$plant_genotypeWt 
  aaa$Wt_control <- aaa$Intercept + aaa$treatmentnoherb + aaa$plant_genotypeWt + aaa$`plant_genotypeWt:treatmentnoherb`
  aaa$Wt_wounding <- aaa$Intercept + aaa$treatmentwounding + aaa$plant_genotypeWt + aaa$`plant_genotypeWt:treatmentwounding`
  
  b.con.h <- hypothesis(aaa, c("Wt_herb = Wt_control",
                               "Wt_wounding = Wt_control",
                               "Wt_herb = Wt_wounding",
                               
                               "RNAi_herb = RNAi_control",
                               "RNAi_wounding = RNAi_control",
                               "RNAi_herb = RNAi_wounding",
                               
                               "Wt_control = RNAi_control",
                               "Wt_herb = RNAi_herb",
                               "Wt_wounding = RNAi_wounding"),
                        alpha = 0.05)
  
  b.con.h <- as.data.frame(b.con.h$hypothesis)
  
  b.con.h$Pmcmc <- c(get.p((aaa$Wt_herb - aaa$Wt_control)),
                     get.p((aaa$Wt_wounding - aaa$Wt_control)),
                     get.p((aaa$Wt_herb - aaa$Wt_wounding)),
                     get.p((aaa$RNAi_herb - aaa$RNAi_control)),
                     get.p((aaa$RNAi_wounding - aaa$RNAi_control)),
                     get.p((aaa$RNAi_herb - aaa$RNAi_wounding)),
                     get.p((aaa$Wt_control - aaa$RNAi_control)),
                     get.p((aaa$Wt_herb - aaa$RNAi_herb)),
                     get.p((aaa$Wt_wounding - aaa$RNAi_wounding)))
  
  makeStars <- function(x){
    stars <- c("****", "***", "**", "*", "ns")
    vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
    i <- findInterval(x, vec)
    stars[i] }
  
  b.con.h$sig <- makeStars(b.con.h$Pmcmc)
  
  b.con.h <- b.con.h[c(1, 9, 10, 2:5)]
  
  return(b.con.h)
}

fit.rst <- bayes.div(div.16S, "Shannon")
summary(fit.rst, probs = c(0.05, 0.95), digits = 3)
bayes.test.posteriors(fit.rst)
```

### Bayesian multivariate analysis

```r
bayes.model <- function(psobject, compartment){
  ps <- prune_samples(sample_data(psobject)$sample_type == compartment, x)
  ps.n <- norm.deseq(ps)
  otus <- as.data.frame(otu_table(ps.n))
  sampledf <- data.frame(sample_data(ps))
  
  data <- otus %>% tibble::rownames_to_column("otu") %>% gather(sample, count, -otu) %>%
  left_join(sampledf %>% rownames_to_column("sample"), by = "sample")
  
  fit1 <- brm(mvbind(count) ~ plant_genotype * treatment, data = data, chains = 4, cores = NTHREADS)
  
  return(fit1)
}

bayes.test.posteriors <- function(x){
  aaa <- as.data.frame(x)
  aaa <- aaa[,c(1:6)]
  colnames(aaa)[1] <- c("Intercept")
  
  colnames(aaa) <- aaa %>% rename_all(~stringr::str_replace(.,"^b_","")) %>% colnames(.) 
  
  aaa$RNAi_herb <- aaa$Intercept
  aaa$RNAi_control <- aaa$Intercept + aaa$treatmentnoherb 
  aaa$RNAi_wounding <- aaa$Intercept + aaa$treatmentwounding 
  aaa$Wt_herb <- aaa$Intercept + aaa$plant_genotypeWt 
  aaa$Wt_control <- aaa$Intercept + aaa$treatmentnoherb + aaa$plant_genotypeWt + aaa$`plant_genotypeWt:treatmentnoherb`
  aaa$Wt_wounding <- aaa$Intercept + aaa$treatmentwounding + aaa$plant_genotypeWt + aaa$`plant_genotypeWt:treatmentwounding`
  
  b.con.h <- hypothesis(aaa, c("Wt_herb = Wt_control",
                               "Wt_wounding = Wt_control",
                               "Wt_herb = Wt_wounding",
                               
                               "RNAi_herb = RNAi_control",
                               "RNAi_wounding = RNAi_control",
                               "RNAi_herb = RNAi_wounding",
                               
                               "Wt_control = RNAi_control",
                               "Wt_herb = RNAi_herb",
                               "Wt_wounding = RNAi_wounding"),
                        alpha = 0.05)
  
  b.con.h <- as.data.frame(b.con.h$hypothesis)
  
  b.con.h$Pmcmc <- c(get.p((aaa$Wt_herb - aaa$Wt_control)),
                     get.p((aaa$Wt_wounding - aaa$Wt_control)),
                     get.p((aaa$Wt_herb - aaa$Wt_wounding)),
                     get.p((aaa$RNAi_herb - aaa$RNAi_control)),
                     get.p((aaa$RNAi_wounding - aaa$RNAi_control)),
                     get.p((aaa$RNAi_herb - aaa$RNAi_wounding)),
                     get.p((aaa$Wt_control - aaa$RNAi_control)),
                     get.p((aaa$Wt_herb - aaa$RNAi_herb)),
                     get.p((aaa$Wt_wounding - aaa$RNAi_wounding)))
  
  makeStars <- function(x){
    stars <- c("****", "***", "**", "*", "ns")
    vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
    i <- findInterval(x, vec)
    stars[i] }
  
  b.con.h$sig <- makeStars(b.con.h$Pmcmc)
  
  b.con.h <- b.con.h[c(1, 9, 10, 2:5)]
  
  return(b.con.h)
}

a <- bayes.model(ps.16S, "rhizosphere")
summary(a, probs = c(0.05, 0.95), digits = 3)
bayes_R2(a)
bayes.test.posteriors(a)
```

### Magnitude of change compared to the control group

The code below selects a specific compartments and compares each treatment to the respective control, returning a list of log2FC values for each OTU for each group.

```r
process.diff.taxa <- function(psobject, compartment, rand){
  ps <- if(rand == "F"){x}else{generate.rand.ps(x)}
  ps.pmv <- prune_samples(sample_data(ps)$sample_type == compartment, ps)
  sampledf <- data.frame(sample_data(ps.pmv))
  GM.diff.taxa <- ps.pmv
  GM.diff.taxa.wt <- subset_samples(GM.diff.taxa, plant_genotype == "Wt")
  GM.diff.taxa.rnai <- subset_samples(GM.diff.taxa, plant_genotype == "RNAi")
  
  dds <- vector("list", length = 4)
  names(dds) <- c("d1", "d2", "d3", "d4")
  
  dds$d1 <- cal.diff.taxa(GM.diff.taxa.wt, "herb", "noherb", factor = "treatment")
  dds$d2 <- cal.diff.taxa(GM.diff.taxa.wt, "wounding", "noherb", factor = "treatment")
  dds$d3 <- cal.diff.taxa(GM.diff.taxa.rnai, "herb", "noherb", factor = "treatment")
  dds$d4 <- cal.diff.taxa(GM.diff.taxa.rnai, "wounding", "noherb", factor = "treatment")
  return(dds)
}

magnitude.treatment <- function(dds){
  d1 <- dds$d1
  d2 <- dds$d2
  d3 <- dds$d3
  d4 <- dds$d4
  
  d1$group <- "Wt_herbivory"
  d2$group <- "Wt_wounding"
  d3$group <- "RNAi_herbivory"
  d4$group <- "RNAi_wounding"
  
  tx <- rbind(d1, d2, d3, d4)
  tx$log2FoldChange <- abs(tx$log2FoldChange)
  return(tx)
}

df.dt.16s.rh <- process.diff.taxa(ps.16S, "rhizosphere", rand = "F")
df.dt.16s.rh2 <- magnitude.treatment(df.dt.16s.rh)

model <- lmer(log2FoldChange ~ group * (1|rn), data = df.dt.its.rh2)
Anova(model)
m1 <- emmeans(model, "group")
pairs(m1)
```

### Identify differentially abundant (frequent) taxa (genes)

```r
list.diff <- function(dataset, psobject, direction, subset){
  dds <- lapply(dataset, selection.diff.taxa, direction = direction)
  '%ni%' <- Negate('%in%')
  
  d1 <- dds$d1
  d2 <- dds$d2
  d3 <- dds$d3
  d4 <- dds$d4
  
  tax.table <- if(deparse(substitute(psobject)) == "ps.genecontent"){
    genes.table}else{
      as.data.frame(tax_table(psobject))[,c(6,7)]
    }
  
  dx <- if(subset == "Wt_herb"){d1[which(d1$rn %ni% d2$rn)]
  } else if(subset == "Wt_wound"){d2[which(d2$rn %ni% d1$rn)]
  } else if(subset == "RNAi_herb"){d3[which(d3$rn %ni% d4$rn)]
  } else {d4[which(d4$rn %ni% d3$rn)]}
  
  tx <- if(deparse(substitute(psobject)) == "ps.genecontent"){
    tax.table[which(tax.table$locus_tag %in% dx$rn),]}else{
      tax.table[which(row.names(tax.table) %in% dx$rn),]
    }
  
  tx <- cbind(dx, tx)

  tx <- if(deparse(substitute(y)) == "ps.genecontent"){
    tx[order(tx$product),]}else{
      tx[order(tx$Genus),]
    }
  return(tx)
}

p2 <- list.diff.asvs(df.dt.16s.rh, ps.16S, "up", "Wt_herb")
```

### Test if the number of differential taxa (genes) are different from random

```r
generate.rand.ps <- function(psobject){
  psx <- psobject
  ps <- as.data.frame(otu_table(psx))
  set.seed(100)
  otu <- randomizeMatrix(ps,null.model = "frequency",iterations = 1000)
  otu_table(psx) <- otu_table(otu, taxa_are_rows = TRUE)
  return(psx)
}

random.test <- function(dataset, datasetRandom, direction){
  dds <- lapply(dataset, selection.diff.taxa, direction = direction)
  dds.rand <- lapply(datasetRandom, selection.diff.taxa, direction = direction)
  
  d1 <- dds[['d1']]$rn
  d2 <- dds[['d2']]$rn
  d3 <- dds[['d3']]$rn
  d4 <- dds[['d4']]$rn
  
  d1.r <- dds.rand[['d1']]$rn
  d2.r <- dds.rand[['d2']]$rn
  d3.r <- dds.rand[['d3']]$rn
  d4.r <- dds.rand[['d4']]$rn
  
  df <- data.frame(Group = c("Wt_herb", "Wt_wound", "RNAi_herb", "RNAi_wound"),
                   Observed_OTUs = c(length(d1), 
                                     length(d2), 
                                     length(d3),
                                     length(d4)),
                   Random_OTUs =   c(length(d1.r), 
                                     length(d2.r), 
                                     length(d3.r),
                                     length(d4.r)))
  
  df <- df %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      chi_sq = if(Observed_OTUs == 0 & Random_OTUs == 0){NA}else{chisq.test(c(Observed_OTUs, Random_OTUs))$statistic},
      p_val = if(Observed_OTUs == 0 & Random_OTUs == 0){NA}else{chisq.test(c(Observed_OTUs, Random_OTUs))$p.value})
  return(df)
}

df.dt.16s.rh <- process.diff.taxa(ps.16S, "rhizosphere", rand = "F")
df.dt.16s.rh.rand <- process.diff.taxa(ps.16S, "rhizosphere", rand = "T")
rnd.test <- random.test(df.dt.16s.rh, df.dt.16s.rh.rand, "none")
```

### Number of OTUs unique to each group

```r
nrotusfun <- function(psobject, compartment, treatment, genotype){
  ks <- sample_data(psobject)[["sample_type"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  ks <- sample_data(ps)[["treatment"]] %in% paste0(treatment)
  ps <- prune_samples(samples = ks, ps)
  ks <- sample_data(ps)[["plant_genotype"]] %in% paste0(genotype)
  ps <- prune_samples(samples = ks, ps)

  ps.nrotus <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
  list1 <- row.names(as.data.frame(otu_table(ps.nrotus)))
  
  ks <- sample_data(psobject)[["sample_type"]] %in% paste0(compartment)
  ps <- prune_samples(samples = ks, psobject)
  ks <- sample_data(ps)[["treatment"]] %in% "noherb"
  ps <- prune_samples(samples = ks, ps)
  ks <- sample_data(ps)[["plant_genotype"]] %in% paste0(genotype)
  ps <- prune_samples(samples = ks, ps)

  ps.nrotus <- filter_taxa(ps, function (x) {sum(x > 0) > 1}, prune=TRUE)
  list2 <- row.names(as.data.frame(otu_table(ps.nrotus)))
  
  list <- list1[!(list1 %in% list2)]
  
  return(list)
}


df <- data.frame(compartment = c("rhizosphere", "rhizosphere", "root", "root"),
                 treatment = c("herbivory", "wounding", "herbivory", "wounding"),
                 NIL = c(length(nrotusfun(ps.16S, "rhizosphere", "herb", "Wt")),
                         length(nrotusfun(ps.16S, "rhizosphere", "wounding", "Wt")),
                         length(nrotusfun(ps.16S, "root", "herb", "Wt")),
                         length(nrotusfun(ps.16S, "root", "wounding", "Wt"))),
                 RNAi = c(length(nrotusfun(ps.16S, "rhizosphere", "herb", "RNAi")),
                          length(nrotusfun(ps.16S, "rhizosphere", "wounding", "RNAi")),
                          length(nrotusfun(ps.16S, "root", "herb", "RNAi")),
                          length(nrotusfun(ps.16S, "root", "wounding", "RNAi"))))
  
df <- df %>%
    dplyr::rowwise() %>% 
    dplyr::mutate(
      chi_sq = chisq.test(c(NIL, RNAi))$statistic,
      p_val = chisq.test(c(NIL, RNAi))$p.value)
```

### Test for differential abundance at the genus level

```r
tax.fun <- function(psobject, compartment, plant_genotype, treatment, plant_genotype2, treatment2, test){
  glom <- microbiome::aggregate_taxa(psobject, "Genus")
  dat <- psmelt(glom)
  dat <- dat %>% dplyr::group_by(Sample, sample_type, plant_genotype, treatment, Genus) %>% dplyr::summarize(cs = mean(Abundance)) %>% dplyr::mutate(cs = cs/sum(cs)) 
  filt.gen <- dat %>% dplyr::group_by(Genus) %>% dplyr::summarize(mean = mean(cs)) %>% dplyr::filter(mean <= 0.01)
  dat <- subset(dat, !(Genus %in% filt.gen$Genus))
  dat <- dat[which(dat$sample_type == sample_type & dat$plant_genotype == plant_genotype & dat$treatment == treatment |
                     dat$sample_type == sample_type & dat$plant_genotype == plant_genotype2 & dat$treatment == treatment2),]
  group.host <- group_by(dat, Genus)
  list.bact <-unique(c(as.character(group.host$Genus)))
  
  model_calculator <- sapply(list.bact,  
                             function(x){
                               data.s <- group.host[which(group.host$Genus==x),]
                               model <- if(test == "plant_genotype"){lm(cs ~ plant_genotype, data = data.s)}else{lm(cs ~ treatment, data = data.s)}
                               aaa <-  Anova(model)
                               aaa$sig = c(rep('',length(aaa$`Pr(>F)`)))
                               makeStars <- function(x){
                                 stars <- c("****", "***", "**", "*", "ns")
                                 vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
                                 i <- findInterval(x, vec)
                                 stars[i] }
                               aaa$sig <- makeStars(aaa$`Pr(>F)`)
                               aaa <- aaa[1,]
                               return(aaa)},
                             simplify = FALSE,USE.NAMES = TRUE)
  res <- do.call(rbind, model_calculator)
  res <- setDT(res, keep.rownames = TRUE)[]
  
  dat2 <- if(test == "plant_genotype"){dat %>% dplyr::group_by(Genus, plant_genotype) %>% dplyr::summarize(cs = mean(cs))}else{dat %>% dplyr::group_by(Genus, treatment) %>% dplyr::summarize(cs = mean(cs))}
  res$group1  <- if(test == "plant_genotype"){dat2[which(dat2$plant_genotype == "Wt"),]$cs}else{dat2[which(dat2$treatment == treatment),]$cs}
  res$group2  <- if(test == "plant_genotype"){dat2[which(dat2$plant_genotype == "RNAi"),]$cs}else{dat2[which(dat2$treatment == treatment2),]$cs}
  
  colnames(res)[7] <- if(test == "plant_genotype"){"Wt"}else{paste0(treatment)}
  colnames(res)[8] <- if(test == "plant_genotype"){"RNAi"}else{paste0(treatment2)}
  
  res <- res[which(res$`Pr(>F)` < 0.05),]
  return(res)
}

df <- tax.fun(ps.16s, "root", "RNAi", "wounding", "RNAi", "noherb", "treatment")
```
