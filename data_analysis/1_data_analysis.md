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

