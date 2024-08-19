
# Original function traitDependentBAMM written by Dan Rabosky & Huateng Huang, 2014 in the BAMMtools package

# Additions 
  # Use the names of the tipStates, tipLambda, and tipMu to check compatibility between tips rather than the tip.label in the phylogeny
  # Allows to use the function to compute STRAPP test for any given time
  # Allows to choose if random sampling of posterior configurations must be done with replacement or not

run_STRAPP_test <- function (BAMM_data, trait_data, reps, rate = "speciation", return.full = FALSE, 
                         method = "spearman", logrates = TRUE, two.tailed = TRUE, replace = FALSE,
                         traitorder = NA, nthreads = 1) 
{
  if (nthreads > 1) {
    if (!"package:parallel" %in% search()) {
      stop("Please load package 'parallel' for using the multi-thread option\n")
    }
  }
  if (BAMM_data$type == "trait") {
    rate <- "trait"
  }
  ratetype.option <- c("speciation", "extinction", "net diversification", 
                       "trait")
  ratetype <- grep(paste0("^", rate), ratetype.option, value = TRUE, 
                   ignore.case = TRUE)
  if (length(ratetype) == 0) {
    stop("Rate must be one of 'speciation', 'extinction', or 'net diversification', only the initial letter is needed.\n")
  }
  if (ratetype == "net diversification" & logrates == TRUE) {
    cat("WARNING: Net diversification might be negative and logged rates would then produce NaNs.\n")
  }
  if (ratetype == "trait" & BAMM_data$type == "diversification") {
    stop("Rate must be either speciation or net diversification if BAMM_data is a diversification analysis.\n")
  }
  if (ratetype %in% c("speciation", "net diversification") & 
      BAMM_data$type == "trait") {
    stop("Rate must be trait if BAMM_data is a trait analysis.\n")
  }
  if (!identical(sort(names(trait_data)), sort(BAMM_data$tip.label))) {
    stop("Species names in the BAMM_data object and trait data are not identical. You need to sort them accordingly")
  }
  method.option <- c("spearman", "pearson", "mann-whitney", 
                     "kruskal")
  method <- grep(paste0("^", method), method.option, ignore.case = TRUE, 
                 value = TRUE)
  if (length(method) == 0) {
    stop("method must be one of 'spearman', 'pearson', 'mann-whitney', or 'kruskal', only the initial letter is needed")
  }
  if (method == "spearman" | method == "pearson") {
    if (!is.numeric(trait_data)) {
      cat(paste0("selected ", method, ", but the trait is not numeric, converted the trait into a numeric vector\n"))
      trait_data <- as.numeric(trait_data)
    }
  }
  else if (method == "mann-whitney" | method == "kruskal") {
    if (length(unique(trait_data[!is.na(trait_data)])) == 1) {
      stop(paste("selected ", method, ", but the trait only has one level\n", 
                 sep = ""))
    }
    if (method == "mann-whitney") {
      if (length(unique(trait_data[!is.na(trait_data)])) > 2) {
        stop(paste("selected ", method, ", but the trait has more than two levels\n", 
                   sep = ""))
      }
    }
  }
  trait.state <- NA
  if (!two.tailed) {
    if (anyNA(traitorder)) {
      stop("selected one-tail test, but traitorder is not specified\n")
    }
    if (method == "kruskal") {
      stop(" currently one-tail test is only available for continous or binary trait")
    }
    if (method == "spearman" | method == "pearson") {
      direction.option <- c("positive", "negative")
      direction <- grep(paste0("^", traitorder), direction.option, 
                        ignore.case = TRUE, value = TRUE)
      if (length(direction) == 0) {
        stop(" for one-tail test with continous trait, traitorder must be either 'positive' or 'negative', only the initial letter is needed")
      }
      else {
        cat(paste0("select one-tailed ", method, " test\nAlternative hypothesis: the trait is ", 
                   direction, "ly correlated with speciation rate\n"))
      }
    }
    else {
      traitorder <- gsub(" ", "", traitorder)
      trait.state <- as.character(unlist(strsplit(x = traitorder, 
                                                  split = ",")))
      if (length(trait.state) != 2) {
        stop("please specify the traitorder for binary trait:\nTwo states separated by comma, and the state that is expected to have lower speciation rate first\n")
      }
      else {
        cat(paste0("selected one-tail ", method, " test\nAlternative hypothesis: species with trait ", 
                   trait.state[2], " has higher speciation rate than those with trait ", 
                   trait.state[1], "\n"))
      }
      for (i in trait.state) {
        if (sum(trait_data == i) == 0) {
          stop(paste("no species with state ", i, " \n", 
                     sep = ""))
        }
      }
    }
  }
  if (ratetype == "speciation") {
    tiprates <- BAMM_data$tipLambda
  }
  else if (ratetype == "extinction") {
    tiprates <- BAMM_data$tipMu
  }
  else if (ratetype == "net diversification") {
    tiprates <- lapply(1:length(BAMM_data$tipLambda), function(i) BAMM_data$tipLambda[[i]] - 
                         BAMM_data$tipMu[[i]])
  }
  else if (ratetype == "trait") {
    tiprates <- BAMM_data$tipLambda
  }
  tipstates <- BAMM_data$tipStates
  trait_data <- trait_data[BAMM_data$tip.label]
  stat.mu <- 0
  if (method == "mann-whitney") {
    trait.stat.count <- table(trait_data)
    trait.stat.count <- trait.stat.count[!is.na(names(trait.stat.count))]
    stat.mu <- prod(trait.stat.count)/2
  }
  if (logrates) {
    tiprates <- lapply(1:length(tiprates), function(x) log(tiprates[[x]]))
  }
  gen <- sample(1:length(tiprates), size = reps, replace = replace)
  gen.tiprates <- list()
  for (l in 1:length(gen)) {
    gen.tiprates[[l]] <- data.frame(rates = tiprates[[gen[l]]], 
                                    states = tipstates[[gen[l]]], stringsAsFactors = FALSE)
  }
  rm("tiprates", "tipstates")
  permute_tiprates <- function(m) {
    tt <- m$states
    tlam <- m$rates
    index <- unique(tt)
    lvec <- numeric(length(index))
    for (k in 1:length(index)) {
      lvec[k] <- tlam[tt == index[k]][1]
    }
    new_index <- sample(index, size = length(index), replace = FALSE)
    x <- rep(0, length(tt))
    for (xx in 1:length(index)) {
      x[which(tt == index[xx])] <- lvec[which(index == 
                                                new_index[xx])]
    }
    x
  }
  if (nthreads > 1) {
    cl <- parallel::makePSOCKcluster(nthreads)
    p.gen.tiprates <- parallel::parLapply(cl, gen.tiprates, 
                                          permute_tiprates)
    parallel::stopCluster(cl)
  }
  else {
    p.gen.tiprates <- lapply(gen.tiprates, permute_tiprates)
  }
  xgen.tiprates <- list()
  for (l in 1:length(gen)) {
    xgen.tiprates[[l]] <- gen.tiprates[[l]]$rates
  }
  gen.tiprates <- xgen.tiprates
  rm("xgen.tiprates")
  cortest <- function(rates, trait_data, method) {
    if (sd(rates, na.rm = TRUE) == 0) {
      return(0)
    }
    else {
      return(cor.test(rates, trait_data, method = method, exact = FALSE)$estimate)
    }
  }
  manntest <- function(rates, trait_data, two.tailed, trait.state) {
    if (two.tailed) {
      return(wilcox.test(rates ~ trait_data, exact = FALSE)$statistic)
    }
    else {
      return(wilcox.test(rates[which(trait_data == trait.state[2])], 
                         rates[which(trait_data == trait.state[1])], exact = FALSE)$statistic)
    }
  }
  kruskaltest <- function(rates, trait_data) {
    testres <- kruskal.test(rates ~ trait_data)
    if (is.na(testres$statistic)) {
      return(qchisq(p = 0.999, df = testres$parameter))
    }
    else {
      return(testres$statistic)
    }
  }
  if (nthreads > 1) {
    cl <- parallel::makePSOCKcluster(nthreads)
    if (method == "spearman" | method == "pearson") {
      obs <- parallel::parLapply(cl, gen.tiprates, cortest, 
                                 trait_data, method)
      permu <- parallel::parLapply(cl, p.gen.tiprates, 
                                   cortest, trait_data, method)
    }
    else if (method == "mann-whitney") {
      obs <- parallel::parLapply(cl, gen.tiprates, manntest, 
                                 trait_data, two.tailed, trait.state)
      permu <- parallel::parLapply(cl, p.gen.tiprates, 
                                   manntest, trait_data, two.tailed, trait.state)
    }
    else {
      obs <- parallel::parLapply(cl, gen.tiprates, kruskaltest, 
                                 trait_data)
      permu <- parallel::parLapply(cl, p.gen.tiprates, 
                                   kruskaltest, trait_data)
    }
    parallel::stopCluster(cl)
  }
  else {
    if (method == "spearman" | method == "pearson") {
      obs <- lapply(gen.tiprates, cortest, trait_data, method)
      permu <- lapply(p.gen.tiprates, cortest, trait_data, 
                      method)
    }
    else if (method == "mann-whitney") {
      obs <- lapply(gen.tiprates, manntest, trait_data, two.tailed, 
                    trait.state)
      permu <- lapply(p.gen.tiprates, manntest, trait_data, 
                      two.tailed, trait.state)
    }
    else {
      obs <- lapply(gen.tiprates, kruskaltest, trait_data)
      permu <- lapply(p.gen.tiprates, kruskaltest, trait_data)
    }
  }
  obs <- unlist(obs)
  permu <- unlist(permu)
  obs <- obs - stat.mu
  permu <- permu - stat.mu
  if (two.tailed) {
    pval <- sum(abs(obs) <= abs(permu))/length(permu)
  }
  else {
    if (method == "spearman" | method == "pearson") {
      if (direction == "positive") {
        pval <- sum(obs <= permu)/length(permu)
      }
      else {
        pval <- sum(obs >= permu)/length(permu)
      }
    }
    else {
      pval <- sum(obs <= permu)/length(permu)
    }
  }
  if (method == "spearman" | method == "pearson") {
    obj <- list(estimate = mean(as.numeric(obs)), p.value = pval, 
                method = method, two.tailed = two.tailed)
  }
  else {
    if (ratetype == "speciation") {
      # ave.tiprate <- getTipRates(BAMM_data)$lambda.avg
      tiprate <- do.call(rbind.data.frame, BAMM_data$tipLambda)
      ave.tiprate <- apply(X = ave.tiprate, MARGIN = 1, FUN = mean, na.rm = F)
    }
    else if (ratetype == "extinction") {
      # ave.tiprate <- getTipRates(BAMM_data)$mu.avg
      tiprate <- do.call(rbind.data.frame, BAMM_data$tipMu)
      ave.tiprate <- apply(X = ave.tiprate, MARGIN = 1, FUN = mean, na.rm = F)
    }
    else if (ratetype == "net diversification") {
      # ave.tiprate <- getTipRates(BAMM_data)$lambda.avg - getTipRates(BAMM_data)$mu.avg
      tiprate_lambda <- do.call(rbind.data.frame, BAMM_data$tipLambda)
      ave.tiprate_lambda <- apply(X = tiprate_lambda, MARGIN = 1, FUN = mean, na.rm = F)
      tiprate_mu <- do.call(rbind.data.frame, BAMM_data$tipMu)
      ave.tiprate_mu <- apply(X = tiprate_mu, MARGIN = 1, FUN = mean, na.rm = F)
      ave.tiprate <- ave.tiprate_lambda - ave.tiprate_mu
    }
    else if (ratetype == "trait") {
      ave.tiprate <- getTipRates(BAMM_data)$beta.avg
    }
    l <- lapply(unique(trait_data[!is.na(trait_data)]), function(x) {
      median(ave.tiprate[which(trait_data == x)], na.rm = TRUE)
    })
    names(l) <- as.character(unique(trait_data[!is.na(trait_data)]))
    obj <- list(estimate = l, p.value = pval, method = method, 
                two.tailed = two.tailed)
  }
  obj$rate <- ratetype
  if (return.full) {
    obj$obs.corr <- as.numeric(obs)
    obj$gen <- gen
    obj$null <- as.numeric(permu)
  }
  return(obj)
}
