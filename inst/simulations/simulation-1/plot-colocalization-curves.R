library(spatstat.geom)
library(spatstat.explore)
library(stringr)

format_spatial_variable <- function(data, spatial_variable, grid, id){
  column <- function(i){data[data[id] == i,][spatial_variable][1] } #calling the outcome function Y
  mat_transpose <- data.frame(grid = grid)
  for(i in sort(unique(data[id])[[1]])){
    mat_transpose <- cbind(mat_transpose, column(i))
  }
  mat <- t(mat_transpose[,-c(1)])
  return(mat)
}

preprocess_data <- function(data,
                            from_cell,
                            to_cell,
                            qc_cellcount_cutoff=0,
                            n_perm=50,
                            perm_yn=FALSE,
                            r_max=200,
                            inc=1,
                            image_dims,
                            summary_function='L',
                            verbose=TRUE){
  image_xmax <- image_dims[2]
  image_ymax <- image_dims[4]
  image_xmin <- image_dims[1]
  image_ymin <- image_dims[3]
  W <- spatstat.geom::owin(c(image_xmin,image_xmax),
                           c(image_ymin,image_ymax))
  
  #########################
  ##Calculate sumfun_data #
  #########################
  
  if(!(summary_function %in% c('K', 'L', 'g'))){
    stop("summary_function must be one of 'K', 'L', or 'g' ")
  }
  
  Kdata <- K_pmean_vec <- NULL
  if(summary_function=='L'){
    L_pmean_vec <- NULL
  }
  if(summary_function=='g'){
    gdata<- g_pmean_vec <- NULL
  }
  for(i in unique(data$image_number)){
    from_count <- sum(data$image_number == i & data$cell_type == from_cell)
    to_count <- sum(data$image_number == i & data$cell_type == to_cell)
    
    if(from_count > qc_cellcount_cutoff &
       to_count > qc_cellcount_cutoff){ #quality control criteria
      qc_data <- data[data$image_number == i,]
      if(perm_yn==TRUE){
        permuted_K <- NULL
        if(summary_function=='L'){
          permuted_L <- NULL
        }
        if(summary_function=='g'){
          permuted_g <- NULL
        }
        for(p in 1:n_perm){
          ppp <- spatstat.geom::ppp(x = qc_data$cell_x,
                                    y = qc_data$cell_y,
                                    marks = factor(sample(qc_data$cell_type,
                                                          size=length(qc_data$cell_type),
                                                          replace = FALSE)
                                    ),
                                    window = W)
          Kdata_temp <- data.frame(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ), image=i)
          permuted_K <- cbind(permuted_K, Kdata_temp$iso)
          if(summary_function=='L'){permuted_L <- sqrt(permuted_K/pi)}
          if(summary_function=='g'){
            gdata_temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ),
                                                           method='c', divisor ="d"), image=i)
            permuted_g <- cbind(permuted_g, gdata_temp$pcf)
          }
        }
        K_pmean <- apply(permuted_K, 1, mean);  K_pmean_vec <- c(K_pmean_vec, K_pmean)
        if(summary_function=='L'){
          L_pmean <- apply(permuted_L, 1, mean); L_pmean_vec <- c(L_pmean_vec, L_pmean)
        }
        if(summary_function=='g'){
          g_pmean <- apply(permuted_g, 1, mean); g_pmean_vec <- c(g_pmean_vec, g_pmean)
        }
        if (verbose) message(paste0('permuted outcome for image ',i,' calculated'))
      }
      ppp <- spatstat.geom::ppp(x = qc_data$cell_x, y = qc_data$cell_y,
                                marks = factor(qc_data$cell_type), window = W)
      Kdata_temp <- data.frame(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ), image=i)
      names(Kdata_temp) <- c('r','K_expect','K_obs','image_number')
      Kdata <- rbind(Kdata, Kdata_temp)
      if(summary_function=='g'){
        gdata_temp <- data.frame(spatstat.explore::pcf(spatstat.explore::Kcross(ppp, i = from_cell, j = to_cell, correction='Ripley', r=seq(0,r_max,by=inc) ),
                                                       method='c', divisor ="d"), image=i)
        names(gdata_temp) <- c('r','g_expect','g_obs','image_number')
        gdata <- rbind(gdata, gdata_temp)
      }
    }else if(from_count > qc_cellcount_cutoff &
             to_count == 0){
      r_grid <- seq(0, r_max, by = inc)
      
      Kdata_temp <- data.frame(
        r = r_grid,
        K_expect = NA_real_,
        K_obs = 0,
        image_number = i
      )
      Kdata <- rbind(Kdata, Kdata_temp)
      
      if(summary_function=='g'){
        gdata_temp <- data.frame(
          r = r_grid,
          g_expect = NA_real_,
          g_obs = 0,
          image_number = i
        )
        gdata <- rbind(gdata, gdata_temp)
      }
      
      if(perm_yn==TRUE){
        K_pmean_vec <- c(K_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        if(summary_function=='L'){
          L_pmean_vec <- c(L_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        }
        if(summary_function=='g'){
          g_pmean_vec <- c(g_pmean_vec, rep(NA_real_, nrow(Kdata_temp)))
        }
      }
    }
  }
  #Calculate the outcome variable
  if(summary_function=='L'){
    sumfun_data <- Kdata
    sumfun_data$K_expect <- sqrt(sumfun_data$K_expect/pi)
    sumfun_data$K_obs <- sqrt(sumfun_data$K_obs/pi)
    names(sumfun_data) <- c('r','L_expect','L_obs','image_number')
    if(perm_yn==TRUE){
      sumfun_data$L_pmean <- L_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$L_pmean) & sumfun_data$L_obs == 0,
                                    0,
                                    sumfun_data$L_obs - sumfun_data$L_pmean)
    }else{
      sumfun_data$L_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$L_expect) & sumfun_data$L_obs == 0,
                                    0,
                                    sumfun_data$L_obs - sumfun_data$L_expect)
    }
    names(sumfun_data) <- c('r','L_expect','L_obs','image_number','L_pmean','outcome')
  }else if(summary_function=='g'){
    sumfun_data <- gdata
    if(perm_yn==TRUE){
      sumfun_data$g_pmean <- g_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$g_pmean) & sumfun_data$g_obs == 0,
                                    0,
                                    sumfun_data$g_obs - sumfun_data$g_pmean)
    }else{
      sumfun_data$g_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$g_expect) & sumfun_data$g_obs == 0,
                                    0,
                                    sumfun_data$g_obs - sumfun_data$g_expect)
    }
    names(sumfun_data) <- c('r','g_expect','g_obs','image_number','g_pmean','outcome')
  }else{
    sumfun_data <- Kdata
    if(perm_yn==TRUE){
      sumfun_data$K_pmean <- K_pmean_vec
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$K_pmean) & sumfun_data$K_obs == 0,
                                    0,
                                    sumfun_data$K_obs - sumfun_data$K_pmean)
    }else{
      sumfun_data$K_pmean <- NA
      sumfun_data$outcome <- ifelse(is.na(sumfun_data$K_expect) & sumfun_data$K_obs == 0,
                                    0,
                                    sumfun_data$K_obs - sumfun_data$K_expect)
    }
    names(sumfun_data) <- c('r','K_expect','K_obs','image_number','K_pmean','outcome')
  }
  covariate_df <- data
  covariate_df$cell_id <- covariate_df$cell_x <- covariate_df$cell_y <- covariate_df$cell_type <- NULL
  covariate_df <- unique(covariate_df)
  sumfun_data <- merge(sumfun_data, covariate_df, by='image_number'
                       , all=FALSE) #all=FALSE:  if an image doesn't meet the qc_cutoff, don't include it
  return(sumfun_data)
}

extract_target_curve <- function(model,
                                 grid,
                                 target,
                                 form,
                                 covar_list = list(),
                                 func) {
  
  #Possible errors
  if (length(target) == 0) stop("No symbols found in target.")
  if (length(form) != length(target)) {
    stop("form must have length equal to the number of variables in target (",
         length(target), "). The length given was ", length(form), ".")
  }
  if (any(!form %in% c("coef", "term"))) {
    stop('Each entry of form must be either "coef" or "term".')
  }
  if (!is.function(func)) stop("func must be a function.")
  if (any(names(formals(func)) != target)) stop("Arguments of func must match target.")
  if (any(names(covar_list) != target)) stop("Names of covar_list must match target.")
  
  #Extract components
  sm <- stats::coef(model, n1 = length(grid))$smterms
  get_beta <- function(name) {
    name <- paste0(name, "(grid)")
    if (name %in% names(sm)) result <- sm[[name]]$coef$value
    if (length(result) > length(grid)) stop("Multiple smooth terms match for ", name)
    if (length(result) == 0)stop("No smooth term found matching ", name)
    if (name == 'Intercept(grid)'){
      result <- result + rep(summary(model)$p.coeff, length(grid)) #include the constant component of beta0
    }
    return(result)
  }
  
  #Prepare the list to perform func on
  components <- vector("list", length(target))
  names(components) <- target
  for (i in seq_along(target)) {
    v <- target[[i]]
    beta_v <- get_beta(v)
    components[[i]] <- if (form[[i]] == "coef") beta_v else beta_v * covar_list[[i]]
  }
  return(unname(do.call(func, components)))
}

get_fixed_and_re <- function(model, n_Im_total, grid, id) {
  Tm <- mgcv::predict.gam(model, type = "terms")
  nm <- colnames(Tm)
  re_col <- grep(id, nm)
  if (length(re_col) == 0) re_col <- NULL #if no RE, delete re_col and return NA for b0
  fx_cols <- grep("^(?!ti\\()", nm, perl = TRUE) #Get all columns that are NOT random effects
  if (length(fx_cols) == 0) stop("Could not find fixed effects.")
  
  b0 <- matrix(Tm[, re_col], nrow = n_Im_total, ncol = length(grid), byrow = TRUE)
  fixed <- matrix(rowSums(Tm[, fx_cols, drop = FALSE]) + as.numeric(summary(model)$p.coeff),
                  nrow = n_Im_total, ncol = length(grid), byrow = TRUE)
  list(fixed = fixed, b0 = b0)
}
#----------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)

set.seed(12345)

nPatients <- 100
nIm <- 1
counts <- seq(100, 400, by = 10)
sigma_vals <- c(20, 40, 60, 80, 100)

all_sumfun <- list()

for (ss in sigma_vals) {
  
  delta <- ss / 5
  
  ## Simulate data from Poisson point process
  window <- owin(xrange = c(0, 1000), yrange = c(0, 1000))
  g1 <- rpois(nPatients / 2, ss)
  g2 <- rpois(nPatients / 2, ss + delta)
  adjustSigma <- c(g1, g2) + 1
  
  x <- y <- cellType <- imageID <- NULL
  
  for (p in 1:nPatients) {
    for (j in 1:nIm) {
      sCount1 <- sample(counts, 1)
      sCount2 <- sample(counts, 1)
      
      a <- rpoispp(sCount1 / 1000^2, win = window)
      aDens <- density(a, sigma = adjustSigma[p], kernel = "disc")
      aDens$v <- pmax(aDens$v, 0) * sCount2 / sCount1
      b <- rpoispp(aDens)
      
      x <- c(x, a$x, b$x)
      y <- c(y, a$y, b$y)
      cellType <- c(cellType, rep("A", a$n), rep("B", b$n))
      imageID <- c(imageID, rep(paste(p, j, sep = "_"), a$n + b$n))
    }
  }
  
  imageID <- factor(imageID)
  
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  
  phenoData <- data.frame(
    imageID = unique(imageID),
    condition = ifelse(
      as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1])) <= nPatients / 2,
      "Group1",
      "Group2"
    ),
    subject = as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1]))
  )
  
  table <- cellExp
  names(table) <- c("Xcoord", "Ycoord", "cell.type", "image")
  table$patient.id <- sapply(table$image, function(x) str_split(x, "_")[[1]][1])
  table$roi <- sapply(table$image, function(x) str_split(x, "_")[[1]][2])
  names(phenoData) <- c("image", "group", "patient.id")
  table <- merge(table, phenoData, by = c("image", "patient.id"), all = TRUE)
  colnames(table) <- c("image_number", "patient_id", "cell_x", "cell_y", "cell_type", "roi", "group")
  
  ## Calculate spatial summary functions
  sumfun_data <- preprocess_data(
    table,
    from_cell = "A",
    to_cell = "B",
    qc_cellcount_cutoff = 0,
    perm_yn = FALSE,
    r_max = 200,
    inc = 1,
    image_dims = c(0, 1000, 0, 1000),
    summary_function = "L"
  )
  
  sumfun_data$group <- ifelse(sumfun_data$group == "Group1", 1, 0)
  sumfun_data$patient_id <- factor(sumfun_data$patient_id)
  sumfun_data$sigma <- ss
  sumfun_data$delta <- delta
  sumfun_data$scenario <- paste0("sigma = ", ss, ", delta = ", delta)
  
  all_sumfun[[as.character(ss)]] <- sumfun_data
}

plot_data <- do.call(rbind, all_sumfun)

ggplot(plot_data, aes(x = r, y = outcome, group = image_number, color = factor(group))) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ scenario, ncol = 2) +
  labs(
    title = "Colocalization curves by simulation setting",
    x = "r",
    y = "L(r) - r",
    color = "Group"
  ) +
  theme_bw()

#-------------------------------------------------------------------------------------------------------------
##Strauss Process

library(ggplot2)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(stringr)

set.seed(12345)

nPatients <- 50
nIm <- 3
tauB_vals <- c(2, 3, 4, 5, 6)

all_sumfun <- list()

for (tauB_base in tauB_vals) {
  
  tauB_diff <- tauB_base *.8
  
  ## Generate data from Strauss process
  x <- y <- cellType <- imageID <- NULL
  tauB_G1 <- tauB_base
  tauB_G2 <- tauB_base + tauB_diff
  
  for (p in 1:nPatients) {
    
    group <- ifelse(p <= nPatients/2, "Group1", "Group2")
    
    # subject-level tauB, shared across that subject's images
    if (group == "Group1") {
      tauB_subj <- abs(rnorm(1, mean = tauB_G1, sd = sqrt(tauB_G1 / 10)))
    } else {
      tauB_subj <- abs(rnorm(1, mean = tauB_G2, sd = sqrt(tauB_G2 / 10)))
    }
    
    for (j in 1:nIm) {
      
      Rhard  <- sample(c(0.02, 0.03, 0.04, 0.05, 0.06), 1)
      beta_p <- sample(seq(50, 100, by = 10), 1)
      muB    <- sample(seq(10, 30, by = 1), 1)
      
      A_unit <- spatstat.random::rHardcore(
        beta = beta_p,
        R    = Rhard,
        W    = owin(c(0,1), c(0,1))
      )
      
      A_x <- A_unit$x * 1000
      A_y <- A_unit$y * 1000
      nA  <- A_unit$n
      
      B_xy <- matrix(numeric(0), ncol = 2)
      
      if (nA > 0) {
        for (i in 1:nA) {
          kB <- rpois(1, muB)
          if (kB > 0) {
            ptsB <- cbind(
              rnorm(kB, A_x[i], tauB_subj),
              rnorm(kB, A_y[i], tauB_subj)
            )
            B_xy <- rbind(B_xy, ptsB)
          }
        }
      }
      
      # remove B cells outside the window
      if (nrow(B_xy) > 0) {
        inside <- B_xy[,1] >= 0 & B_xy[,1] <= 1000 &
          B_xy[,2] >= 0 & B_xy[,2] <= 1000
        B_xy <- B_xy[inside, , drop = FALSE]
      }
      
      ID <- paste(p, j, sep = "_")
      
      if (nA > 0) {
        x <- c(x, A_x)
        y <- c(y, A_y)
        cellType <- c(cellType, rep("A", nA))
        imageID  <- c(imageID, rep(ID, nA))
      }
      
      if (nrow(B_xy) > 0) {
        x <- c(x, B_xy[,1])
        y <- c(y, B_xy[,2])
        cellType <- c(cellType, rep("B", nrow(B_xy)))
        imageID  <- c(imageID, rep(ID, nrow(B_xy)))
      }
    }
  }
  
  imageID <- factor(imageID)
  
  cellExp <- data.frame(
    x = x,
    y = y,
    cellType = factor(cellType),
    imageID = imageID
  )
  
  phenoData <- data.frame(
    imageID  = unique(imageID),
    condition = ifelse(
      as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1])) <= nPatients/2,
      "Group1",
      "Group2"
    ),
    subject = as.numeric(sapply(unique(imageID), function(x) str_split(x, "_")[[1]][1]))
  )
  
  table <- cellExp
  names(table) <- c("Xcoord", "Ycoord", "cell.type", "image")
  table$patient_id <- sapply(table$image, function(x) str_split(x, "_")[[1]][1])
  table$roi        <- sapply(table$image, function(x) str_split(x, "_")[[1]][2])
  names(phenoData) <- c("image", "group", "patient_id")
  
  table <- merge(table, phenoData, by = c("image", "patient_id"), all = TRUE)
  colnames(table) <- c(
    "image_number", "patient_id", "cell_x", "cell_y",
    "cell_type", "roi", "group"
  )
  
  ## Calculate spatial summary functions
  sumfun_data <- preprocess_data(
    table,
    from_cell = "A",
    to_cell = "B",
    qc_cellcount_cutoff = 0,
    perm_yn = FALSE,
    r_max = 200,
    inc = 1,
    image_dims = c(0, 1000, 0, 1000),
    summary_function = "L"
  )
  
  sumfun_data$group <- ifelse(sumfun_data$group == "Group1", 1, 0)
  sumfun_data$patient_id <- factor(sumfun_data$patient_id)
  sumfun_data$tauB <- tauB_base
  sumfun_data$tauB_diff <- tauB_diff
  sumfun_data$scenario <- paste0("tauB = ", tauB_base, ", tauBdiff = ", tauB_diff)
  
  all_sumfun[[as.character(tauB_base)]] <- sumfun_data
}

plot_data <- do.call(rbind, all_sumfun)

ggplot(plot_data, aes(x = r, y = outcome, group = image_number, color = factor(group))) +
  geom_line(alpha = 0.25) +
  facet_wrap(~ scenario, ncol = 5) +
  labs(
    title = "Strauss process colocalization curves by simulation setting",
    x = "r",
    y = "L(r)-r",
    color = "Group"
  ) +
  theme_bw()