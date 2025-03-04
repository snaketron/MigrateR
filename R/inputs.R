

process_input <- function(x) {

  check_x <- function(x) {
    if(missing(x)) {
      stop("missing x")
    }
    if(is.data.frame(x)==FALSE) {
      stop("x must be a data.frame")
    }
    if(nrow(x)<=1) {
      stop("not enough rows in x")
    }
  }

  check_cols <- function(x) {
    if(!"compound" %in% colnames(x)) {
      stop("x does not have a column compound")
    }
    if(!"dose" %in% colnames(x)) {
      stop("x does not have a column dose")
    }
    if(!"sample" %in% colnames(x)) {
      stop("x does not have a column sample")
    }
    if(!"plate" %in% colnames(x)) {
      stop("x does not have a column plate")
    }
    if(!"v" %in% colnames(x)) {
      stop("x does not have a column v")
    }

    if(is.character(x[,"compound"])==FALSE) {
      stop("column compound must be character")
    }
    if(is.character(x[,"dose"])==FALSE) {
      stop("column dose must be character")
    }
    if(is.character(x[,"sample"])==FALSE) {
      stop("column sample must be character")
    }
    if(is.character(x[,"plate"])==FALSE) {
      stop("column plate must be character")
    }
    if(is.numeric(x[,"v"])==FALSE) {
      stop("column v must be numeric")
    }

    if(any(is.na(x[,"compound"]))) {
      stop("column compound contains NAs")
    }
    if(any(is.na(x[,"dose"]))) {
      stop("column dose contains NAs")
    }
    if(any(is.na(x[,"plate"]))) {
      stop("column plate contains NAs")
    }
    if(any(is.na(x[,"sample"]))) {
      stop("column sample contains NAs")
    }
    if(any(is.na(x[,"v"]))) {
      stop("column v contains NAs")
    }
  }

  check_x(x=x)
  check_cols(x=x)

  x$group <- paste0(x$compound, '|', x$dose)
  x$group_id <- as.numeric(as.factor(x$group))
  x$plate_id <- as.numeric(as.factor(x$plate))
  # x$sv <- x$v/max(x$v)
  x$sv <- scale(x = x$v, center = TRUE, scale = TRUE)
  x$sample <- paste0(x$plate, '|', x$sample, '|', x$group)
  x$well_id <- as.numeric(as.factor(x$sample))
  x$plate_group <- paste0(x$plate, '|', x$group)
  x$plate_group_id <- as.numeric(as.factor(x$plate_group_id))
  x$well_id <- as.numeric(as.factor(x$sample))
  x$N <- nrow(x)
  x$N_well <- max(x$well_id)
  x$N_plate <- max(x$plate_id)
  x$N_plate_group <- max(x$plate_group_id)
  x$N_group <- max(x$group_id)
  
  
  # transform data
  # well_id -> 
  q <- x[, c("well_id", "plate_id", "plate_group_id")]
  q <- q[duplicated(q)==F, ]
  q <- q[order(q$well_id, decreasing = F),]
  
  # group_id -> 
  z <- x[, c("plate_group_id", "group_id")]
  z <- z[duplicated(z)==F, ]
  z <- z[order(z$plate_group_id, decreasing = F),]
  
  
  return(list(d = x, x = x, map_w = q, map_pg = z))
}

process_control <- function(control_in) {
  control <- list(mcmc_warmup = 500,
                  mcmc_steps = 1500,
                  mcmc_chains = 4,
                  mcmc_cores = 1,
                  mcmc_algorithm = "NUTS",
                  adapt_delta = 0.95,
                  max_treedepth = 12)

  # if missing control_in -> use default values
  if(missing(control_in) || is.null(control_in)) {
    return(control)
  }
  if(is.list(control_in) == FALSE) {
    stop("control must be a list")
  }
  if(all(names(control_in) %in% names(control)) == FALSE) {
    stop("unrecognized elements found in control")
  }

  ns <- names(control_in)
  for (i in seq_len(length(control_in))) {
    control[[ns[i]]] <- control_in[[ns[i]]]
  }

  check_mcmc_steps(mcmc_steps = control$mcmc_steps,
                   mcmc_warmup = control$mcmc_warmup)

  check_mcmc_chains(mcmc_chains = control$mcmc_chains)

  check_mcmc_cores(mcmc_cores = control$mcmc_cores)

  return(control)
}


# MCMC Iterations check
check_mcmc_steps <- function(mcmc_steps,
                             mcmc_warmup) {

  if(length(mcmc_steps) != 1 |
     length(mcmc_warmup) != 1) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }

  if(!is.numeric(mcmc_steps) |
     !is.numeric(mcmc_warmup)) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }


  if(is.finite(x = mcmc_steps)==FALSE |
     is.finite(x = mcmc_warmup)==FALSE) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }


  if(as.integer(x = mcmc_steps) < 500 |
     as.integer(x = mcmc_warmup) < 100) {
    stop("mcmc_steps >= 500 & mcmc_warmup >= 100.")
  }


  if(as.integer(x = mcmc_steps) <= as.integer(x = mcmc_warmup)) {
    stop("mcmc_steps > mcmc_warmup")
  }
}

# MCMC Chain number check
check_mcmc_chains <- function(mcmc_chains) {
  if(length(mcmc_chains) != 1) {
    stop("mcmc_chains must be a positive integer > 0")
  }

  if(!is.numeric(mcmc_chains)) {
    stop("mcmc_chains must be a positive integer > 0")
  }

  if(is.finite(x = mcmc_chains) == FALSE) {
    stop("mcmc_chains must be a positive integer > 0")
  }

  if(as.integer(x = mcmc_chains) <= 0) {
    stop("mcmc_chains must be a positive integer > 0")
  }
}

# MCMC Cores number check
check_mcmc_cores <- function(mcmc_cores) {
  if(length(mcmc_cores) != 1) {
    stop("mcmc_cores must be a positive integer > 0")
  }

  if(is.numeric(mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }

  if(is.finite(x = mcmc_cores) == FALSE) {
    stop("mcmc_cores must be a positive integer > 0")
  }

  if(as.integer(x = mcmc_cores) <= 0) {
    stop("mcmc_cores must be a positive integer > 0")
  }
}

