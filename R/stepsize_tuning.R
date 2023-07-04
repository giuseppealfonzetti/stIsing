#'@export
stepsize_tuning <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS'),
        METHOD,
        CPP_CONTROL = list(),
        INIT = NULL,
        STEPSIZE_GRID = NULL,
        VERBOSEFLAG = 0,
        NCL_SAMPLE_SIZE = 1000
){

    out <- list()
    start_time <- Sys.time()

    p <- ncol(DATA_LIST$DATA)
    d <- p + p*(p-1)/2

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
        out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at zero vector.')
        out$theta_init <-  rep(0, d)
    }

    # Check if method entry is correct
    if(!(METHOD %in% c('ucminf','GD', 'standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
    out$method <- METHOD

    # Check holdout
    if(is.null(DATA_LIST$HOLDOUT)){
        CPP_CONTROL$HOLDOUTFLAG <- F
        DATA_LIST$HOLDOUT <- matrix(0, 1, p)
    }else{
        if(ncol(DATA_LIST$HOLDOUT)!=ncol(DATA_LIST$DATA)) stop(' "DATA" and "HOLDOUT" must have the same number of nodes.')
        CPP_CONTROL$HOLDOUTFLAG <- T

    }

    # Check if method entry is correct
    if(!(METHOD %in% c('ucminf','GD', 'OSGD', 'CSGD_bernoulli', 'CSGD_hyper'))) stop('Method not available.')
    out$method <- METHOD

    # Check stochastic control parameters
    cpp_ctrl <- check_stoc_args(CPP_CONTROL, N = n, D = d)

    # Guarantee reproducibility stochastic optimisation
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

    #args$SEED <- NULL

    if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}
    Rwr_ncl <- function(par){ ncl(DATA_LIST$HOLDOUT, par, DATA_LIST$CONSTRAINTS)$nll }

    fun_grid <- pbapply::pbsapply(STEPSIZE_GRID, function(x){
        args$STEPSIZE <- x
        fit <- do.call(isingGraph, args)
        # fun <- mean(fit$path_ncl[CPP_CONTROL$BURN:length(fit$path_ncl)])
        #fun <- fit$path_ncl[length(fit$path_ncl)]
        fun <- Rwr_ncl(fit$path_av_theta[nrow(fit$path_av_theta),])

        return(fun)
    })

    out <- list(
        'tab' = cbind('stepsize' = STEPSIZE_GRID, 'ncl' = fun_grid),
        'chosen_eta' = STEPSIZE_GRID[which.min(fun_grid)]
    )

    return(out)
}

#'@export
stepsize_tuning3 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS', 'HOLDOUT'),
        METHOD,
        CPP_CONTROL = list(),
        STEPSIZE_GRID = NULL,
        INIT = NULL,
        VERBOSEFLAG = 0
){
    out <- list()
    start_time <- Sys.time()

    p <- ncol(DATA_LIST$DATA)
    d <- p + p*(p-1)/2

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
        out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at zero vector.')
        out$theta_init <-  rep(0, d)
    }

    # Check if method entry is correct
    if(!(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
    out$method <- METHOD

    # Check holdout
    if(is.null(DATA_LIST$HOLDOUT)){
        CPP_CONTROL$HOLDOUTFLAG <- F
        DATA_LIST$HOLDOUT <- matrix(0, 1, p)
    }else{
        if(ncol(DATA_LIST$HOLDOUT)!=ncol(DATA_LIST$DATA)) stop(' "DATA" and "HOLDOUT" must have the same number of nodes.')
        CPP_CONTROL$HOLDOUTFLAG <- T

    }

    # Check stochastic control parameters
    cpp_ctrl <- check_stoc_args(CPP_CONTROL, N = n, D = d)

    # Guarantee reproducibility stochastic optimisation with bernoulli sampling
    # For the other schemes the seed is passed directly to the cpp function
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

    if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}

    # objective function on holdout data
    Rwr_ncl <- function(par){ ncl(DATA_LIST$HOLDOUT, par, DATA_LIST$CONSTRAINTS)$nll/nrow(DATA_LIST$HOLDOUT) }

    fun_grid <- pbapply::pblapply(STEPSIZE_GRID, function(x){
        args$STEPSIZE <- x
        fit <- do.call(isingGraph3, args)

        tib <- tibble(iter = fit$iter_idx)
        tib$path <- fit$path_av_theta
        tib$stepsize <- x
        tib %>%
            mutate(nll = map_dbl(path, ~Rwr_ncl(.x)))

    })

    grid = reduce(fun_grid, rbind)
    best = grid %>% filter(nll == min(nll))

    return(list(grid = grid,
                best_step = best %>% pluck('stepsize', 1),
                best_theta = best %>% pluck('path', 1),
                best_nll = best %>% pluck('nll', 1))
    )
}

#'@export
stepsize_tuning4 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS', 'HOLDOUT'),
        METHOD,
        CPP_CONTROL = list(),
        STEPSIZE_INIT = NULL,
        LENGTH = 0.5,
        INIT = NULL,
        VERBOSEFLAG = 0,
        MAXATTEMPT = 10
){
    out <- list()
    start_time <- Sys.time()

    p <- ncol(DATA_LIST$DATA)
    d <- p + p*(p-1)/2

    # Check Initialisation
    if(is.vector(INIT)){
        if(length(INIT)!=d)
            stop(paste0('init vector has length ', length(INIT), ' instead of ', d ,'.'))
        else
            message('1. Initialising at init vector.')
        out$theta_init <-  INIT
    }else{
        if(is.null(INIT))
            message('1. Initialising at zero vector.')
        out$theta_init <-  rep(0, d)
    }

    # Check if method entry is correct
    if(!(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
    out$method <- METHOD

    # Check holdout
    if(is.null(DATA_LIST$HOLDOUT)){
        CPP_CONTROL$HOLDOUTFLAG <- F
        DATA_LIST$HOLDOUT <- matrix(0, 1, p)
    }else{
        if(ncol(DATA_LIST$HOLDOUT)!=ncol(DATA_LIST$DATA)) stop(' "DATA" and "HOLDOUT" must have the same number of nodes.')
        CPP_CONTROL$HOLDOUTFLAG <- T

    }

    # Check stochastic control parameters
    cpp_ctrl <- check_stoc_args(CPP_CONTROL, N = n, D = d)

    # Guarantee reproducibility stochastic optimisation with bernoulli sampling
    # For the other schemes the seed is passed directly to the cpp function
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

    if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}

    # objective function on holdout data
    Rwr_ncl <- function(par){ ncl(DATA_LIST$HOLDOUT, par, DATA_LIST$CONSTRAINTS)$nll/nrow(DATA_LIST$HOLDOUT) }

    cond <- T
    args$STEPSIZE <- STEPSIZE_INIT
    fit <- do.call(isingGraph3, args)
    nll <- Rwr_ncl(fit$path_av_theta[[length(fit$path_av_theta)]])
    tib <- tibble(stepsize = args$STEPSIZE, hnll = nll)
    tib$theta <- list(fit$path_av_theta[[length(fit$path_av_theta)]])
    tib$path_theta <- list(tibble(iter = fit$iter_idx, path_av_theta = fit$path_av_theta))

    res <- tib
    new_args <- args
    for (i in 1:MAXATTEMPT) {
        new_args$STEPSIZE <- new_args$STEPSIZE*LENGTH
        #cat(paste0('Using stepsize ', new_args$STEPSIZE, ':\n'))

        new_fit <- do.call(isingGraph3, new_args)
        new_nll <- Rwr_ncl(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])

        #cat(paste0('\n', round(new_nll, 4), '\n'))

        tib <- tibble(stepsize = new_args$STEPSIZE,  hnll = new_nll)
        tib$theta <- list(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])
        tib$path_theta <- list(tibble(iter = new_fit$iter_idx, path_av_theta = new_fit$path_av_theta))
        res <- res %>% bind_rows(tib)

        cond <- (new_nll < nll)
        # if(!is.numeric(new_nll)){break}else{
        #     cond <- (new_nll < nll)
        #     nll <- new_nll
        # }
        if(!is.logical(cond)){stop('likelihood differene undefined')}
        if(!cond) break;
        nll <- new_nll
        #path_theta <- tibble(iter = new_fit$iter_idx, path_av_theta = new_fit$path_av_theta)

    }
    # while (cond) {
    #     new_args$STEPSIZE <- new_args$STEPSIZE*LENGTH
    #     #cat(paste0('Using stepsize ', new_args$STEPSIZE, ':\n'))
    #
    #     new_fit <- do.call(isingGraph3, new_args)
    #     new_nll <- Rwr_ncl(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])
    #
    #     #cat(paste0('\n', round(new_nll, 4), '\n'))
    #
    #     tib <- tibble(stepsize = new_args$STEPSIZE,  hnll = new_nll)
    #     tib$theta <- list(new_fit$path_av_theta[[length(new_fit$path_av_theta)]])
    #     res <- res %>% bind_rows(tib)
    #     if(!is.numeric(new_nll)){break}else{
    #         cond <- (new_nll < nll)
    #         nll <- new_nll
    #     }
    #     if(!is.logical(cond)){stop('likelihood differene undefined')}
    #
    #
    # }


    best = res %>% filter(hnll == min(hnll))

    return(list(
        best_theta = best %>% pluck('theta', 1),
        grid = res,
        best_step = best %>% pluck('stepsize', 1),
        best_nll = best %>% pluck('hnll', 1),
        best_path = best %>% pluck('path_theta', 1)
        )
    )
}
