#'@export
fit_isingGraph <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS'),
        METHOD,
        CPP_CONTROL = list(
            MAXT = 1000,
            BURN = 500,
            STEPSIZE = .01,
            STEPSIZE0 = NULL,
            NU = 1,
            SEED = 123
        ),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
        INIT = NULL,
        ITERATIONS_SUBSET = NULL,
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
    if(!(METHOD %in% c('ucminf','GD', 'OSGD', 'CSGD_bernoulli', 'CSGD_hyper'))) stop('Method not available.')
    out$method <- METHOD

    # Numerical optimisation
    if(METHOD == 'ucminf'){

        message('2. Optimising with ucminf...')
        # R wrapper of cpp function for negative composite log-likelihood
        Rwr_ncl <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$nll }
        #as.matrix(data), par, Q, VERBOSEFLAG = F
        # R wrapper of cpp function for negative composite score
        Rwr_ngr <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$ngradient }

        # list of ucminf args
        args <- list(
            'par' = out$theta_init,
            'fn' = Rwr_ncl,
            'gr' = Rwr_ngr,
            'control' = UCMINF_CONTROL$ctrl,
            'hessian' = UCMINF_CONTROL$hessian)

        # optimisation
        opt <- do.call(ucminf::ucminf, args)
        out$fit <- opt

        out$control <- UCMINF_CONTROL
        out$theta   <- opt$par

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('3. Done! (', round(out$time,2),' secs)')
        return(out)
    }

    # Stochastic approximation of numerical optimiser
    if(METHOD == 'OSGD' | METHOD == 'CSGD_bernoulli' | METHOD == 'CSGD_hyper'){

        message(paste0('2. Optimising with ', METHOD, '...'))

        # Check stochastic control parameters
        cpp_ctrl <- check_SCSD_args(CPP_CONTROL, N = n, D = d)

        # Check iterations selected
        if(!is.null(ITERATIONS_SUBSET)){
            out$iterations_subset <- c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT)
        }else{
            out$iterations_subset <- 0:cpp_ctrl$MAXT
        }

        # Guarantee reproducibility stochastic optimisation
        set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

        #args$SEED <- NULL

        if(METHOD == 'OSGD'){args$METHODFLAG <- 1} else if(METHOD == 'CSGD_bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'CSGD_hyper'){args$METHODFLAG <- 3}
        #args$METHODFLAG <- dplyr::if_else(METHOD == 'OSGD', 1, 2)


        fit <- do.call(isingGraph, args)
        fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,]
        fit$path_av_theta <- fit$path_av_theta[out$iterations_subset + 1,]
        fit$path_grad     <- fit$path_grad[out$iterations_subset,]
        #fit$path_ncl      <- fit$path_ncl[out$iterations_subset]


        fit$methodflag <- NULL



        out$control <- cpp_ctrl
        out$fit <- fit
        out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]
        if('RcppClock'%in% (.packages())) out$clock <- summary(clock, units = 's')

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

        message('\n3. Done! (', round(out$time,2),' secs)')
        return(out)

    }
}


#'@export
fit_isingGraph2 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS'),
        METHOD,
        CPP_CONTROL = list(
            MAXT = 1000,
            BURN = 500,
            STEPSIZE = .01,
            STEPSIZE0 = NULL,
            NU = 1,
            SEED = 123
        ),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
        INIT = NULL,
        ITERATIONS_SUBSET = NULL,
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
    if(!(METHOD %in% c('ucminf','GD', 'standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper'))) stop('"METHOD" not available.')
    out$method <- METHOD

    # Numerical optimisation
    if(METHOD == 'ucminf'){

        message('2. Optimising with ucminf...')
        # R wrapper of cpp function for negative composite log-likelihood
        Rwr_ncl <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$nll }
        #as.matrix(data), par, Q, VERBOSEFLAG = F
        # R wrapper of cpp function for negative composite score
        Rwr_ngr <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$ngradient }

        # list of ucminf args
        args <- list(
            'par' = out$theta_init,
            'fn' = Rwr_ncl,
            'gr' = Rwr_ngr,
            'control' = UCMINF_CONTROL$ctrl,
            'hessian' = UCMINF_CONTROL$hessian)

        # optimisation
        opt <- do.call(ucminf::ucminf, args)
        out$fit <- opt

        out$control <- UCMINF_CONTROL
        out$theta   <- opt$par

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('3. Done! (', round(out$time,2),' secs)')
        return(out)
    }

    # Stochastic approximation of numerical optimiser
    if(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper')){

        message(paste0('2. Optimising with ', METHOD, '...'))

        # Check stochastic control parameters
        cpp_ctrl <- check_SCSD_args(CPP_CONTROL, N = n, D = d)

        # Check iterations selected
        if(!is.null(ITERATIONS_SUBSET)){
            out$iterations_subset <- c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT)
        }else{
            out$iterations_subset <- 0:cpp_ctrl$MAXT
        }

        # Guarantee reproducibility stochastic optimisation
        set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

        #args$SEED <- NULL

        if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}
        #args$METHODFLAG <- dplyr::if_else(METHOD == 'OSGD', 1, 2)


        fit <- do.call(isingGraph2, args)
        fit$path_theta    <- fit$path_theta[out$iterations_subset + 1,]
        fit$path_av_theta <- fit$path_av_theta[out$iterations_subset + 1,]
        fit$path_grad     <- fit$path_grad[out$iterations_subset,]
        #fit$path_ncl      <- fit$path_ncl[out$iterations_subset]


        fit$methodflag <- NULL



        out$control <- cpp_ctrl
        out$fit <- fit
        out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]
        if('RcppClock'%in% (.packages())) out$clock <- summary(clock, units = 's')

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

        message('\n3. Done! (', round(out$time,2),' secs)')
        return(out)

    }
}

#'@export
fit_isingGraph3 <- function(
        DATA_LIST = list('DATA', 'CONSTRAINTS', 'HOLDOUT'),
        METHOD,
        CPP_CONTROL = list(
            MAXT = 1000,
            BURN = 500,
            STEPSIZE = .01,
            STEPSIZE0 = NULL,
            NU = 1,
            SEED = 123
        ),
        UCMINF_CONTROL = list('ctrl' = list(), 'hessian' = 0 ),
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

    # Numerical optimisation
    if(METHOD == 'ucminf'){

        message('2. Optimising with ucminf...')
        # R wrapper of cpp function for negative composite log-likelihood
        Rwr_ncl <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$nll }
        #as.matrix(data), par, Q, VERBOSEFLAG = F
        # R wrapper of cpp function for negative composite score
        Rwr_ngr <- function(par){ ncl(DATA_LIST$DATA, par, DATA_LIST$CONSTRAINTS)$ngradient }

        # list of ucminf args
        args <- list(
            'par' = out$theta_init,
            'fn' = Rwr_ncl,
            'gr' = Rwr_ngr,
            'control' = UCMINF_CONTROL$ctrl,
            'hessian' = UCMINF_CONTROL$hessian)

        # optimisation
        opt <- do.call(ucminf::ucminf, args)
        out$fit <- opt

        out$control <- UCMINF_CONTROL
        out$theta   <- opt$par

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])
        message('3. Done! (', round(out$time,2),' secs)')
        return(out)
    }

    # Stochastic approximation of numerical optimiser
    if(METHOD %in% c('standard', 'bernoulli', 'hyper', 'recycle_standard', 'recycle_hyper')){

        message(paste0('2. Optimising with ', METHOD, '...'))

        # Check stochastic control parameters
        cpp_ctrl <- check_SCSD_args(CPP_CONTROL, N = n, D = d)

        # # Check iterations selected
        # if(!is.null(ITERATIONS_SUBSET)){
        #     out$iterations_subset <- c(0, ITERATIONS_SUBSET[ITERATIONS_SUBSET < cpp_ctrl$MAXT], cpp_ctrl$MAXT)
        # }else{
        #     out$iterations_subset <- 0:cpp_ctrl$MAXT
        # }

        # Guarantee reproducibility stochastic optimisation with bernoulli sampling
        # For the other schemes the seed is passed directly to the cpp function
        set.seed(cpp_ctrl$SEED)

        # Collect and rearrange arguments to pass to cpp function
        args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

        #args$SEED <- NULL

        if(METHOD == 'standard'){args$METHODFLAG <- 1} else if(METHOD == 'bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'hyper'){args$METHODFLAG <- 3}else if(METHOD == 'recycle_standard'){args$METHODFLAG <- 4}else if(METHOD == 'recycle_hyper'){args$METHODFLAG <- 5}
        #args$METHODFLAG <- dplyr::if_else(METHOD == 'OSGD', 1, 2)


        fit <- do.call(isingGraph3, args)

        message(paste0('\n3. Storing results'))

        fit$path_theta    <- reduce(fit$path_theta, rbind)
        fit$path_av_theta <- reduce(fit$path_av_theta, rbind)
        fit$path_grad     <- reduce(fit$path_grad, rbind)
        #fit$path_ncl      <- fit$path_ncl[out$iterations_subset]


        fit$methodflag <- NULL



        out$control <- cpp_ctrl
        out$fit <- fit
        out$theta <- fit$path_av_theta[nrow(fit$path_av_theta),]
        if('RcppClock'%in% (.packages())) out$clock <- summary(clock, units = 's')

        end_time <- Sys.time()
        out$time <- as.numeric(difftime(end_time, start_time, units = 'secs')[1])

        message('\n4. Done! (', round(out$time,2),' secs)')
        return(out)

    }
}
