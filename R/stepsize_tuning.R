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
    if(!(METHOD %in% c('ucminf','GD', 'OSGD', 'CSGD_bernoulli', 'CSGD_hyper'))) stop('Method not available.')
    out$method <- METHOD

    # Check stochastic control parameters
    cpp_ctrl <- check_SCSD_args(CPP_CONTROL, N = n)



    # Guarantee reproducibility stochastic optimisation
    set.seed(cpp_ctrl$SEED)

    # Collect and rearrange arguments to pass to cpp function
    args <- append(list( 'THETA_INIT' = out$theta_init), c( DATA_LIST, cpp_ctrl) )

    #args$SEED <- NULL

    sub_data <- DATA_LIST$DATA[sample(1:nrow(DATA_LIST$DATA), NCL_SAMPLE_SIZE),]
    if(METHOD == 'OSGD'){args$METHODFLAG <- 1} else if(METHOD == 'CSGD_bernoulli'){args$METHODFLAG <- 2}else if(METHOD == 'CSGD_hyper'){args$METHODFLAG <- 3}
    Rwr_ncl <- function(par){ ncl(sub_data, par, DATA_LIST$CONSTRAINTS)$nll }

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
