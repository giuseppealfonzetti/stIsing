#'@export
check_SCSD_args <- function(ARGS, N, D){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N,0)
    if(is.null(ARGS$BURN)) out$BURN <- 0
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1
    if(is.null(ARGS$NU)) out$NU <- 1
    if(is.null(ARGS$SEED)) out$SEED <- 123
    if(is.null(ARGS$SCALEVEC)) out$SCALEVEC <- rep(1, D)

    return(out)
}

#'@export
check_stoc_args <- function(ARGS, N, D){

    out <- ARGS

    if(is.null(ARGS$MAXT)) out$MAXT <- N
    if(is.null(ARGS$BURN)) out$BURN <- N*.25
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1
    if(is.null(ARGS$NU)) out$NU <- 1
    if(is.null(ARGS$SEED)) out$SEED <- 123
    if(is.null(ARGS$SCALEVEC)) out$SCALEVEC <- rep(1, D)
    if(is.null(ARGS$EACH)) out$EACH <- 1
    if(is.null(ARGS$HOLDOUTFLAG)) out$HOLDOUTFLAG <- F

    return(out)
}

#'@export
check_GD_args <- function(ARGS, N){

    out <- list('MAXT' = ARGS$MAXT, 'STEPSIZE' = ARGS$STEPSIZE)

    if(is.null(ARGS$MAXT)) out$MAXT <- round(N^.75,0)
    if(is.null(ARGS$STEPSIZE)) out$STEPSIZE <- 1e-2

    return(out)
}

#'@export
get_tidy_path <- function(MOD_OBJ, PATH_LAB){
    iters <- MOD_OBJ$iterations_subset
    path  <- MOD_OBJ$fit[[PATH_LAB]]

    if(PATH_LAB%in%c('path_nll')){
        out <- dplyr::tibble(iter = iters) %>%
            dplyr::mutate(
                path_chosen = c(path,NA)
            )
    }else if(PATH_LAB%in%c('path_grad')){
        out <- dplyr::tibble(iter = iters) %>%
            dplyr::mutate(
                path_chosen = split(t(rbind(path, rep(NA, ncol(path)))), rep(1:(nrow(path)+1), each = ncol(path)))
            )

    }else{
        out <- dplyr::tibble(iter = iters) %>%
            dplyr::mutate(
                path_chosen = split(t(path), rep(1:nrow(path), each = ncol(path)))
            )
    }


    colnames(out) <- c('iter', PATH_LAB)

    return(out)
}

#'@export
ising_from_theta_to_emat <- function(par, p){
    emat <- matrix(0, p, p)
    counter <- p+1
    for (col in 1:(p-1)) {
        for (row in (col+1):p) {
            #cat('row:', row, 'col:', col, ' counter:', counter, '\n')
            emat[row, col] <- par[counter]
            counter <- counter + 1
        }
    }
    return(emat)
}

#'@export
ising_from_graph_to_theta <- function(graph, intercepts){
    p <- length(intercepts)
    par <- intercepts

    counter <- p+1
    for (col in 1:(p-1)) {
        for (row in (col+1):p) {
            #cat('row:', row, 'col:', col, ' counter:', counter, '\n')
            par[counter] <- graph[row, col]
            counter <- counter + 1
        }
    }

    return(par)
}
