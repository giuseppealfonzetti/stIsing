library(tidyverse)
##### functions #####
ising_graph_lik <- function(par){
    ncl(as.matrix(data), par, Q, VERBOSEFLAG = F)$nll
}
ising_graph_gradient <- function(par){
    ncl(as.matrix(data), par, Q, VERBOSEFLAG = F)$ngradient
}
##### setup parameters ######
p <- 10
d <- p + p*(p-1)/2
n <- 4000

parNodes <- rep(c(-1,1), p/2)
parEdgesMat <- matrix(0,p, p)
counter <- 1
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        if(abs(col-row)==1 & (min(col,row)!=p/2)) parEdgesMat[row, col] <- .5
        if(abs(col-row)==p/2) parEdgesMat[row, col] <- -.5
        # parEdgesMat[row, col] <- counter
        # counter <- counter + 1
    }
}
#parEdgesMat <- parEdgesMat + t(parEdgesMat) - diag(diag(parEdgesMat))
parEdgesMat
corrplot::corrplot(parEdgesMat)
true_theta <- parNodes
for (col in 1:(p-1)) {
    for (row in (col+1):p) {
        true_theta <- c(true_theta, parEdgesMat[row, col])
    }
}
true_theta
length(true_theta)
Q <- rep(TRUE, length(true_theta))

##### generate data ####
graph_mat <- ising_from_theta_to_emat(true_theta, p)
graph_mat <- graph_mat + t(graph_mat)
thr_vec <- true_theta[1:p]

seed <- 1
set.seed(seed)
data <- IsingSampler::IsingSampler(n, graph_mat, thr_vec, 1, method = "direct")
##### optimisation #######
theta_init <- rep(0, length(true_theta))
# theta_init <- true_theta + runif(length(true_theta), -1, 1)
#
# #numDeriv::grad(ising_graph_lik, true_theta) - ising_graph_gradient(true_theta)
# system.time(
#     opt <- ucminf::ucminf(theta_init, fn = ising_graph_lik, gr = ising_graph_gradient)
# )
# opt$par
# true_theta

# set.seed(1)
# Opt <- isingGraph(
#     DATA = as.matrix(data),
#     THETA_INIT = theta_init,
#     CONSTRAINTS = Q,
#     MAXT = 10000,
#     BURN = 500,
#     STEPSIZE = 1,
#     NU = 1,
#     METHODFLAG = 1,
#     VERBOSEFLAG = F
# )
# H <- sampleH(THETA =theta_init, DATA = as.matrix(data)/1, CONSTRAINTS = Q, INVERTFLAG = F)
# diag(H)
# Hnum <- numDeriv::jacobian(ising_graph_gradient, theta_init)/n
# diag(Hnum)
# diag(H)-diag(Hnum)
# mean((Opt$path_av_theta[nrow(Opt$path_av_theta),]-true_theta)^2)
# mean((opt$par-true_theta)^2)
# mean((rep(0, length(true_theta))-true_theta)^2)
cpp_ctrl <- list(
    MAXT = 5000,
    BURN = 1000,
    NU = 1,
    SEED = 5
)
stepsize_tuning(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'OSGD',
    CPP_CONTROL = cpp_ctrl,
    INIT = theta_init,
    STEPSIZE_GRID = seq(0.5,5,.5),
    VERBOSEFLAG = 0
)

fit_uc <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'ucminf',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)

cpp_ctrl <- list(
    MAXT = 10000,
    BURN = 500,
    STEPSIZE = 1.5,
    NU = 1,
    SEED = 1
)
fit_sgd <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'OSGD',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
#fit_sgd$fit$sublik_pool
# path <- fit_sgd$fit$path_av_theta[(nrow(fit_sgd$fit$path_av_theta)-1000):nrow(fit_sgd$fit$path_av_theta),]
# nll_vec <- map_dbl(
#     as.list(data.frame(t(path))),
#     ~ ncl(as.matrix(data), .x, Q, VERBOSEFLAG = F)$nll
# )
# mean(nll_vec)



cpp_ctrl <- list(
    MAXT = 10000,
    BURN = 500,
    STEPSIZE = 1.5,
    NU = 1,
    SEED = 1
)

fit_scsd <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'CSGD_bernoulli',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)

fit_hyper <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'CSGD_hyper',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
# fit_hyper$theta
# fit_scsd$theta
#which(fit_hyper$fit$sampling_weights ==1)
#fit$theta
mean((fit_uc$theta-true_theta)^2)
mean((fit_sgd$theta-true_theta)^2)
mean((fit_scsd$theta-true_theta)^2)
mean((fit_hyper$theta-true_theta)^2)

mean((theta_init-true_theta)^2)


gg <- get_tidy_path(fit_sgd, 'path_av_theta') %>%
    mutate( mod = 'OSGD') %>%
    bind_rows(
        get_tidy_path(fit_scsd, 'path_av_theta') %>%
            mutate( mod = 'CSGD_bernoulli')
    ) %>%
    bind_rows(
        get_tidy_path(fit_hyper, 'path_av_theta') %>%
            mutate( mod = 'CSGD_hyper')
    ) %>%
    mutate(
        mse = map_dbl(path_av_theta, ~mean((.x-true_theta)^2))
    ) %>%
    ggplot( aes(x = iter, y = (mse), col = mod))+
    geom_line()+
    geom_hline(yintercept = (mean((fit_uc$theta-true_theta)^2)), linetype = 'dashed')+
    theme_minimal()
plotly::ggplotly(gg)

fit_uc$time; fit_sgd$time; fit_scsd$time; fit_hyper$time
fit_sgd$clock
fit_scsd$clock
fit_hyper$clock
# eta <- 5
# ncl(as.matrix(data[1,]), theta_init, Q, VERBOSEFLAG = F)$nll
# gr <- ncl(as.matrix(data[1,]), theta_init, Q, VERBOSEFLAG = F)$ngr
# ncl(as.matrix(data[1,]), theta_init-eta*gr, Q, VERBOSEFLAG = F)$nll

# mat <- matrix(rep(c(1,2,3), 5), 5,3, byrow = T)
# purrr::transpose(.l = mat)
# as.list(data.frame(t(mat)))
