library(tidyverse)
library(RcppClock)
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
n <- 10000

parNodes <- rep(c(-.5,.5), p/2)
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
    MAXT = n*2,
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
fit_sgd$theta
fit_sgd2$theta
#fit_sgd$fit$sublik_pool
# path <- fit_sgd$fit$path_av_theta[(nrow(fit_sgd$fit$path_av_theta)-1000):nrow(fit_sgd$fit$path_av_theta),]
# nll_vec <- map_dbl(
#     as.list(data.frame(t(path))),
#     ~ ncl(as.matrix(data), .x, Q, VERBOSEFLAG = F)$nll
# )
# mean(nll_vec)



cpp_ctrl <- list(
    MAXT = 1000,
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
fit_scsd$theta
fit_scsd2 <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'CSGD_bernoulli',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
fit_scsd2$theta

fit_hyper2 <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'CSGD_hyper',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
fit_hyper$theta
fit_hyper2$theta
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

#### simulation test ####
theta_init <- rep(0, length(true_theta))
scls <- diag(sampleH(true_theta, DATA = as.matrix(data), CONSTRAINTS = Q, INVERTFLAG = T))
# scls <- 1/diag(sampleH(true_theta, DATA = as.matrix(data), CONSTRAINTS = Q, INVERTFLAG = F))
# scls <- rep(1,d)

fit_uc <- fit_isingGraph(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)

library(tidyverse)
sim_settings <- expand_grid(
    mod = c('OSGD', 'CSGD_bernoulli'),
    stepsize = c(.01, 1, 10, 20),
    stoc_seed = 1:5,
    maxiter = 4000,
    burn = 500
)

custom_est_fun <- function(MOD, STEPSIZE, SEED, MAXT, BURN){
    ctrl <- list(
        MAXT = MAXT,
        BURN = BURN,
        STEPSIZE = STEPSIZE,
        PAR1 = 1,
        PAR2 = STEPSIZE,
        PAR3 = 3/4,#.501,
        SCALEVEC = scls,
        NU = 1,
        SEED = SEED,
        STEPSIZEFLAG = 1
    )
    mod_obj <- fit_isingGraph(
        DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
        METHOD = MOD,
        CPP_CONTROL = ctrl,
        INIT = theta_init,
        ITERATIONS_SUBSET = seq(0, 4000, 100),
        VERBOSEFLAG = 0
    )

    return(mod_obj)
}


est_tab <- sim_settings %>%
    mutate(
        mod_obj = purrr::pmap(
            list(mod, stepsize, stoc_seed, maxiter, burn),
            function(mod_, stepsize_, stoc_seed_, maxiter_, burn_){
                custom_est_fun(
                    MOD = mod_,
                    STEPSIZE = stepsize_,
                    SEED = stoc_seed_,
                    MAXT = maxiter_,
                    BURN = burn_)
            }
        )
    )

test_mod <- est_tab %>% pluck('mod_obj', 1)
test_mod$fit$path_av_theta %>% dim()
test_mod$fit$path_grad %>% dim()
#get_tidy_path(test_mod, 'path_grad')
get_tidy_path(test_mod, 'path_av_theta')

metrics_tab <- est_tab %>%
    mutate(
        path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta')),
        #path_nll = map(mod_obj, ~get_tidy_path(.x, 'path_nll')),
        #path_grad = map(mod_obj, ~get_tidy_path(.x, 'path_grad'))
    ) %>%
    select(-mod_obj) %>%
    unnest('path_av_theta') %>%
    mutate(
        mse = map_dbl(path_av_theta, ~mean((.x-true_theta)^2)),
        #grad_norm = map_dbl(path_grad, ~norm(as.matrix(.x)))
    ) %>%
    gather(key = 'performance', value = 'val', mse)


gg <- metrics_tab %>%
    ggplot( aes(x = iter, y = val, col = factor(stepsize), group = interaction(mod, stepsize, stoc_seed)))+
    geom_line(aes(linetype = mod))+
    #geom_hline(yintercept = log(mean((Opt_u$theta-repar_theta)^2)), linetype = 'dashed')+
    facet_wrap(vars(performance), scales = 'free')+
    theme_bw()+
    scale_color_viridis_d()
plotly::ggplotly(gg, dynamicTicks = T)

true_tib <- tibble(par = 1:length(true_theta), true_val = true_theta)
num_tib <- tibble(par = 1:length(true_theta), num_val = fit_uc$theta)

av_par_tab <- est_tab %>%
    mutate(
        path_av_theta = map(mod_obj, ~get_tidy_path(.x, 'path_av_theta'))
    ) %>%
    unnest(c(path_av_theta)) %>%
    mutate(
        path_av_theta = lapply(path_av_theta, function(x){
            tib <- tibble(
                par = 1:length(x),
                val = x
            )
            tib
        })
    ) %>%
    unnest(c(path_av_theta)) %>%
    select(mod, stepsize, iter, par, val) %>%
    group_by(mod, stepsize, iter, par) %>%
    summarise(av_val = mean(val)) %>%
    mutate(
        par = as.factor(par)
    )

gg1 <- av_par_tab  %>%
    filter(iter %in% seq(0, 5000, 100)) %>%
    #mutate(av_val = map2_dbl(av_val, par_type, ~if_else(.y == 'correlation', rofz_cpp(.x), .x))) %>%
    ggplot(aes(x = iter, y = av_val))+
    geom_line(aes(linetype = mod,  col = factor(stepsize), group = interaction(mod, stepsize, par))) +
    geom_point(data = num_tib
               , aes(x = 4020, y = num_val, group = par), col = 'red', shape = 4, size = 2)+
    geom_point(data = true_tib, aes(x = 4040, y = true_val, group = par), col = 'blue', shape = 4, size = 2)+
    facet_wrap(vars(stepsize), scales = 'free') +
    theme_bw()+
    scale_color_viridis_d()
plotly::ggplotly(gg1, dynamicTicks = T)

gg1



#####
index_to_component(P=10, N=5, INDEX=49)
hyper_sampling(K = 10, N = 5, SEED = 13)
unit_sampling(N = 5, SEED = 5)
bernoulli_sampling(K = 10, N = 5, PROB = 1/5)
components_given_unit(0, 10)

theta_init <- rep(0, length(true_theta))
subtraj <- seq(0, n*5, 1000)
fit_uc <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'ucminf',
    CPP_CONTROL = list(),
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
mean((fit_uc$theta-true_theta)^2)

cpp_ctrl <- list(
    MAXT = n*5,
    BURN = 1,
    STEPSIZE = 1,
    NU = 1,
    SEED = 1,
    STEPSIZEFLAG = 0,
    SAMPLING_WINDOW = 100#n*.25
)
fit_sgd <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'standard',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = subtraj,
    VERBOSEFLAG = 0
)
fit_sgd2 <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'recycle_standard',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = subtraj,
    VERBOSEFLAG = 0
)
cpp_ctrl$EACH <- 1000

fit_sgd3 <- fit_isingGraph3(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'recycle_standard',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,

    VERBOSEFLAG = 0
)
fit_sgd3$clock
fit_sgd3$theta

fit_sgd3$fit$path_theta
fit_sgd3$fit$path_av_theta
as_tibble(fit_sgd3$fit$path_grad)

fit_sgd3$fit$iter_idx %>% length()
get_tidy_path(fit_sgd3, 'path_av_theta')
get_tidy_path3(fit_sgd3, 'path_av_theta')
get_tidy_path3(fit_sgd3, 'path_grad')

fit_sgd3$control$STEPSIZE
fit_sgd$theta
fit_sgd2$theta
fit_sgd3$theta
mean((fit_sgd$theta-true_theta)^2)
mean((fit_sgd2$theta-true_theta)^2)
mean((fit_sgd3$theta-true_theta)^2)

fit_sgd$clock
fit_sgd2$clock



cpp_ctrl <- list(
    MAXT = n*5,
    BURN = 1,
    STEPSIZE = 1,
    NU = 1,
    SEED = 1,
    STEPSIZEFLAG = 0,
    SAMPLING_WINDOW = 100#n*.25
)
fit_hyper <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'hyper',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = subtraj,
    VERBOSEFLAG = 0
)
fit_hyper2 <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'recycle_hyper',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = subtraj,
    VERBOSEFLAG = 0
)
fit_hyper$theta
fit_hyper2$theta
mean((fit_hyper$theta-true_theta)^2)
mean((fit_hyper2$theta-true_theta)^2)
fit_hyper$clock
fit_hyper2$clock

fit_be <- fit_isingGraph2(
    DATA_LIST = list(DATA = as.matrix(data), CONSTRAINTS = Q),
    METHOD = 'bernoulli',
    CPP_CONTROL = cpp_ctrl,
    #UCMINF_CONTROL = list(),
    INIT = theta_init,
    ITERATIONS_SUBSET = NULL,
    VERBOSEFLAG = 0
)
mean((fit_be$theta-true_theta)^2)
fit_be$clock
