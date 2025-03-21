simulate_data = function(n = 1e4, p = 10, k = 100) {

    sample_labels = 1:k

    indices = sample(1:k, size = n, replace = TRUE)

    sample_map = purrr::map(sample_labels, ~which(indices == .x)) %>%
                    purrr::set_names(sample_labels)

    naive_count = purrr::map_int(sample_map, length)

    x = matrix(rnorm(n * p), n, p)

    intercept = -0.8
    beta = rnorm(p, 0, 2)
    eta = x %*% beta
    z = rbinom(n, 1, plogis(eta + intercept))

    true_count = purrr::map_dbl(sample_map, ~sum(z[.x])) %>%
        purrr::set_names(sample_labels)

    beta1 = 1.0
    y = 60 + 10 * log2(true_count) * beta1 + rnorm(k, 0, 5)

    return(list("x" = x, "y" = y, "z" = z, "sample_map" = sample_map, true_count = true_count, naive_count = naive_count))
}

wrap_simulation = function(nsims = 10, .progress = TRUE) {

    out = purrr::map_dfr(
            1:nsims, 
            function(x) {
            dat = simulate_data()
            fitted_model = fit_model(dat$x, dat$y, dat$sample_map, iters = 100, size = NULL)
            return(
                tibble::tibble(
                    cor_gem    = cor(dat$true_count, fitted_model$fitted, method = "spearman"),
                    cor_burden = cor(dat$true_count, dat$naive_count, method = "spearman")
                )

            )
        },
        .progress = .progress
    )

    return(out)
}

#' Fit a GEM model
#'
#' @param x A numeric matrix of predictors
#' @param y A numeric response vector
#' @param sample_map A list of sample indices
#' @param lr Learning rate
#' @param iters Number of iterations
#' @param num_threads Number of threads to set
#' @param size Optional subsample size (to reduce compute burden)
#' @return A list containing model parameters and fitted values
#' @export
fit_model = function(x, y, sample_map, lr = 0.1, iters = 180, num_threads = 4L, size = NULL, verbose = FALSE) {

    verbose && tictoc::tic()

    torch::torch_set_num_threads(num_threads)

    stopifnot(is.matrix(x))
    stopifnot(is.vector(y))
    stopifnot(length(sample_map) == length(y))

    if(!is.null(size)) {
        ids = sample(seq_along(y), size = size, replace = FALSE)
        all_nwds = names(sample_map)
        nwds = all_nwds[ids]
        sample_map = purrr::map(ids, ~sample_map[[.x]]) %>%
                        purrr::set_names(nwds)
    } else {
        ids = seq_along(y)
        nwds = names(sample_map)
    }

    verbose && logger::log_info("now selecting {length(unlist(sample_map))} mutations")

    # n x p
    x = torch_tensor(x, dtype = torch_float64())
    # k x 1 
    y = torch_tensor((y[ids] - 60) / 11, dtype = torch_float64())

    verbose && logger::log_info("done converting to torch_tensors")

    # p x 1
    beta = torch_zeros(ncol(x), dtype = torch_float64(), requires_grad = TRUE)
    # scalar
    mu = torch_tensor(log(0.005), requires_grad = TRUE)
    # scalar
    log_sigma = torch_tensor(log(0.1), requires_grad = TRUE)
    beta_prior = distr_normal(0, 1.0)
    log_sigma_prior = distr_normal(log(0.8), 2.0)
    mu_prior = distr_normal(log(0.005), 1.0)
    intercept = torch_tensor(-0.1, requires_grad = TRUE)
    inner_intercept = torch_tensor(-1.4, requires_grad = TRUE)
    pred = torch_zeros(length(y))
    loss_tracker = torch_tensor(1000)

    likelihood = torch_tensor(0.0, requires_grad = TRUE)
    optimizer <- optim_adam(
                    c(
                      beta,
                      mu,
                      log_sigma,
                      intercept,
                      inner_intercept
                    ),
                    lr
                )

    for(j in 1:iters) {

        verbose && logger::log_info(glue::glue("Iteration: {j}"))

        optimizer$zero_grad()

        likelihood = 0
        for(i in seq_along(ids)) {
            expected = intercept + torch_exp(mu) * torch_log2(torch_sum(
                torch_sigmoid(
                    inner_intercept + torch_matmul(x[sample_map[[i]], ], beta)
                ))
            ) 
            pred[i] = expected

            likelihood = likelihood + distr_normal(
                                        expected,
                                        torch_exp(log_sigma)
                                    )$log_prob(y[i])
        }

        prior = torch_sum(beta_prior$log_prob(beta)) +
            log_sigma_prior$log_prob(log_sigma) +
            mu_prior$log_prob(mu)

        loss = -1.0 * (likelihood + prior)

        loss_tracker = torch_cat(list(loss_tracker, loss))

        if(j %% 3 == 0 && verbose) {
            logger::log_info(
                glue::glue(
                    "loss = {round(loss$item(), 3)}, likelihood = {round(likelihood$item(), 3)}, prior = {round(prior$item(), 3)}"
                    )
                )

        }

        loss$backward()
        optimizer$step()
    }
    # n x p
    torch_set_num_threads(1L)

    verbose && tictoc::toc()

    return(
        list(
            "beta" = as.numeric(beta),
            "sigma" = as.numeric(torch_exp(log_sigma)),
            "intercept" = as.numeric(intercept),
            "inner_intercept" = as.numeric(inner_intercept),
            "x" = as.matrix(x),
            "mu" = as.numeric(mu),
            "fitted" = as.numeric(pred),
            "y" = as.numeric(y),
            "sample_map" = sample_map,
            "nwds" = nwds,
            "indicator" = as.numeric(torch_sigmoid(inner_intercept + torch_matmul(x, beta))),
            "loss" = as.numeric(loss_tracker)
        )
    )
}
