simulate_data = function(n = 1e4, p = 10, k = 100) {

    sample_labels = sample(1:k, size = n)


}


torch_model = function(x, y, sample_map, lr = 0.1, iters = 180, num_threads = 20L, size = 2000) {

    tictoc::tic()

    torch_set_num_threads(num_threads)

    stopifnot(is.matrix(x))
    stopifnot(is.vector(y))
    stopifnot(ncol(x) == length(torch_predictors()))
    stopifnot(length(sample_map) == length(y))

    ids = sample(1:length(y), size = size, replace = FALSE)
    all_nwds = names(sample_map)
    nwds = all_nwds[ids]
    sample_map = purrr::map(ids, ~sample_map[[.x]]) %>%
                    purrr::set_names(nwds)
        
    logger::log_info("now selecting {length(unlist(sample_map))} mutations")

    # n x p
    x = torch_tensor(x, dtype = torch_float64())
    # k 
    y = torch_tensor((y[ids] - 60) / 11, dtype = torch_float64())

    logger::log_info("done converting to torch_tensors")

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
    # optimizer <- optim_lbfgs(c(beta, mu, log_sigma, intercept), lr, max_iter = 5)


    for(j in 1:iters) {

        logger::log_info(glue::glue("Iteration: {j}"))

        optimizer$zero_grad()
        # active = torch_sigmoid(torch_matmul(x, beta))

        # pred$zero_()
        likelihood = 0
        for(i in 1:length(ids)) {
            # pred[i] = exp(mu) * torch_sum(active[sample_map[[i]]]) 
            # pred[i] = torch_exp(mu) * torch_sum(
            #     torch_sigmoid(
            #         torch_matmul(x[sample_map[[i]], ], beta)
            #     )
            # ) 
            # expected = torch_exp(mu) * torch_sum(torch_sigmoid(
            #                                         torch_matmul(x[sample_map[[i]], ], beta)
            #                                     )
            #                             ) 
            # expected = torch_exp(mu) + torch_sum(torch_matmul(x[sample_map[[i]], ], beta)) 
            expected = intercept + torch_exp(mu) * torch_log2(torch_sum(
                torch_sigmoid(
                    inner_intercept + torch_matmul(x[sample_map[[i]], ], beta)
                    # torch_matmul(x[sample_map[[i]], ], beta)
                ))
            ) 
            pred[i] = expected

            if(i == 1) {

                logger::log_info(glue::glue("prediction = {as.numeric(expected)}, actual = {as.numeric(y[i])}"))
            }
            likelihood = likelihood + distr_normal(
                # torch_exp(mu) + torch_sum(torch_matmul(x[sample_map[[i]], ], beta)), 
                # torch_exp(mu) * torch_sum(torch_sigmoid(torch_matmul(x[sample_map[[i]], ], beta))), 
                expected, 
                torch_exp(log_sigma))$log_prob(y[i])
        }
        # browser()

        # print(pred[1])

        # y_dist = distr_normal(pred, torch_exp(log_sigma))

        # likelihood = torch_sum(y_dist$log_prob(y))
        prior = torch_sum(beta_prior$log_prob(beta)) + 
            log_sigma_prior$log_prob(log_sigma) + 
            mu_prior$log_prob(mu)

        loss = -1.0 * (likelihood + prior)
        # loss = -1.0 * (likelihood / length(y) + prior)
        # loss = -1.0 * (likelihood / length(y))

        loss_tracker = torch_cat(list(loss_tracker, loss))
        # logger::log_info(glue::glue("loss = {loss}, likelihood = {likelihood}, prior = {prior}"))
        if(j %% 3 == 0) {
            logger::log_info("loss, likelihood, prior")
            print(loss)
            print(likelihood)
            print(prior)

            print("beta[1]")
            print(beta[2])
            print("mu")
            print(torch_exp(mu))
        }

        # loss$backward(retain_graph = TRUE)
        loss$backward()
        # logger::log_info("backward()")
        optimizer$step()
        # logger::log_info("step()")
    }
    # n x p
    torch_set_num_threads(1L)
    tictoc::toc()

    return(
        list(
            "beta" = as.numeric(beta) %>% setNames(torch_predictors()),
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
