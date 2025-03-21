test_that("GEM torch model", {

  sim = GEM:::simulate_data(5000, 5, 40)
  model = GEM::fit_model(sim$x, sim$y, sim$sample_map, 0.1, 120, 2L, NULL)

  expect_type(model, "list")
  expect_type(model$fitted, "double")  
  expect_type(model$intercept, "double")
  expect_type(model$beta, "double")
  expect_gt(cor(sim$y, model$fitted, method = "spearman"), 0.20)
  expect_gt(cor(sim$true_count, model$fitted, method = "spearman"), 0.30)
  
})
