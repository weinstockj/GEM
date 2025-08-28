# GEM Package

The GEM package provides tools for simulating data and fitting a GEM model using Torch in R. It includes functions to generate synthetic datasets and fit a model based on those data.

## Installation

Install with:

```r
renv::install("weinstockj/GEM")
```

## Usage

#### Simulate Data
Use `simulate_data` to generate a synthetic dataset:

This function returns a list containing:

 - x: a numeric matrix of predictors with one row for each candidate somatic mutation and one column for each feature
 - y: a numeric response vector containing zcore'd age at blood draw
 - z: an indicator for whether a candidate mutation is a "real"
 - sample_map: a list mapping samples to indices
 - true_count and naive_count for model assessment

#### Fit a GEM Model
Fit the GEM model on your data using the fit_model function:

```r
model <- fit_model(x, y, sample_map, lr = 0.1, iters = 180, num_threads = 4L)
```

The returned list includes:

 - Model coefficients (e.g. beta),
 - Fitted mutation burden estimates (fitted),
 - Loss history, and other parameters.

### Vignettes

Please see the "articles" portion of the package website for a brief vignette. 

### Testing
To run the tests for GEM, execute:

```r
devtools::test()
```

## Citation

Please cite our [manuscript](https://www.medrxiv.org/content/10.1101/2024.08.22.24312319v2) if you use this package in your research. 
