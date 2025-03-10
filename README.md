
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IVDML

<!-- badges: start -->

[![R-CMD-check](https://github.com/cyrillsch/IVDML/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cyrillsch/IVDML/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The IVDML package implements an instrumental variable (IV) estimator for
potentially heterogeneous treatment effects in the presence of
endogeneity as presented in Scheidegger, Guo and Bühlmann (2025). The
estimator is based on double/debiased machine learning (DML)
(Chernozhukov et al., 2018) and uses efficient machine learning
instruments (MLIV) and kernel smoothing.

We assume that we observe $N$ i.i.d. copies of the model
$$Y_i = \beta(A_i)D_i + g(X_i) + \epsilon_i.$$ $Y_i$ is the response
variable and $D_i$ is the treatment variable. $X_i$ are the observed
covariates and $A_i$ is a univariate continuous covariate with respect
to which we want to consider heterogeneity (usually $A_i$ is a component
of $X_i$). $\epsilon_i$ is an error term that satisfies
$\mathbb E[\epsilon_i|X_i] = 0$, but we have endogeneity, meaning that
$\mathbb E[\epsilon_i|D_i, X_i]\neq 0$. To deal with the endogeneity, we
assume that we have access to an instrumental variable $Z_i$ that
satisfies $\mathbb E[\epsilon_i|Z_i, X_i] = 0$. To goal is to conduct
inference on the heterogeneous treatment effect $\beta(a)$ for some
specific value $a$ of $A_i$.

We quickly describe the main idea of our estimator. Assume for a moment
that the treatment effect $\beta(\cdot) = \beta$ is constant (i.e. does
not depend on $A_i$). Then, $\beta$ is identified via
$$\beta = \frac{\mathbb E[(Y_i - \mathbb E[Y_i|X_i])(\mathbb E[D_i|Z_i, X_i] - \mathbb E[D_i|X_i])]}{\mathbb E[(D_i - \mathbb E[D_i|X_i])(\mathbb E[D_i|Z_i, X_i] - \mathbb E[D_i|X_i])]}.$$
Hence, an estimate $\hat\beta$ of $\beta$ can be obtained by estimating
the so-called nuisance functions
$f(Z_i, X_i) = \mathbb E[D_i|Z_i, X_i]$,
$\phi(X_i) = \mathbb E[D_i|X_i]$ and $l(X_i) = \mathbb E[Y_i|X_i]$ using
arbitrary user-chosen machine-learning and calculating a sample version
of the identification given above. In practice, this is done using a
cross-fitting scheme (Chernozhukov et al. 2018).

If we allow for the treatment effect $\beta(\cdot)$ to depend on the
univariate continuous quantity $A_i$, and $A_i$ is a component of $X_i$,
we can make the identification
$$\beta(a) = \frac{\mathbb E[(Y_i - \mathbb E[Y_i|X_i])(\mathbb E[D_i|Z_i, X_i] - \mathbb E[D_i|X_i])|A_i = a]}{\mathbb E[(D_i - \mathbb E[D_i|X_i])(\mathbb E[D_i|Z_i, X_i] - \mathbb E[D_i|X_i])|A_i = a]}.$$
conditional on $A_i = a$. In the sample version, we then use a kernel
function $K(\cdot)$ (for example the Gaussian kernel
$K(t) = \exp(-t^2/2)/\sqrt{2 \pi}$ ) to estimate the conditional
expectations given $A_i = a$.

For a detailed discussion of the method, we refer to Scheidegger, Guo
and Bühlmann (2025). We now demonstrate, how the IVDML package is used
in practice.

## Installation

You can install the development version of IVDML from
[GitHub](https://github.com/) with

``` r
devtools::install_github("cyrillsch/IVDML")
```

## Homogeneous Treatment Effect

This is a basic example presenting the functionality of the IVDML
package for the estimation and inference for homogeneous treatment
effects. We first simulate a dataset with $N = 200$ observations.

``` r
set.seed(1)
N <- 200
Z <- rnorm(N)
X <- Z + rnorm(N)
H <- rnorm(N) # hidden confounding
D <- sin(2 * Z) - 2 * exp(-X^2) + H + 0.5 * rnorm(N)
Y <- -2 * D + tanh(X) - H + 0.5 * rnorm(N)
```

We see that the treatment effect $\beta = -2$ is homogeneous. Moreover,
there is a hidden confounding variable $H$, which is why an instrument
$Z$ is needed. An IVDML model can be fit using the following code.

``` r
library(IVDML)
fitted_ivdml <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, ml_method = "gam", S_split = 10)
```

Here, `ml_method = "gam"` indicates that the nuisance functions should
be estimated uing a generalized additive model. `S_split = 10` indicates
that the cross-fitting is repeated 10 times with 10 random sample
splits. We can obtain a coefficient estimate and its associated standard
error using the `coef()` and `se()` methods.

``` r
print(coef(fitted_ivdml, iv_method = "mlIV"))
#> [1] -2.062829
print(se(fitted_ivdml, iv_method = "mlIV"))
#> [1] 0.1256957
```

Indicating `iv_method = "mlIV"` gives the coefficient and standard error
estimate for the estimator based on the identification presented above.
Alternatively, we can also use an estimator that does not estimate the
nuisance function $f(Z_i, X_i) = \mathbb E[D_i|Z_i, X_i]$ and uses the
instruments linearly (i.e. using $Z_i - \mathbb E[Z_i|X_i]$ instead of
$\mathbb E[D_i|Z_i, X_i] - \mathbb E[D_i|X_i])]$). This estimator was
considered in Chernozhukov et al. (2018) and in Emmenegger and Bühlmann
(2021).

``` r
print(coef(fitted_ivdml, iv_method = "linearIV"))
#> [1] -2.192405
print(se(fitted_ivdml, iv_method = "linearIV"))
#> [1] 0.1991072
```

Observe that using `iv_method = "mlIV"` leads to smaller estimated
standard error than `iv_method = "linearIV"`. There are also methods for
confidence intervals for the estimated treatment effect, both standard
confidence intervals based on the point estimate and its standard error
and a robust confidence interval that is more robust with respect to
weak instrumental variables.

``` r
# mlIV
print(standard_confint(fitted_ivdml, iv_method = "mlIV", level = 0.95))
#> $CI
#>     lower     upper 
#> -2.309188 -1.816470 
#> 
#> $level
#> [1] 0.95
#> 
#> $heterogeneous_parameters
#> NULL
print(robust_confint(fitted_ivdml, iv_method = "mlIV", level = 0.95, CI_range = c(-10, 10)))
#> $CI
#>     lower     upper 
#> -2.258353 -1.761513 
#> 
#> $level
#> [1] 0.95
#> 
#> $message
#> [1] "The interval is contained in CI_range."
#> 
#> $heterogeneous_parameters
#> NULL

# linearIV
print(standard_confint(fitted_ivdml, iv_method = "linearIV", level = 0.95))
#> $CI
#>     lower     upper 
#> -2.582648 -1.802162 
#> 
#> $level
#> [1] 0.95
#> 
#> $heterogeneous_parameters
#> NULL
print(robust_confint(fitted_ivdml, iv_method = "linearIV", level = 0.95, CI_range = c(-10, 10)))
#> $CI
#>     lower     upper 
#> -2.533527 -1.596524 
#> 
#> $level
#> [1] 0.95
#> 
#> $message
#> [1] "The interval is contained in CI_range."
#> 
#> $heterogeneous_parameters
#> NULL
```

Here, `CI_range = c(-10, 10)` specifies an a priori range for the robust
confidence interval.

## Heterogeneous Treatment Effect

We now present a basic example on how to estimate a heterogeneous
treatment effect using IVDML. We first simulate a dataset with $N = 200$
observations.

``` r
set.seed(1)
N <- 200
Z <- rnorm(N)
X <- Z + rnorm(N)
A <- X
H <- rnorm(N) # hidden confounding
D <- sin(2 * Z) - 2 * exp(-X^2) + H + 0.5 * rnorm(N)
Y <- -2 * cos(A) * D + tanh(X) - H + 0.5 * rnorm(N)
```

We see that the treatment effect $\beta(A_i) = -2 \cos(A_i)$ is
heterogeneous. We fit an IVDML model using the following code.

``` r
library(IVDML)
fitted_ivdml <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = "gam", S_split = 10)
```

The arguments of the `fit_IVDML()` functions are the same as before with
the addition that we now also specify `A = A`. This is however not
strictly necessary, as `A` can also be provided when the coefficient,
standard error and confidence intervals are estimated, provided that `A`
is a deterministic function of `X` (usually a component). The `coef()`,
`se()`, `standard_confint()` and `robust_confint()` methods now need
additional arguments `a`, `kernel_name` and `bandwidth`. The normal
reference rule bandwidth according to Silverman (1986) can be obtained
using

``` r
h_normal <- bandwidth_normal(A = A)
print(h_normal)
#> [1] 0.4319743
```

For asymptotically valid inference, one needs to choose a slightly
smaller bandwidth (undersmoothing). We use the following code for point
estimate, standard errors and confidence intervals for the heterogeneous
treatment effect $\beta(0)$ evaluated at $a = 0$.

``` r
print(coef(fitted_ivdml, iv_method = "mlIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1))
#> [1] -2.003535
print(se(fitted_ivdml, iv_method = "mlIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1))
#> [1] 0.218914
print(standard_confint(fitted_ivdml, iv_method = "mlIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1), level = 0.95)
#> $CI
#>     lower     upper 
#> -2.432598 -1.574471 
#> 
#> $level
#> [1] 0.95
#> 
#> $heterogeneous_parameters
#> $heterogeneous_parameters$a
#> [1] 0
#> 
#> $heterogeneous_parameters$kernel_name
#> [1] "gaussian"
#> 
#> $heterogeneous_parameters$bandwidth
#> [1] 0.254305
print(robust_confint(fitted_ivdml, iv_method = "mlIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1, level = 0.95, CI_range = c(-10, 10)))
#> $CI
#>     lower     upper 
#> -2.362179 -1.318388 
#> 
#> $level
#> [1] 0.95
#> 
#> $message
#> [1] "The interval is contained in CI_range."
#> 
#> $heterogeneous_parameters
#> $heterogeneous_parameters$a
#> [1] 0
#> 
#> $heterogeneous_parameters$kernel_name
#> [1] "gaussian"
#> 
#> $heterogeneous_parameters$bandwidth
#> [1] 0.254305
```

As for the homogeneous treatment effect estimator, we can also calculate
theses quantities for `iv_method = "linearIV"` instead of
`iv_method = "mlIV"`.

``` r
print(coef(fitted_ivdml, iv_method = "linearIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1))
#> [1] -1.765916
print(se(fitted_ivdml, iv_method = "linearIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1))
#> [1] 0.4572597
print(standard_confint(fitted_ivdml, iv_method = "linearIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1), level = 0.95)
#> $CI
#>      lower      upper 
#> -2.6621283 -0.8697031 
#> 
#> $level
#> [1] 0.95
#> 
#> $heterogeneous_parameters
#> $heterogeneous_parameters$a
#> [1] 0
#> 
#> $heterogeneous_parameters$kernel_name
#> [1] "gaussian"
#> 
#> $heterogeneous_parameters$bandwidth
#> [1] 0.254305
print(robust_confint(fitted_ivdml, iv_method = "linearIV", a = 0, A = A, kernel_name = "gaussian", bandwidth = h_normal/N^0.1, level = 0.95, CI_range = c(-10, 10)))
#> $CI
#>      lower      upper 
#> -2.4715430  0.5675936 
#> 
#> $level
#> [1] 0.95
#> 
#> $message
#> [1] "The interval is contained in CI_range."
#> 
#> $heterogeneous_parameters
#> $heterogeneous_parameters$a
#> [1] 0
#> 
#> $heterogeneous_parameters$kernel_name
#> [1] "gaussian"
#> 
#> $heterogeneous_parameters$bandwidth
#> [1] 0.254305
```

## More Examples

More examples can be found in Scheiegger, Guo and Bühlmann (2025) and
the associated GithHb repository
[IVDML_Application](https://github.com/cyrillsch/IVDML_Application).

## References

Victor Chernozhukov, Denis Chetverikov, Mert Demirer, Esther Duflo,
Christian Hansen, Whitney Newey, and James Robins. Double/debiased
machine learning for treatment and structural parameters. *The
Econometrics Journal*, 21(1): C1–C68, 2018.

Corinne Emmenegger and Peter Bühlmann. Regularizing double machine
learning in partially linear endogenous models. *Electronic Journal of
Statistics*, 15(2):6461–6543, 2021.

Cyrill Scheidegger, Zijian Guo and Peter Bühlmann. Inference for
heterogeneous treatment effects with efficient instruments and machine
learning. Preprint, arXiv:2503.03530, 2025.

Bernard W. Silverman. Density estimation for statistics and data
analysis. Chapman & Hall/CRC monographs on statistics and applied
probability. Chapman and Hall, London, 1986.
