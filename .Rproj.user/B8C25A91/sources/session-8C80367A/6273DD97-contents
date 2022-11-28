
<!-- README.md is generated from README.Rmd. Please edit that file -->

# havok

**NOTE: This package is in active development (beta), updates regularly,
and may have bugs. Please report any bugs to the package maintainer.**

This package allows for modeling of chaotic systems as intermittently
forced linear systems through the use of Hankel Alternative View of
Koopman (HAVOK) analysis (Brunton, Brunton, Proctor, Kaiser, & Kutz,
2017). This package has additional functionality for the SINDy algorithm
(Brunton, Proctor, & Kutz, 2016).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("RobertGM111/havok")
```

<!--
``` r
install.packages("havok")
```
-->

## Example

### Lorenz attractor

Simulate data from a Lorenz attractor.

``` r
library(havok)
library(deSolve)


#Generate Data
##Set Lorenz Parameters
parameters <- c(s = 10, r = 28, b = 8/3)
n <- 3
state <- c(X=-8, Y=8, Z=27) ##Inital Values

dt<-0.001
tspan<-seq(dt,200,dt)
N<-length(tspan)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- s * (Y - X)
    dY <- X * (r - Z) - Y
    dZ <- X * Y - b * Z
    list(c(dX, dY, dZ))
  })
}

out <- ode(y = state, times = tspan, func = Lorenz, parms = parameters, rtol = 1e-12, atol = 1e-12)
xdat <- out[,"X"]
t <- out[,"time"]

# Run HAVOK Analysis
hav <- havok(xdat = xdat, dt = dt)
```

To plot the resulting time series and forcing term use `plot(hav)`

``` r
plot(hav, what = "both")
```

<img src="man/figures/README-plotting-1.png" width="100%" />

## References

-   Brunton, S. L., Brunton, B. W., Proctor, J. L., Kaiser, E., &
    Kutz, J. N. (2017). Chaos as an intermittently forced linear system.
    Nature communications, 8(1), 19.

-   Brunton, S. L., Proctor, J. L., & Kutz, J. N. (2016). Discovering
    governing equations from data by sparse identification of nonlinear
    dynamical systems. Proceedings of the National Academy of Sciences,
    113(15), 3932-3937.
