# LLfreq
## Overviwe
The LLfreq package provides Laplace approximation Bayesian Gaussian Process regressionb for one dimentional data. Althoug this package was created to correctly estimate the age frequency spectrum of detrital zircon, it can be widely applied to various data instead of histogram.

## Installation

```r
install.packages("mvtnorm")
install.packages("devtools")
library(devtools)
install_github("Tan-Furukawa/MCMCfreq")
```

## Examples

* The function `lalgp()` is effective for all one dimensional data.
```r
library(LLfreq)
d <- Osayama
e <- auto_cmpt_freq(d)
freq_graph(e, hist = T, ylab = "Frequency")
```

## Author
Tan Furukawa (古川旦)

e-mail: rpackagetan@gmail.com

## Reference

* Furukawa, T. (2019). Bayesian statistical evaluation method for detrital zircon geochronology. In Abstract of annual meeting of JpGU, Chiba, Japan, 2019.

* Riihimäki, J., and Vehtari, A. (2014). Laplace approximation for logistic Gaussian process density estimation and regression. Bayesian analysis, 9(2), 425-448.

* Tsujimori, T. (1998). Geology of the Osayama serpentinite melange in the central Chugoku Mountains,
southwestern Japan: 320 Ma blueschist-bearing serpentinite melange beneath the Oyama ophiolite.
Jour. Geol. Soc. Japan, 104, 213-231.


