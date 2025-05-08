# The cBB package
An R package to implement the methodology described in *Modeling Bounded Count Environmental Data Using a Contaminated Beta-Binomial Regression Model* (2025).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/cBB")
```
## Example
Code to reproduce Example 5.1: Mule Deer Mortality
```{r}
library(cBB)
data("MuleDeer")

#intercept only cBB-RM
mbb <- ml.mbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ 1, data=MuleDeer, reltol = 1e-10, method="Nelder-Mead")
cmbb<- ml.cmbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ 1, data=MuleDeer, reltol = 1e-10, method = "Nelder-Mead")

#state as covariate
mbb <- ml.mbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ State, data=MuleDeer, reltol = 1e-10, method = "BFGS")
cmbb <- ml.cmbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ State, data=MuleDeer, reltol = 1e-10, method = "BFGS")
```
