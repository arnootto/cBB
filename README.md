# The cBB package
An R package to implement the methodology described in *Modeling Bounded Count Environmental Data Using a Contaminated
Beta-Binomial Regression Model* (2024).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/cBB-RM")
```
## Example
Code to reproduce Example 5.1: Number of visits to a doctor in a year
```{r}
library(cBB)
data("MuleDeer")

b1 <- ml.mbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ 1, data=as.data.frame(df), reltol = 1e-10, hessian=T)
b3 <- ml.opt.cmbb(formula = cbind(Winter_malnutrition_n, Radiocollared_fawns - Winter_malnutrition_n) ~ 1, data=as.data.frame(df), reltol = 1e-10, method = "BFGS")
```
