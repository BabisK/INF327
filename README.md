# Timeseries and Forecasting Methods
>Charalampos Kaidos

## Loading and exploring the data

The first thing to do is to load the data in the R workspace.
```R
header <- c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN", "MA", "RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf")
dates = seq(from = as.Date("1990-04-01", format='%Y-%m-%d'), to = as.Date("2005-12-01", format='%Y-%m-%d'), by = 'month')
as_data <- read.table("data_assignment.txt")
colnames(as_data) <- header
rownames(as_data) <- dates
```

Now we are going to plot all the dependent variables:
```R
par(mfrow=c(4,2))
for (col in c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN")){
  plot(ts(data = as_data[[paste(col)]], start = c(1990, 4), end = c(2004, 12), frequency = 12), ylab = paste(col))
}
par(mfrow=c(1,1))
```

![Dependent Variables](images/DependentVariables.png)

We are going to use the HFRI and DS variables for the analysis.

## Timeseries models

In the analysis and model fitting phase we will use the data from April 1990 to December 2004.

```R
HFRI_data <- ts(data = as_data$HFRI, start = c(1990, 4), end = c(2004, 12), frequency = 12)
DS_data <- ts(data = as_data$DS, start = c(1990, 4), end = c(2004, 12), frequency = 12)
```

### Stationarity

Before wokring on the models we need to make sure that the timeseries we are using are stationary processes. We will use the `ar()` function to automatically fit an AR model and perform an Augmented Dickey-Fuller test for the order decided by `ar()`.

```R
ar <- ar(HFRI_data)
adfTest(HFRI_data, lags = ar$order, type = "c")
```

```
Title:
 Augmented Dickey-Fuller Test

Test Results:
  PARAMETER:
    Lag Order: 1
  STATISTIC:
    Dickey-Fuller: -7.8606
  P VALUE:
    0.01

Description:
 Mon Feb 20 16:55:34 2017 by user: ckaidos

Warning message:
In adfTest(HFRI_data, lags = ar$order, type = "c") :
  p-value smaller than printed p-value
```

As the ADF test indicates, the HFRI series is stationary (for lag 1). P-value is less than 0.01 thus the H_0 hypothesis that HFRI is non-stationary is rejected.

```R
ar <- ar(DS_data)
adfTest(DS_data, lags = ar$order, type = "c")
```

```
Title:
 Augmented Dickey-Fuller Test

Test Results:
  PARAMETER:
    Lag Order: 2
  STATISTIC:
    Dickey-Fuller: -6.2674
  P VALUE:
    0.01

Description:
 Mon Feb 20 16:57:57 2017 by user: ckaidos

Warning message:
In adfTest(DS_data, lags = ar$order, type = "c") :
  p-value smaller than printed p-value
```

Again, as the ADF test indicated the DS series is stationary too. For DS the lag order is 2.

Since both timeseries are stationary processes we can apply ARMA models on them.

### Timeseries exploration

First we are going to work with the HFRI timeseries.

Bellow are the plots of the timeseries, the autocorrelation and partial autocorrelation plots.

```R
par(mfrow = c(3,1))
plot(HFRI_data)
acf(HFRI_data)
pacf(HFRI_data)
par(mfrow = c(1,1))
```

![HFRI timeseries](images/HFRIexploration.png)

From the ACF plot we observe correlation on lag 1 while the PACF plot indicates partial autocorrelation on lag 1 and lag 15.

The Box-Pierce test bellow confirms these observations.

```R
Box.test(HFRI_data,1,type="Box-Pierce")
```
```
	Box-Pierce test

data:  HFRI_data
X-squared = 12.244, df = 1, p-value = 0.0004667
```

As does the Ljung-Box test.

```R
Box.test(HFRI_data,1,type="Ljung-Box")
```
```
	Box-Ljung test

data:  HFRI_data
X-squared = 12.453, df = 1, p-value = 0.0004173
```

Next up is the DS timeseries.

```R
par(mfrow = c(3,1))
plot(DS_data)
acf(DS_data)
pacf(DS_data)
par(mfrow = c(1,1))
```

![DS timeseries](images/DSexploration.png)

From the ACF plot we observe correlation on lag 2 while the PACF plot indicates partial autocorrelation on lag 1 and lag 15.

The Box-Pierce test bellow confirms these observations.

```R
Box.test(DS_data,2,type="Box-Pierce")
```
```
	Box-Pierce test

data:  DS_data
X-squared = 51.85, df = 2, p-value = 5.506e-12
```

As does the Ljung-Box test.

```R
Box.test(DS_data,2,type="Ljung-Box")
```
```
	Box-Ljung test

data:  DS_data
X-squared = 52.762, df = 2, p-value = 3.49e-12
```

### Model fitting

Given the ACF and PACF plots above we choose to fit:

#### HFRI

A MA(1) model for HFRI:
```R
modelHFRI = arima(HFRI_data, order = c(0,0,1))
```
```
Call:
arima(x = HFRI_data, order = c(0, 0, 1))

Coefficients:
         ma1  intercept
      0.2357     0.0081
s.e.  0.0664     0.0018

sigma^2 estimated as 0.0003721:  log likelihood = 447.64,  aic = -889.28
```

#### DS

A MA(2) model for DS:
```R
modelDS = arima(DS_data, order = c(0,0,2))
```
```
Call:
arima(x = DS_data, order = c(0, 0, 2))

Coefficients:
         ma1     ma2  intercept
      0.6033  0.2295     0.0088
s.e.  0.0766  0.0745     0.0021

sigma^2 estimated as 0.0002294:  log likelihood = 490.28,  aic = -972.56
```

## Linear Regression

In this section we will try to use linear regression to fit the dependent variables based on the independent ones. 

![Expression](https://latex.codecogs.com/gif.latex?Y_t%20%3D%20a%20&plus;%20b_1X_%7B1%2Ct-1%7D%20&plus;%20b_1X_%7B2%2Ct-1%7D%20&plus;%20...%20&plus;%20b_kX_%7Bk%2C1%7D%20&plus;%20%5Cvarepsilon_t%20%5Chfil%20%5Cvarepsilon_t%20%5Csim%20N%280%2C%20%5Csigma%20%5E2%29)

To achieve this regression on previous values of the independent variables, we need to created new columns containing the shifted values

```R
shifted_data <- shift.column(data = as_data, columns = c("RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf"))
```

### HFRI linear model fit

We will work with the HFRI series first. We will fit the model on all the independent variables and then use a stepwise algorithm based on AIC to find the model with the largest (absolute) AIC value.
```R
fitHFRI <- lm(HFRI ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted, shifted_data)
stepHFRI <- stepAIC(fitHFRI, direction="both", trace = 0)
```
The summary of the fitted model:
```
Call:
lm(formula = HFRI ~ RUS_1_Rf_1.Shifted + MEM_Rf.Shifted + DEFSPR.Shifted + 
    VIX.Shifted, data = shifted_data)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.030908 -0.006944  0.000595  0.006436  0.055811 

Coefficients:
                     Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.0066077  0.0009076   7.280 9.53e-12 ***
RUS_1_Rf_1.Shifted  0.3248679  0.0237401  13.684  < 2e-16 ***
MEM_Rf.Shifted      0.0300763  0.0165486   1.817  0.07078 .  
DEFSPR.Shifted     -2.0108122  0.8067856  -2.492  0.01358 *  
VIX.Shifted         0.0867212  0.0307636   2.819  0.00535 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01237 on 183 degrees of freedom
Multiple R-squared:  0.6134,	Adjusted R-squared:  0.605 
F-statistic: 72.59 on 4 and 183 DF,  p-value: < 2.2e-16
```
The selected model has 4 independent variable parameters. The model was chosen based on AIC but the P-value is also significant for all but one variable.

The process and changes on the AIC value are presented with ANOVA:
```
Stepwise Model Path 
Analysis of Deviance Table

Initial Model:
HFRI ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + 
    MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + 
    SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + 
    FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted

Final Model:
HFRI ~ RUS_1_Rf_1.Shifted + MEM_Rf.Shifted + DEFSPR.Shifted + 
    VIX.Shifted


                 Step Df     Deviance Resid. Df Resid. Dev       AIC
1                                           172 0.02746477 -1628.283
2   - SBGC_Rf.Shifted  1 4.176713e-07       173 0.02746518 -1630.280
3       - MOM.Shifted  1 5.586272e-07       174 0.02746574 -1632.276
4       - HML.Shifted  1 2.821390e-06       175 0.02746856 -1634.257
5  - GSCI__Rf.Shifted  1 9.506713e-06       176 0.02747807 -1636.192
6   - SBWG_Rf.Shifted  1 1.349719e-05       177 0.02749157 -1638.100
7        - Rf.Shifted  1 1.996187e-05       178 0.02751153 -1639.963
8    - RUS_Rf.Shifted  1 3.988954e-05       179 0.02755142 -1641.691
9   - MXUS_Rf.Shifted  1 2.502533e-05       180 0.02757644 -1643.520
10      - SMB.Shifted  1 6.957733e-05       181 0.02764602 -1645.046
11  - FRBI_Rf.Shifted  1 1.551728e-04       182 0.02780119 -1645.994
12   - LHY_Rf.Shifted  1 1.871318e-04       183 0.02798833 -1646.733
```

### DS linear model fit

Following the same pattern, we fit a linear model for the DS series using all independent variables and then use the stepwise algorithm to select the optimal model.

```R
fitDS <- lm(DS ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted, shifted_data)
stepDS <- stepAIC(fitDS, direction="both", trace = 0)
```

The summary of the selected model:
```
Call:
lm(formula = DS ~ RUS_1_Rf_1.Shifted + MEM_Rf.Shifted + HML.Shifted + 
    MOM.Shifted + DEFSPR.Shifted + VIX.Shifted + Rf.Shifted, 
    data = shifted_data)

Residuals:
      Min        1Q    Median        3Q       Max 
-0.047743 -0.007798 -0.000913  0.008027  0.049501 

Coefficients:
                    Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.011292   0.002747   4.111 5.99e-05 ***
RUS_1_Rf_1.Shifted  0.129167   0.028231   4.575 8.83e-06 ***
MEM_Rf.Shifted      0.042008   0.020214   2.078  0.03912 *  
HML.Shifted         0.091155   0.032111   2.839  0.00505 ** 
MOM.Shifted         0.056104   0.021174   2.650  0.00877 ** 
DEFSPR.Shifted     -4.634487   0.949650  -4.880 2.32e-06 ***
VIX.Shifted         0.107258   0.035786   2.997  0.00311 ** 
Rf.Shifted         -1.173343   0.742434  -1.580  0.11577    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.01423 on 180 degrees of freedom
Multiple R-squared:  0.3696,	Adjusted R-squared:  0.3451 
F-statistic: 15.08 on 7 and 180 DF,  p-value: 1.951e-15
```

The stepwise algorithm has created a much more complex model with seven independent variables. The selection process follows:
```
Stepwise Model Path 
Analysis of Deviance Table

Initial Model:
DS ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + 
    MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + 
    SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + 
    FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted

Final Model:
DS ~ RUS_1_Rf_1.Shifted + MEM_Rf.Shifted + HML.Shifted + MOM.Shifted + 
    DEFSPR.Shifted + VIX.Shifted + Rf.Shifted


                Step Df     Deviance Resid. Df Resid. Dev       AIC
1                                          172 0.03589531 -1577.955
2  - SBWG_Rf.Shifted  1 2.165473e-06       173 0.03589747 -1579.944
3 - GSCI__Rf.Shifted  1 7.334253e-06       174 0.03590481 -1581.905
4      - SMB.Shifted  1 2.228979e-05       175 0.03592710 -1583.789
5   - LHY_Rf.Shifted  1 6.192050e-05       176 0.03598902 -1585.465
6  - SBGC_Rf.Shifted  1 4.392285e-05       177 0.03603294 -1587.236
7  - FRBI_Rf.Shifted  1 1.049138e-04       178 0.03613785 -1588.689
8  - MXUS_Rf.Shifted  1 1.738925e-04       179 0.03631175 -1589.787
9   - RUS_Rf.Shifted  1 1.317118e-04       180 0.03644346 -1591.106
```
