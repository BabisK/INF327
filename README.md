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
