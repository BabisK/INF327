install.packages("fUnitRoots")
install.packages("useful")
library("fUnitRoots")
library("useful")
library("MASS")

# Load and name the data
header <- c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN", "MA", "RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf")
dates = seq(from = as.Date("1990-04-01", format='%Y-%m-%d'), to = as.Date("2005-12-01", format='%Y-%m-%d'), by = 'month')
as_data <- read.table("data_assignment.txt")
colnames(as_data) <- header
rownames(as_data) <- dates

# Plot all time series
par(mfrow=c(4,2))
for (col in c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN")){
  plot(ts(data = as_data[[paste(col)]], start = c(1990, 4), end = c(2004, 12), frequency = 12), ylab = paste(col))
}
par(mfrow=c(1,1))

# Analyze HFRI
par(mfrow = c(3,1))
HFRI_data <- ts(data = as_data$HFRI, start = c(1990, 4), end = c(2004, 12), frequency = 12)
plot(HFRI_data)
ar <- ar(HFRI_data)
adfTest(HFRI_data, lags = ar$order, type = "ct")

acf(HFRI_data)
pacf(HFRI_data)
par(mfrow = c(1,1))

boxPierce <- Box.test(HFRI_data,1,type="Box-Pierce")
boxLjung <- Box.test(HFRI_data,1,type="Ljung-Box")

modelHFRI = arima(HFRI_data, order = c(0,0,1))

# Analyze DS
par(mfrow = c(3,1))
DS_data <- ts(data = as_data$DS, start = c(1990, 4), end = c(2004, 12), frequency = 12)
plot(DS_data)
ar <- ar(DS_data)
adfTest(DS_data, lags = ar$order, type = "ct")

acf(DS_data)
pacf(DS_data)
par(mfrow = c(1,1))

boxPierce <- Box.test(DS_data,2,type="Box-Pierce")
boxLjung <- Box.test(DS_data,2,type="Ljung-Box")

modelDS = arima(DS_data, order = c(0,0,2))

#######################################################################################

shifted_hfri_data <- shift.column(data = as_data, columns = 1, up = FALSE, newNames = "HFRI_1")

fitHFRI <- lm(HFRI_1 ~ RUS_Rf + RUS_1_Rf_1 + MXUS_Rf + MEM_Rf + SMB + HML + MOM + SBGC_Rf + SBWG_Rf + LHY_Rf + DEFSPR + FRBI_Rf + GSCI__Rf + VIX + Rf, shifted_hfri_data)
stepHFRI <- stepAIC(fitHFRI, direction="both")

shifted_ds_data <- shift.column(data = as_data, columns = 7, up = FALSE, newNames = "DS_1")

fitDS <- lm(DS_1 ~ RUS_Rf + RUS_1_Rf_1 + MXUS_Rf + MEM_Rf + SMB + HML + MOM + SBGC_Rf + SBWG_Rf + LHY_Rf + DEFSPR + FRBI_Rf + GSCI__Rf + VIX + Rf, shifted_ds_data)
stepDS <- stepAIC(fitDS, direction="both")

######################################################################################

# ACF/PAC on residuals: https://onlinecourses.science.psu.edu/stat510/node/72