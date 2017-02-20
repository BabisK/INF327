install.packages("fUnitRoots")
install.packages("useful")
install.packages("lmtest")# for Breush Pagan Test (bptest) to detect hoscedasticity
install.packages("forecast")
install.packages("fGarch")
library("fUnitRoots")
library("useful")
library("MASS")
library("lmtest") # for Breush Pagan Test (bptest) to detect hoscedasticity
library("forecast")
library("fGarch")

# Load and name the data
header <- c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN", "MA", "RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf")
dates = seq(from = as.Date("1990-04-01", format='%Y-%m-%d'), to = as.Date("2005-12-01", format='%Y-%m-%d'), by = 'month')
as_data <- read.table("data_assignment.txt")
colnames(as_data) <- header
rownames(as_data) <- dates

# Plot all time series
png(filename = "images/dependentVariables.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(4,2))
for (col in c("HFRI ", "EH", "M", "RVA", "ED", "CA", "DS", "EMN")){
  plot(ts(data = as_data[[paste(col)]], start = c(1990, 4), end = c(2004, 12), frequency = 12), ylab = paste(col))
}
par(mfrow=c(1,1))
dev.off()

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

# Fit HFRI
shifted_hfri_data <- shift.column(data = as_data, columns = 1, up = FALSE, newNames = "HFRI_1")

fitHFRI <- lm(HFRI_1 ~ RUS_Rf + RUS_1_Rf_1 + MXUS_Rf + MEM_Rf + SMB + HML + MOM + SBGC_Rf + SBWG_Rf + LHY_Rf + DEFSPR + FRBI_Rf + GSCI__Rf + VIX + Rf, shifted_hfri_data)
stepHFRI <- stepAIC(fitHFRI, direction="both")
par(mfrow=c(2,2))
plot(stepHFRI)
par(mfrow=c(1,1))

shapiro.test(residuals(stepHFRI))
par(mfrow=c(2,1))
acf(residuals(stepHFRI))
pacf(residuals(stepHFRI))
par(mfrow=c(1,1))

fixHFRI<-arima(shifted_hfri_data$`HFRI `, order = c(20,0,0), xreg = cbind.data.frame(shifted_hfri_data$RUS_1_Rf_1, shifted_hfri_data$MEM_Rf, shifted_hfri_data$DEFSPR, shifted_hfri_data$VIX))
par(mfrow=c(2,1))
acf(residuals(fixHFRI))
pacf(residuals(fixHFRI))
par(mfrow=c(1,1))


# FIT DS
shifted_ds_data <- shift.column(data = as_data, columns = 7, up = FALSE, newNames = "DS_1")

fitDS <- lm(DS_1 ~ RUS_Rf + RUS_1_Rf_1 + MXUS_Rf + MEM_Rf + SMB + HML + MOM + SBGC_Rf + SBWG_Rf + LHY_Rf + DEFSPR + FRBI_Rf + GSCI__Rf + VIX + Rf, shifted_ds_data)
stepDS <- stepAIC(fitDS, direction="both")
par(mfrow=c(2,2))
plot(stepDS)
par(mfrow=c(1,1))

######################################################################################

# ACF/PAC on residuals: https://onlinecourses.science.psu.edu/stat510/node/72
residualsHFRI <- ts(stepHFRI$residuals, start = c(1990, 4), end = c(2004, 12), frequency = 12)
par(mfrow=c(3,1))
plot(residualsHFRI)
acf(residualsHFRI)
pacf(residualsHFRI)
par(mfrow=c(1,1))
adfTest(residualsHFRI, lags = 10, type = "ct")
