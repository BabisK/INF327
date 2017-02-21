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
install.packages("rugarch")
library("rugarch")

#######################################################################################
# Initial configuration
#######################################################################################

# Load and name the data
header <- c("HFRI", "EH", "M", "RVA", "ED", "CA", "DS", "EMN", "MA", "RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf")
dates = seq(from = as.Date("1990-04-01", format='%Y-%m-%d'), to = as.Date("2005-12-01", format='%Y-%m-%d'), by = 'month')
as_data <- read.table("data_assignment.txt")
colnames(as_data) <- header
rownames(as_data) <- dates

# Plot all time series
png(filename = "images/DependentVariables.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(4,2))
for (col in c("HFRI", "EH", "M", "RVA", "ED", "CA", "DS", "EMN")){
  plot(ts(data = as_data[[paste(col)]], start = c(1990, 4), end = c(2004, 12), frequency = 12), ylab = paste(col))
}
par(mfrow=c(1,1))
dev.off()

#######################################################################################
# First part
#######################################################################################

# Analyze HFRI
HFRI_data <- ts(data = as_data$HFRI, start = c(1990, 4), end = c(2004, 12), frequency = 12)
ar <- ar(HFRI_data)
adfTest(HFRI_data, lags = ar$order, type = "c")

png(filename = "images/HFRIexploration.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(3,1))
plot(HFRI_data)
acf(HFRI_data)
pacf(HFRI_data)
par(mfrow = c(1,1))
dev.off()

Box.test(HFRI_data,1,type="Box-Pierce")
Box.test(HFRI_data,1,type="Ljung-Box")

modelHFRI = arima(HFRI_data, order = c(0,0,1))

png(filename = "images/HFRIarmaacf.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelHFRI))
pacf(residuals(modelHFRI))
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelHFRI),1,type="Box-Pierce")
Box.test(residuals(modelHFRI),1,type="Ljung-Box")

png(filename = "images/HFRIarmaqq.png", width = 2400, height = 2000, res = 250)
qqnorm(residuals(modelHFRI))
qqline(residuals(modelHFRI))
dev.off()

shapiro.test(residuals(modelHFRI))

png(filename = "images/HFRIarmaacf2.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelHFRI)^2)
pacf(residuals(modelHFRI)^2)
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelHFRI)^2,1,type="Box-Pierce")
Box.test(residuals(modelHFRI)^2,1,type="Ljung-Box")

# Analyze DS
DS_data <- ts(data = as_data$DS, start = c(1990, 4), end = c(2004, 12), frequency = 12)
ar <- ar(DS_data)
adfTest(DS_data, lags = ar$order, type = "c")

png(filename = "images/DSexploration.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(3,1))
plot(DS_data)
acf(DS_data)
pacf(DS_data)
par(mfrow = c(1,1))
dev.off()

Box.test(DS_data,2,type="Box-Pierce")
Box.test(DS_data,2,type="Ljung-Box")

modelDS = arima(DS_data, order = c(0,0,2))

png(filename = "images/DSarmaacf.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelDS))
pacf(residuals(modelDS))
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelDS),1,type="Box-Pierce")
Box.test(residuals(modelDS),1,type="Ljung-Box")

png(filename = "images/DSarmaqq.png", width = 2400, height = 2000, res = 250)
qqnorm(residuals(modelDS))
qqline(residuals(modelDS))
dev.off()

shapiro.test(residuals(modelDS))

png(filename = "images/DSarmaacf2.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelDS)^2)
pacf(residuals(modelDS)^2)
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelDS)^2,1,type="Box-Pierce")
Box.test(residuals(modelDS)^2,1,type="Ljung-Box")

#######################################################################################
# Second part
#######################################################################################

# Shift variables
shifted_data <- shift.column(data = as_data, columns = c("RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf"))

# Fit HFRI
fitHFRI <- lm(HFRI ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted, shifted_data)
stepHFRI <- stepAIC(fitHFRI, direction="both", trace = 0)
summary(stepHFRI)
stepHFRI$anova

# Fit DS
fitDS <- lm(DS ~ RUS_Rf.Shifted + RUS_1_Rf_1.Shifted + MXUS_Rf.Shifted + MEM_Rf.Shifted + SMB.Shifted + HML.Shifted + MOM.Shifted + SBGC_Rf.Shifted + SBWG_Rf.Shifted + LHY_Rf.Shifted + DEFSPR.Shifted + FRBI_Rf.Shifted + GSCI__Rf.Shifted + VIX.Shifted + Rf.Shifted, shifted_data)
stepDS <- stepAIC(fitDS, direction="both", trace = 0)
summary(stepDS)
stepDS$anova

#######################################################################################
# Third part
#######################################################################################

# Analyse HFRI
png(filename = "images/HFRIlm.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,2))
plot(stepHFRI)
par(mfrow=c(1,1))
dev.off()

# Normality
shapiro.test(residuals(stepHFRI))

# Autocorrelation
png(filename = "images/HFRIlmacf.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,1))
acf(residuals(stepHFRI))
pacf(residuals(stepHFRI))
par(mfrow=c(1,1))
dev.off()

Box.test(residuals(stepHFRI),1,type="Box-Pierce")
Box.test(residuals(stepHFRI),1,type="Ljung-Box")

# Heteroscedasticity
png(filename = "images/HFRIlm2acf.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,1))
acf(residuals(stepHFRI)^2)
pacf(residuals(stepHFRI)^2)
par(mfrow=c(1,1))
dev.off()

bptest(stepHFRI)

#####

fixHFRI<-arima(shifted_data$HFRI, order = c(0,0,10), fixed=c(NA, 0, 0, 0, 0, 0, 0, 0, 0, NA, NA, NA, NA, NA, NA), xreg = cbind(shifted_data$RUS_1_Rf_1.Shifted, shifted_data$MEM_Rf.Shifted, shifted_data$DEFSPR.Shifted, shifted_data$VIX.Shifted)) 
par(mfrow=c(2,1))
acf(residuals(fixHFRI))
pacf(residuals(fixHFRI))
par(mfrow=c(1,1))

par(mfrow=c(2,1))
acf(residuals(fixHFRI)^2)
pacf(residuals(fixHFRI)^2)
par(mfrow=c(1,1))

Box.test(residuals(stepHFRI),1,type="Box-Pierce")


######

# Analyse DS
png(filename = "images/DSlm.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,2))
plot(stepDS)
par(mfrow=c(1,1))
dev.off()

# Normality
shapiro.test(residuals(stepDS))

# Autocorrelation
png(filename = "images/DSlmacf.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,1))
acf(residuals(stepDS))
pacf(residuals(stepDS))
par(mfrow=c(1,1))
dev.off()

Box.test(residuals(stepDS),1,type="Box-Pierce")
Box.test(residuals(stepDS),1,type="Ljung-Box")

# Heteroscedasticity
png(filename = "images/DSlm2acf.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,1))
acf(residuals(stepDS)^2)
pacf(residuals(stepDS)^2)
par(mfrow=c(1,1))
dev.off()

bptest(stepDS)

######################################################################################

# ACF/PAC on residuals: https://onlinecourses.science.psu.edu/stat510/node/72
residualsHFRI <- ts(stepHFRI$residuals, start = c(1990, 4), end = c(2004, 12), frequency = 12)
par(mfrow=c(3,1))
plot(residualsHFRI)
acf(residualsHFRI)
pacf(residualsHFRI)
par(mfrow=c(1,1))
adfTest(residualsHFRI, lags = 10, type = "ct")


fit.spec <- ugarchspec(variance.model = list(model = "fGARCH",  garchOrder = c(1, 1), submodel="GARCH"), mean.model = list(armaOrder = c(0, 10), include.mean = TRUE, external.regressors = cbind(shifted_data$RUS_1_Rf_1.Shifted, shifted_data$MEM_Rf.Shifted, shifted_data$DEFSPR.Shifted, shifted_data$VIX.Shifted)), distribution.model = "norm")
fit <- ugarchfit(data = shifted_data$HFRI, spec = fit.spec)
