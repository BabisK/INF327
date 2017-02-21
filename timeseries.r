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

# MA(1)
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

# MA(12)
modelHFRI12 = arima(HFRI_data, order = c(0,0,12))

png(filename = "images/HFRI12armaacf.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelHFRI12))
pacf(residuals(modelHFRI12))
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelHFRI12),12,type="Box-Pierce")
Box.test(residuals(modelHFRI12),12,type="Ljung-Box")

png(filename = "images/HFRI12armaqq.png", width = 2400, height = 2000, res = 250)
qqnorm(residuals(modelHFRI12))
qqline(residuals(modelHFRI12))
dev.off()

shapiro.test(residuals(modelHFRI12))

png(filename = "images/HFRI12armaacf2.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelHFRI12)^2)
pacf(residuals(modelHFRI12)^2)
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelHFRI12)^2,1,type="Box-Pierce")
Box.test(residuals(modelHFRI12)^2,1,type="Ljung-Box")

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

# MA(2)
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

# MA(12)
modelDS12 = arima(DS_data, order = c(0,0,12))

png(filename = "images/DS12armaacf.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelDS12))
pacf(residuals(modelDS12))
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelDS12),1,type="Box-Pierce")
Box.test(residuals(modelDS12),1,type="Ljung-Box")

png(filename = "images/DS12armaqq.png", width = 2400, height = 2000, res = 250)
qqnorm(residuals(modelDS12))
qqline(residuals(modelDS12))
dev.off()

shapiro.test(residuals(modelDS12))

png(filename = "images/DS12armaacf2.png", width = 2400, height = 2000, res = 250)
par(mfrow = c(2,1))
acf(residuals(modelDS12)^2)
pacf(residuals(modelDS12)^2)
par(mfrow = c(1,1))
dev.off()

Box.test(residuals(modelDS12)^2,1,type="Box-Pierce")
Box.test(residuals(modelDS12)^2,1,type="Ljung-Box")

#######################################################################################
# Second part
#######################################################################################

# Shift variables
shifted_data <- shift.column(data = as_data, up = FALSE, columns = c("RUS_Rf", "RUS_1_Rf_1", "MXUS_Rf", "MEM_Rf", "SMB", "HML", "MOM", "SBGC_Rf", "SBWG_Rf", "LHY_Rf", "DEFSPR", "FRBI_Rf", "GSCI__Rf", "VIX", "Rf"))
shifted_data_new <- shifted_data[c(177,178,179,180,181,182,183,184,185,186,187,188),]
shifted_data <- shifted_data[-c(188,187,186,185,184,183,182,181,180,179, 178, 177),]

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

Box.test(residuals(stepHFRI)^2,2,type="Box-Pierce")
Box.test(residuals(stepHFRI)^2,2,type="Ljung-Box")

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

Box.test(residuals(stepDS),8,type="Box-Pierce")
Box.test(residuals(stepDS),8,type="Ljung-Box")

# Heteroscedasticity
png(filename = "images/DSlm2acf.png", width = 2400, height = 2000, res = 250)
par(mfrow=c(2,1))
acf(residuals(stepDS)^2)
pacf(residuals(stepDS)^2)
par(mfrow=c(1,1))
dev.off()

Box.test(residuals(stepDS)^2,1,type="Box-Pierce")
Box.test(residuals(stepDS)^2,1,type="Ljung-Box")

#######################################################################################
# Bellow ARMA models are fitted on the residuals of regression. This is not needed according to
# the tests performed but is left in the code for the sake of completeness
#######################################################################################

# Correcting HFRI
#fixHFRI<-arima(shifted_data$HFRI, order = c(0,0,1), xreg = cbind(shifted_data$MEM_Rf.Shifted)) 
#png(filename = "images/HFRIlmarmaacf.png", width = 2400, height = 2000, res = 250)
#par(mfrow=c(2,1))
#acf(residuals(fixHFRI))
#pacf(residuals(fixHFRI))
#par(mfrow=c(1,1))
#dev.off()

#Box.test(residuals(fixHFRI),1,type="Box-Pierce")
#Box.test(residuals(fixHFRI),1,type="Ljung-Box")

#png(filename = "images/HFRIlmarmaacf2.png", width = 2400, height = 2000, res = 250)
#par(mfrow=c(2,1))
#acf(residuals(fixHFRI)^2)
#pacf(residuals(fixHFRI)^2)
#par(mfrow=c(1,1))
#dev.off()

#Box.test(residuals(fixHFRI)^2,1,type="Box-Pierce")
#Box.test(residuals(fixHFRI)^2,1,type="Ljung-Box")

#png(filename = "images/HFRIlmarmaqq.png", width = 2400, height = 2000, res = 250)
#qqnorm(residuals(fixHFRI))
#qqline(residuals(fixHFRI))
#dev.off()

#shapiro.test(residuals(fixHFRI))

# Correcting DS
#fixDS<-arima(shifted_data$DS, order = c(0,0,1), xreg = cbind(shifted_data$RUS_Rf.Shifted, shifted_data$RUS_1_Rf_1.Shifted, shifted_data$HML.Shifted, shifted_data$MOM.Shifted, shifted_data$LHY_Rf.Shifted, shifted_data$DEFSPR.Shifted, shifted_data$VIX.Shifted, shifted_data$Rf.Shifted)) 
#png(filename = "images/DSlmarmaacf.png", width = 2400, height = 2000, res = 250)
#par(mfrow=c(2,1))
#acf(residuals(fixDS))
#pacf(residuals(fixDS))
#par(mfrow=c(1,1))
#dev.off()

#Box.test(residuals(fixDS),1,type="Box-Pierce")
#Box.test(residuals(fixDS),1,type="Ljung-Box")

#png(filename = "images/DSlmarmaacf2.png", width = 2400, height = 2000, res = 250)
#par(mfrow=c(2,1))
#acf(residuals(fixDS)^2)
#pacf(residuals(fixDS)^2)
#par(mfrow=c(1,1))
#dev.off()

#Box.test(residuals(fixDS)^2,1,type="Box-Pierce")
#Box.test(residuals(fixDS)^2,1,type="Ljung-Box")

#png(filename = "images/DSlmarmaqq.png", width = 2400, height = 2000, res = 250)
#qqnorm(residuals(fixDS))
#qqline(residuals(fixDS))
#dev.off()

#shapiro.test(residuals(fixDS))

# ACF/PAC on residuals: https://onlinecourses.science.psu.edu/stat510/node/72
#residualsHFRI <- ts(stepHFRI$residuals, start = c(1990, 4), end = c(2004, 12), frequency = 12)
#par(mfrow=c(3,1))
#plot(residualsHFRI)
#acf(residualsHFRI)
#pacf(residualsHFRI)
#par(mfrow=c(1,1))
#adfTest(residualsHFRI, lags = 10, type = "ct")

#fit.spec <- ugarchspec(variance.model = list(model = "fGARCH",  garchOrder = c(1, 1), submodel="GARCH"), mean.model = list(armaOrder = c(0, 10), include.mean = TRUE, external.regressors = cbind(shifted_data$RUS_1_Rf_1.Shifted, shifted_data$MEM_Rf.Shifted, shifted_data$DEFSPR.Shifted, shifted_data$VIX.Shifted)), distribution.model = "norm")
#fit <- ugarchfit(data = shifted_data$HFRI, spec = fit.spec)

#######################################################################################
# Forecasting
#######################################################################################

# HFRI
forecastMA1HFRI <- predict(modelHFRI, 12)
forecastMA12HFRI <- predict(modelHFRI12, 12)
forecastLMHFRI <- predict(stepHFRI, shifted_data_new, se.fit = TRUE)

png(filename = "images/forecastMA1HFRI.png", width = 2400, height = 2000, res = 250)
UL <- forecastMA1HFRI$pred + forecastMA1HFRI$se
LL <- forecastMA1HFRI$pred - forecastMA1HFRI$se
ts.plot(ts(shifted_data$HFRI, start = c(1990, 5), end = c(2004, 12), frequency = 12), forecastMA1HFRI$pred, ts(shifted_data_new$HFRI, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="MA(1) HFRI")
lines(forecastMA1HFRI$pred, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

png(filename = "images/forecastMA12HFRI.png", width = 2400, height = 2000, res = 250)
UL <- forecastMA12HFRI$pred + forecastMA12HFRI$se
LL <- forecastMA12HFRI$pred - forecastMA12HFRI$se
ts.plot(ts(shifted_data$HFRI, start = c(1990, 5), end = c(2004, 12), frequency = 12), forecastMA12HFRI$pred, ts(shifted_data_new$HFRI, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="MA(12) HFRI")
lines(forecastMA12HFRI$pred, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

png(filename = "images/forecastLMHFRI.png", width = 2400, height = 2000, res = 250)
UL <- forecastLMHFRI$fit + forecastLMHFRI$se.fit
LL <- forecastLMHFRI$fit - forecastLMHFRI$se.fit
p <- ts(forecastLMHFRI$fit, start = c(2005, 1), end = c(2005, 12), frequency = 12)
ts.plot(ts(shifted_data$HFRI, start = c(1990, 5), end = c(2004, 12), frequency = 12), p, ts(shifted_data_new$HFRI, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="LM HFRI")
lines(p, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

msfeMA1HFRI <- sum((forecastMA1HFRI$pred - shifted_data_new$HFRI)^2)/12
msfeMA12HFRI <- sum((forecastMA12HFRI$pred - shifted_data_new$HFRI)^2)/12
msfeLMHFRI <- sum((forecastLMHFRI$fit - shifted_data_new$HFRI)^2)/12

hitMA1HFRI <- sum((forecastMA1HFRI$pred * shifted_data_new$HFRI) > 0)/12
hitMA12HFRI <- sum((forecastMA12HFRI$pred * shifted_data_new$HFRI) > 0)/12
hitLMHFRI <- sum((forecastLMHFRI$fit * shifted_data_new$HFRI) > 0)/12

# DS

forecastMA1DS <- predict(modelDS, 12)
forecastMA12DS <- predict(modelDS12, 12)
forecastLMDS <- predict(stepDS, shifted_data_new, se.fit = TRUE)

png(filename = "images/forecastMA1DS.png", width = 2400, height = 2000, res = 250)
UL <- forecastMA1DS$pred + forecastMA1DS$se
LL <- forecastMA1DS$pred - forecastMA1DS$se
ts.plot(ts(shifted_data$DS, start = c(1990, 5), end = c(2004, 12), frequency = 12), forecastMA1DS$pred, ts(shifted_data_new$DS, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="MA(2) DS")
lines(forecastMA1DS$pred, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

png(filename = "images/forecastMA12DS.png", width = 2400, height = 2000, res = 250)
UL <- forecastMA12DS$pred + forecastMA12DS$se
LL <- forecastMA12DS$pred - forecastMA12DS$se
ts.plot(ts(shifted_data$DS, start = c(1990, 5), end = c(2004, 12), frequency = 12), forecastMA12DS$pred, ts(shifted_data_new$DS, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="MA(12) DS")
lines(forecastMA12DS$pred, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

png(filename = "images/forecastLMDS.png", width = 2400, height = 2000, res = 250)
UL <- forecastLMDS$fit + forecastLMDS$se.fit
LL <- forecastLMDS$fit - forecastLMDS$se.fit
p <- ts(forecastLMDS$fit, start = c(2005, 1), end = c(2005, 12), frequency = 12)
ts.plot(ts(shifted_data$DS, start = c(1990, 5), end = c(2004, 12), frequency = 12), p, ts(shifted_data_new$DS, start = c(2005, 1), end = c(2005, 12), frequency = 12), main="LM DS")
lines(p, col="red")
lines(UL, col="blue", lty="dashed")
lines(LL, col="blue", lty="dashed")
dev.off()

msfeMA1DS <- sum((forecastMA1DS$pred - shifted_data_new$DS)^2)/12
msfeMA12DS <- sum((forecastMA12DS$pred - shifted_data_new$DS)^2)/12
msfeLMDS <- sum((forecastLMDS$fit - shifted_data_new$DS)^2)/12

hitMA1DS <- sum((forecastMA1DS$pred * shifted_data_new$DS) > 0)/12
hitMA12DS <- sum((forecastMA12DS$pred * shifted_data_new$DS) > 0)/12
hitLMDS <- sum((forecastLMDS$fit * shifted_data_new$DS) > 0)/12
