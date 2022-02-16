# Supplementary code to:

# How to estimate carabid biomass? – An valuation of size-weight models for ground beetles (Coleoptera: Carabidae) and perspectives for further improvement

# submitted to Journal of Insect Conservation
# Fabio Weiss, Andreas Linde
# University for Sustainable Development Eberswalde, fabio.weiss@hnee.de

# compiled with R version 4.1.2. 


# Packagelist
library(lme4)
library(DHARMa)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(nlme)
library(writexl)
library(lmerTest)
library(car)
library(effects)
library(multcomp)

# load assembled table
weight_data <- read.csv2("compiled_data.csv")

# data can be found here: 
# Schultz, R. (1996) Die Laufkäfer als Indikatoren der Renaturierung des Salzgrünlandes im Ostseebereich Vorpommerns. Culliver-Verlag, Göttingen.
# Booij, K., den Nijs, L., Heijermann, T., Jorritsma, I., Lock, C., Noorlander, J. (1994) Size and weight of carabid beetles: ecological applications. Proc Exper Appl Entomol.

weight_data$ln_size <- log(weight_data$av_size_mm)
weight_data$ln_weight <- log(weight_data$weight_mg)

exclude_subfamily <- c("Omophroninae", "Loricerinae", "Broscinae", "Cicindelinae")

weight_data2 <- weight_data[!(weight_data$subfamily %in% exclude_subfamily) ,]
weight_data$subfamily <- as.factor(weight_data$subfamily)

#figure 1
par(mfrow= c(1,2))
plot(weight_data$weight_mg ~ weight_data$av_size_mm, pch=20, main= "a", ylab="weight (mg)", xlab="body length (mm)", col="grey",cex.lab=1.4, cex.axis=1.1, cex=1.2)
pred_size <- seq(min(weight_data$av_size_mm), max(weight_data$av_size_mm), length.out = 100)
lines(pred_size, szyszko_weight(pred_size), lty="dotted", lwd=2)
lines(pred_size, booij_weight(pred_size), lty="twodash", lwd=2)

plot(weight_data$ln_weight ~ weight_data$ln_size, pch=20, main= "b", ylab="ln(weight)", xlab="ln(body length)",col="grey",cex.lab=1.4, cex.axis=1.1, cex=1.2)
pred_ln_size <- log(pred_size)
lines(pred_ln_size, log(szyszko_weight(pred_size)),lty="dotted", lw=2)
lines(pred_ln_size, log(booij_weight(pred_size)),lty="twodash", lw=2)

weights_booij <- weight_data[weight_data$source_weight == "booij",]
weights_schultz <- weight_data[weight_data$source_weight == "schultz",]

# figure 2
par(mfrow= c(2,2))
hist(weights_booij$av_size_mm, breaks=c(1:35), main= "Data of Booij et al. (1994), (training set)", ylim = c(0,25), xlab= "body length (mm)",cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
plot(weights_booij$av_size_mm, weights_booij$weight_mg, ylab="weight (mg)", xlab = "body length (mm)", ylim= c(0,900), xlim = c(0,35),cex.lab=1.4, cex.axis=1.4)

hist(weights_schultz$av_size_mm, breaks=c(1:35), main= "Data of Schultz (1996), (validation set)", ylim = c(0,25), xlab= "body length (mm)",cex.lab=1.4, cex.axis=1.4, cex.main=1.3)
plot(weights_schultz$av_size_mm, weights_schultz$weight_mg, ylab="weight (mg)", xlab = "body length (mm)", ylim= c(0,900), xlim = c(0,35),cex.lab=1.4, cex.axis=1.4)
par(mfrow= c(1,1))

# fitting models
m0<- lm(ln_weight ~ ln_size, data= weights_booij)
m1<- lm(ln_weight ~ ln_size + subfamily, data= weights_booij)
m1r <-lmer(ln_weight ~ ln_size + (1|subfamily), weights_booij, REML = F)
m2 <- lm(ln_weight ~ ln_size * subfamily, data= weights_booij)
# in the manuscript the models are referred to as m_base, m_fixed, m_mixed, m_inter (in this order)

tab_model(m0,m1,m1r,m2, dv.labels= c("m_base","m_fixed","m_mixed","m_inter"), show.df=T, show.dev=T, digits = 4)

# Diagnostics with DHARMa package
resids <- simulateResiduals(fittedModel = m0, plot=F)
plotQQunif(resids)

resids <- simulateResiduals(fittedModel = m1, plot=F)
plotQQunif(resids)

resids <- simulateResiduals(fittedModel = m1r, plot=F)
plotQQunif(resids)

resids <- simulateResiduals(fittedModel = m2, plot=F)
plotQQunif(resids)

# implementing Syzsko (1983) as a function
szyszko_weight<- function(size){
  ln_size <- log(size)
  ln_weight <- -8.92804283 + 2.5554921* ln_size
  weight_mg <- exp(ln_weight)*1000
  return(weight_mg)
}

# implementing Booij et al. (1994) as a function
booij_weight<- function(size){
  log_size <- log10(size)
  log_weight <- (-1.3 + 2.95 *log_size)
  weight_mg <- 10^(log_weight)
  return(weight_mg)
}

# predicting weights for the Schultz (1996) dataset
 
weights_schultz$szyszko_pred <- szyszko_weight(weights_schultz$av_size_mm)

weights_schultz$booij_pred <- booij_weight(weights_schultz$av_size_mm)

pred<- predict(m0, newdata = weights_schultz, type="response")
weights_schultz$m0_pred <- exp(pred)

pred<- predict(m1, newdata = weights_schultz, type="response")
weights_schultz$m1_pred <- exp(pred)

pred<- predict(m1r, newdata = weights_schultz, re.form=NA)
weights_schultz$m1r_pred <- exp(pred)

pred<- predict(m2, newdata = weights_schultz, type="response")
weights_schultz$m2_pred <- exp(pred)


# calculated error in percent bodyweight
weights_schultz$szyszko_off <- (weights_schultz$szyszko_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

weights_schultz$booij_off <- (weights_schultz$booij_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

weights_schultz$m0_off <- (weights_schultz$m0_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

weights_schultz$m1_off <- (weights_schultz$m1_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

weights_schultz$m1r_off <- (weights_schultz$m1r_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

weights_schultz$m2_off <- (weights_schultz$m2_pred - weights_schultz$weight_mg)/weights_schultz$weight_mg *100

# Figure 3, relative deviation graphs (sensu Mitchell, 1997)
par(mfrow= c(3,2))

plot(weights_schultz$szyszko_off ~ weights_schultz$av_size_mm,ylim=c(-70,220), ylab= "", xlab="",pch=20, cex.lab=1.5, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0,  col="red")
title("m_Szyszko", line = 1, cex=1.4)

plot(weights_schultz$booij_off ~ weights_schultz$av_size_mm, ylim=c(-70,220),  ylab= "", xlab="",pch=20, cex.lab=1.4, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0, col="red")
title("m_Booij", line = 1, cex=1.4)

plot(weights_schultz$m0_off ~ weights_schultz$av_size_mm, ylim=c(-70,220), ylab= "percent of error", xlab="",pch=20, cex.lab=2, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0, col="red")
title("m_base", line = 1, cex=1.4)
?title

plot(weights_schultz$m1_off ~ weights_schultz$av_size_mm, ylim=c(-70,220), ylab= "", xlab="",pch=20, cex.lab=1.4, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0, col="red")
title("m_fixed", line = 1, cex=1.4)

plot(weights_schultz$m1r_off ~ weights_schultz$av_size_mm,ylim=c(-70,220), ylab= "", xlab="body length (mm)",pch=20, cex.lab=2, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0, col="red")
title("m_mixed", line = 1, cex=1.4)

plot(weights_schultz$m2_off ~ weights_schultz$av_size_mm,ylim=c(-70,220), ylab= "", xlab="",pch=20, cex.lab=1.4, cex=1.5, cex.axis=1.5, cex.main=1.2)
abline(h=0, col="red")
title("m_inter", line = 1, cex=1.4)

# observed vs. predicted regression (sensu Pineiro et al., 2008)
ref_weights <- data.frame(
  model= "reference",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$weight_mg, 
  weight_off= 0)

szyzsko_preds <- data.frame(
  model= "szyszko",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$szyszko_pred, 
  weight_off= weights_schultz$szyszko_off)

booij_preds <- data.frame(
  model= "booij",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$booij_pred, 
  weight_off= weights_schultz$booij_off)

m0_preds <- data.frame(
  model= "m_base",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$m0_pred, 
  weight_off= weights_schultz$m0_off)

m1_preds <- data.frame(
  model= "m_fixed",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$m1_pred, 
  weight_off= weights_schultz$m1_off)

m1r_preds <- data.frame(
  model= "m_mixed",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$m1r_pred, 
  weight_off= weights_schultz$m1r_off)

m2_preds <- data.frame(
  model= "m_inter",
  obs_weight=weights_schultz$weight_mg,
  pred_weight= weights_schultz$m2_pred, 
  weight_off= weights_schultz$m2_off)

resids_data <- rbind(ref_weights,szyzsko_preds, booij_preds, m0_preds, m1_preds, m1r_preds, m2_preds)
resids_data$ln_pred_weight <- log(resids_data$pred_weight)
resids_data$ln_obs_weight <- log(resids_data$obs_weight)
resids_data$model <- as.factor(resids_data$model)
resids_data <- within(resids_data, model <- relevel(model, ref = "reference"))

err_mod <- lm(obs_weight ~ pred_weight * model, data=resids_data)

resids <- simulateResiduals(fittedModel = err_mod, plot=T)
plotQQunif(resids)

log_err_mod <- lm(ln_obs_weight ~ ln_pred_weight * model, data=resids_data)

resids <- simulateResiduals(fittedModel = log_err_mod, plot=T)
plotQQunif(resids)

# Figure 4
ref_data <- resids_data[resids_data$model== "reference",]
par(mfrow= c(6,2))

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="szyszko")
preds       <- predict(log_err_mod, newdata = pred_data, intervals="confidence")

data <- resids_data[resids_data$model== "szyszko",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=2)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("log-transformed", 3, line=2,cex=1.3, at=3.5)
mtext("m_Szyszko", 3, line=-2,cex=1.3, at=1.5)
mtext("intercept: ***", 3, line=-4,cex=1, at=1.5)
mtext("    slope: ***", 3, line=-5.5,cex=1, at=1.5)
mtext("R² = 0.9515", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="szyszko")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "szyszko",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=4)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("not transformed", 3, line=2,cex=1.3, at=450)
mtext("m_Szyszko", 3, line=-2,cex=1.3, at=150)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=150)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=150)
mtext("R² = 0.8823", 1, line=-3, cex = 1.3, at=700)

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="booij")
preds       <- predict(log_err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "booij",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=2)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("m_Booij", 3, line=-2,cex=1.3, at=1)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=1)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=1)
mtext("R² = 0.9515", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="booij")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "booij",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=2)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("m_Booij", 3, line=-2,cex=1.3, at=150)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=170)
mtext("    slope: ***", 3, line=-5.5,cex=1, at=160)
mtext("R² = 0.8539", 1, line=-3, cex = 1.3, at=1000)

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="m_base")
preds       <- predict(log_err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_base",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="log(observed weights)", xlab="",  cex.lab=2, cex= 1.4, cex.axis=1.5, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=3)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("m_base", 3, line=-2,cex=1.3, at=1)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=1)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=1)
mtext("R² = 0.9515", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="m_base")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_base",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="observed weights (mg)", xlab="",  cex.lab=2, cex= 1.4, cex.axis=1.5, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=2)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("m_base", 3, line=-2,cex=1.3, at=150)
mtext("intercept: *", 3, line=-4,cex=1, at=150)
mtext("    slope: ***", 3, line=-5.5,cex=1, at=150)
mtext("R² = 0.8516", 1, line=-3, cex = 1.3, at=1000)

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="m_fixed")
preds       <- predict(log_err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_fixed",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=3)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("m_fixed", 3, line=-2,cex=1.3, at=1)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=1)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=1)
mtext("R² = 0.9516", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="m_fixed")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_fixed",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=1.5, cex= 1.4, cex.axis=1.5, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=2)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("m_fixed", 3, line=-2,cex=1.3, at=150)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=150)
mtext("    slope: ***", 3, line=-5.5,cex=1, at=150)
mtext("R² = 0.8584", 1, line=-3, cex = 1.3, at=900)

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="m_mixed")
preds       <- predict(log_err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_mixed",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=2, cex= 1.4, cex.axis=1.5, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=3)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("m_mixed", 3, line=-2,cex=1.3, at=1)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=1)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=1)
mtext("R² = 0.9515", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="m_mixed")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_mixed",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="", xlab="", cex.lab=2, cex= 1.4, cex.axis=1.5, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=2)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("m_mixed", 3, line=-2,cex=1.3, at=150)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=150)
mtext("    slope: ***", 3, line=-5.5,cex=1, at=150)
mtext("R² = 0.8558", 1, line=-3, cex = 1.3, at=900)

ln_pred_weight <- seq(min(resids_data$ln_pred_weight), max(resids_data$ln_pred_weight), length.out = 100)
ln_obs_weight  <- seq(min(resids_data$ln_obs_weight), max(resids_data$ln_obs_weight), length.out = 100)
pred_data   <- data.frame(ln_pred_weight=ln_pred_weight, model="m_inter")
preds       <- predict(log_err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_inter",]
r2_mod <- lm(ln_obs_weight ~ ln_pred_weight, data=data)
summary(r2_mod)

plot(ln_obs_weight ~ ln_pred_weight,pch=20, data=data, type="p", ylab="", xlab="log(predicted weights)", cex.lab=2, cex= 1.4, cex.axis=1.4, col="grey")
lines(-1:9 , -1:9,lty= 2,lwd=3)
lines(preds$fit ~ ln_pred_weight, lwd=2, col="red")
mtext("m_inter", 3, line=-2,cex=1.3, at=1)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=1)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=1)
mtext("R² = 0.9520", 1, line=-3, cex = 1.3, at=6)

pred_weight <- seq(min(resids_data$pred_weight), max(resids_data$pred_weight), length.out = 100)
obs_weight  <- seq(min(resids_data$obs_weight), max(resids_data$obs_weight), length.out = 100)
pred_data   <- data.frame(pred_weight=pred_weight, model="m_inter")
preds       <- predict(err_mod, newdata = pred_data, type="response", se.fit = T)

data <- resids_data[resids_data$model== "m_inter",]
r2_mod <- lm(obs_weight ~ pred_weight, data=data)
summary(r2_mod)

plot(obs_weight ~ pred_weight,pch=20, data=data, type="p", ylab="", xlab="predicted weights (mg)", cex.lab=2, cex= 1.4, cex.axis=1.4, col="black")
lines(1:1200 , 1:1200,lty= 2,lwd=2)
lines(preds$fit ~ pred_weight, lwd=2, col="red")
mtext("m_inter", 3, line=-2,cex=1.3, at=150)
mtext("intercept: n.s.", 3, line=-4,cex=1, at=150)
mtext("    slope: n.s.", 3, line=-5.5,cex=1, at=150)
mtext("R² = 0.9052", 1, line=-3, cex = 1.3, at=700)

pred_ln_size <- seq(min(weights_clean$ln_size), max(weights_clean$ln_size), length.out = 100)

pred_trechinae <- c(rep("Trechinae",100))
pred_harpalinae <- c(rep("Harpalinae", 100))
pred_carabinae <- c(rep("Carabinae", 100))
pred_scaritinae <- c(rep("Scaritinae ", 100))
pred_elaphrinae <- c(rep("Elaphrinae", 100))
pred_nebriinae <- c(rep("Nebriinae", 100))

pred_trechinae <- data.frame(ln_size=pred_ln_size, subfamily=pred_trechinae)
pred_harpalinae <- data.frame(ln_size=pred_ln_size, subfamily=pred_harpalinae)
pred_carabinae <- data.frame(ln_size=pred_ln_size, subfamily=pred_carabinae)
pred_scaritinae <-data.frame(ln_size=pred_ln_size, subfamily=pred_scaritinae)
pred_elaphrinae <- data.frame(ln_size=pred_ln_size, subfamily=pred_elaphrinae)
pred_nebriinae <-data.frame(ln_size=pred_ln_size, subfamily=pred_nebriinae)


# figure S1
par(mfrow= c(1,2))
preds<- predict(m0, newdata = pred_harpalinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="weight (mg)", xlab="body length (mm)", main="", col="grey", cex=0.8,cex.lab=1.4, cex.axis=1.1)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="ln(weight)", xlab="ln(body length)", main="", col="grey", cex=0.8,cex.lab=1.4, cex.axis=1.1)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

# figure S2
par(mfrow= c(6,2))
preds<- predict(m1, newdata = pred_harpalinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Harpalinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m1, newdata = pred_carabinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Carabinae", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m1, newdata = pred_elaphrinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="weight (mg)", xlab="", main="Elaphrinae", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="ln(weight)", xlab="", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m1, newdata = pred_nebriinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Nebriinae", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m1, newdata = pred_scaritinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Scaritinae", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m1, newdata = pred_trechinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="body length (mm)", main="Trechinae", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="ln(body length)", main="", col="grey", cex=0.8,cex.lab=1.6,cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

# figure S3
par(mfrow= c(2,1))
preds<- predict(m1r, newdata = pred_harpalinae, re.form=NA)

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="weight (mg)", xlab="body length (mm)", main="", col="grey", cex=0.8,cex.lab=1.4, cex.axis=1.1)
mtext("a", 3, line=-2,cex=1.3, at=7)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="ln(weight)", xlab="ln(body length)", main="", col="grey", cex=0.8,cex.lab=1.4, cex.axis=1.1)
mtext("b", 3, line=-2,cex=1.3, at=1.5)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds,col="red", lw=2)

# figure S4
par(mfrow= c(6,2))
preds<- predict(m2, newdata = pred_harpalinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Harpalinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m2, newdata = pred_carabinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Carabinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m2, newdata = pred_elaphrinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="weight (mg)", xlab="", main="Elaphrinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="ln(weight)", xlab="", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m2, newdata = pred_nebriinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Nebriinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m2, newdata = pred_scaritinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="", main="Scaritinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)

preds<- predict(m2, newdata = pred_trechinae, interval = "confidence")

plot(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, ylab="", xlab="body length (mm)", main="Trechinae", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=7)
polygon(c(rev(exp(pred_ln_size)), exp(pred_ln_size)), c(rev(exp(preds[ ,3])), exp(preds[ ,2])), col = 'mistyrose2', border = NA)
points(weights_schultz$weight_mg ~ weights_schultz$av_size_mm, pch=20, cex=0.8)
points(weights_booij$weight_mg ~ weights_booij$av_size_mm, pch=20, cex=0.8)
lines(exp(pred_ln_size), exp(preds[,1]), col="red", lw=2)

plot(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, ylab="", xlab="ln(body length)", main="", col="grey", cex=0.8,cex.lab=1.5, cex.axis=1.4, cex.main=1.4)
mtext("", 3, line=-2,cex=1.3, at=1.5)
polygon(c(rev(pred_ln_size), pred_ln_size), c(rev(preds[ ,3]), preds[ ,2]), col = 'mistyrose2', border = NA)
points(weights_schultz$ln_weight ~ weights_schultz$ln_size, pch=20, cex=0.8)
points(weights_booij$ln_weight ~ weights_booij$ln_size, pch=20, cex=0.8)
lines(pred_ln_size, preds[,1],col="red", lw=2)


# compile Table S1
pred<- predict(m0, newdata = weight_data, type="response")
weight_data$m0_pred <- exp(pred)
pred<- predict(m1, newdata = weight_data, type="response")
weight_data$m1_pred <- exp(pred)
pred<- predict(m1r, newdata = weight_data, re.form=NA)
weight_data$m1r_pred <- exp(pred)
pred<- predict(m2, newdata = weight_data, type="response")
weight_data$m2_pred <- exp(pred)

colnames(weight_data)
appendix_table <- weight_data[,c(6,1,13,12)]
appendix_table$dataset <- "training"
appendix_table$m_Szyszko <- round(szyszko_weight(weight_data$av_size_mm), digits = 2)
appendix_table$m_Booij <- round(booij_weight(weight_data$av_size_mm), digits = 2)
appendix_table$m_base <- round(weight_data$m0_pred, digits = 2)
appendix_table$m_fixed <- round(weight_data$m1_pred, digits = 2)
appendix_table$m_mixed <- round(weight_data$m1r_pred, digits = 2)
appendix_table$m_inter <- round(weight_data$m2_pred, digits = 2)

pred<- predict(m0, newdata = weight_data, type="response")
weight_data$m0_pred <- exp(pred)

pred<- predict(m1r, newdata = weight_data, re.form=NA)
weight_data$m1r_pred <- exp(pred)

colnames(weight_data)
appendix_table2 <- weight_data[,c(6,1,13,12)]
appendix_table2$dataset <- "training"
appendix_table2$m_Szyszko <- round(szyszko_weight(weight_data$av_size_mm), digits = 2)
appendix_table2$m_Booij <- round(booij_weight(weight_data$av_size_mm), digits = 2)
appendix_table2$m_base <- round(weight_data$m0_pred, digits = 2)
appendix_table2$m_fixed <- NA
appendix_table2$m_mixed <- round(weight_data$m1r_pred, digits = 2)
appendix_table2$m_inter <- NA

appendix_table <- rbind(appendix_table,appendix_table2)

appendix_table$dataset[appendix_table$source_weight == "schultz"]<- "validation"
exclude_subfamily <- c("Omophroninae", "Loricerinae", "Broscinae", "Cicindelinae")
appendix_table$dataset[appendix_table$subfamily %in% exclude_subfamily] <- "excluded"
colnames(species_double)
appendix_table$double[appendix_table$species %in% species_double & appendix_table$source_weight == "schultz"] <- "double"

appendix_table$source_weight[appendix_table$source_weight == "schultz"] <- "1)"
appendix_table$source_weight[appendix_table$source_weight == "booij"] <- "2)"
appendix_table$source_seize[appendix_table$source_seize == "booij"] <- "2)"
appendix_table$source_seize[appendix_table$source_seize == "freude"] <- "3)"

write_xlsx(appendix_table,"Table_S1")
