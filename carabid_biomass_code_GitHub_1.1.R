# Supplementary code to:

# How to estimate carabid biomass - An evaluation of different size-weight models for ground-beetles (Coleoptera:Carabidae) and a proposition for further improvement.

# submitted to Journal of Insect Conservation
# Fabio Weiss, Andreas Linde
# University for Sustainable Development Eberswalde, fabio.weiss@hnee.de

# compiled with R version 4.1.0. 


## Packagelist ####
library(lme4)
library(DHARMa)
library(MuMIn)
library(sjPlot)
library(sjmisc)
library(sjlabelled)
library(nlme)
library(writexl)

## Data preparation ####

data <- read.csv2("compiled_data.csv")

# data can be found here: 
# Schultz, R. (1996) Die Laufkäfer als Indikatoren der Renaturierung des Salzgrünlandes im Ostseebereich Vorpommerns. Culliver-Verlag, Göttingen.
# Booij, K., den Nijs, L., Heijermann, T., Jorritsma, I., Lock, C., Noorlander, J. (1994) Size and weight of carabid beetles: ecological applications. Proc Exper Appl Entomol.

data$ln_size <- log(weights_clean$av_size_mm)
data$ln_weight <- log(weights_clean$weight_mg)

single_subfamily <- c("Omophroninae", "Loricerinae", "Broscinae", "Cicindelinae" )

data <- data[!(data$subfamily %in% single_subfamily) ,]
data$subfamily <- as.factor(weights_clean2$subfamily)

data$size_section <- 1 
data$size_section[data$av_size_mm > max(data$av_size_mm)/3]   <- 2 
data$size_section[data$av_size_mm > max(data$av_size_mm)/3*2] <- 3 

max(data$av_size_mm[data$size_section==1])
max(data$av_size_mm)/3
max(data$av_size_mm)/3*2

data <- within(data, subfamily <- relevel(subfamily, ref = "Harpalinae"))

## model fitting ####

m0<- lm(ln_weight ~ ln_size, data= data)

m1<- lm(ln_weight ~ ln_size + subfamily, data= data)

m1r <-lmer(ln_weight ~ ln_size + (1|subfamily), data= data)

m2 <-lm(ln_weight ~ ln_size * subfamily, data= data)

tab_model(m0,m1,m1r,m2, dv.labels= c("m0","m1","m1r",  "m2"))

r.squaredGLMM(m1r)

## implementing existing size-weight models ####

# Szyszko, J. (1983) Methods of macrofauna investigations., in: Szyszko J (Eds) The Process of Forest Soil Macrofauna Formation after Afforestation Farmland. Warsaw Agricultural University Press, Warsaw, pp. 10-16.
# ln(weight) = -8.92804283 + 2.5554921 * ln(size)

syszko_weight<- function(size){
  ln_size <- log(size)
  ln_weight <- -8.92804283 + 2.5554921* ln_size
  weight_mg <- exp(ln_weight)*1000
  return(weight_mg)
}

# Booij, K., den Nijs, L., Heijermann, T., Jorritsma, I., Lock, C., Noorlander, J. (1994) Size and weight of carabid beetles: ecological applications. Proc Exper Appl Entomol.
# log(weight) = -1,3 + 2.95 * log(size)    

booij_weight<- function(size){
  log_size <- log10(size)
  log_weight <- (-1.3 + 2.95 *log_size)
  weight_mg <- 10^(log_weight)
  return(weight_mg)
}

## leave one out cross-validation ####

#m0
data$m0_loocv <- 0

its <- c(1:nrow(data))

for (i in its){
  loopdata <- data[-i,] 
  beetle <-  data[i,]
  m <-  lm(ln_weight ~ ln_size, data= loopdata)
  pred<- predict.lm(m, newdata = beetle, type="response")
  data[i,]$m0_loocv <- exp(pred)
}

# m1 #
data$m1_loocv <- 0

its <- c(1:nrow(data))

for (i in its){
  loopdata <- data[-i,] 
  beetle <-  data[i,]
  m <-  lm(ln_weight ~ ln_size + subfamily, data= loopdata)
  pred<- predict.lm(m, newdata = beetle, type="response")
  data[i,]$m1_loocv <- exp(pred)
}

# m1r#
data$m1r_loocv <- 0

its <- c(1:nrow(data))

for (i in its){
  loopdata <- data[-i,] 
  beetle <-  data[i,]
  m <-  lmer(ln_weight ~ ln_size + (1|subfamily), data= loopdata)
  pred<- predict(m, newdata = beetle, re.form=NA)
  data[i,]$m1r_loocv <- exp(pred)
}

# m2 #
data$m2_loocv <- 0

its <- c(1:nrow(data))

for (i in its){
  loopdata <- data[-i,] 
  beetle <-  data[i,]
  m <-  lm(ln_weight ~ ln_size * subfamily, data= loopdata)
  pred<- predict(m, newdata = beetle, type="response")
  data[i,]$m2_loocv <- exp(pred)
}


## predict with models of Szyszko and Booij et al. ####

data$syszko_weight <- syszko_weight(data$av_size_mm)

data$booij_weight <- booij_weight(data$av_size_mm)


## model evaluation ####

# calculate percent offset from real values 

data$m0_off <- (data$m0_loocv - data$weight_mg)/data$weight_mg *100

data$m1_off <- (data$m1_loocv - data$weight_mg)/data$weight_mg *100

data$m1r_off <- (data$m1r_loocv - data$weight_mg)/data$weight_mg *100

data$m2_off <- (data$m2_loocv - data$weight_mg)/data$weight_mg *100

data$syszko_off <- (data$syszko_weight - data$weight_mg)/data$weight_mg *100

data$booij_off <- (data$booij_weight - data$weight_mg)/data$weight_mg *100



# mean deviation per section 

section <- c(1:3)

f <- c(
  mean(abs(data$m0_off[data$size_section == 1])),
  mean(abs(data$m0_off[data$size_section == 2])),
  mean(abs(data$m0_off[data$size_section == 3]))
)

g <- c(
  mean(abs(data$m1_off[data$size_section == 1])),
  mean(abs(data$m1_off[data$size_section == 2])),
  mean(abs(data$m1_off[data$size_section == 3]))
)

z <- c(
  mean(abs(data$m1r_off[data$size_section == 1])),
  mean(abs(data$m1r_off[data$size_section == 2])),
  mean(abs(data$m1r_off[data$size_section == 3]))
)

h <- c(
  mean(abs(data$m2_off[data$size_section == 1])),
  mean(abs(data$m2_off[data$size_section == 2])),
  mean(abs(data$m2_off[data$size_section == 3]))
)

k<- c(
  mean(abs(data$syszko_off[data$size_section == 1])),
  mean(abs(data$syszko_off[data$size_section == 2])),
  mean(abs(data$syszko_off[data$size_section == 3]))
)

l <- c(
  mean(abs(data$booij_off[data$size_section == 1])),
  mean(abs(data$booij_off[data$size_section == 2])),
  mean(abs(data$booij_off[data$size_section == 3]))
)

model_performance <- data.frame(section=section , m0=f , m1=g , m2=h ,m1r=z , syszko=k , booij=l)

model_performance$section <- as.numeric(model_performance$section) 


# evaluation plots

# m0
par(mfrow= c(2,1))

plot(data$m0_off ~ data$av_size_mm, main="m0", ylim= c(-100,220), ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(data$av_size_mm)/3, max(data$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$m0, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,3, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$m0, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")

# m1
plot(data$m1_off ~ data$av_size_mm, main="m1", ylim= c(-100,220), ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(data$av_size_mm)/3, max(data$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$m1, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,4, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$m1, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")

#  m1r
plot(data$m1r_off ~ data$av_size_mm, main="m1r", ylim= c(-100,220), ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(data$av_size_mm)/3, max(data$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$m1r, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,4, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$m1r, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")

# m2 

plot(data$m2_off ~ data$av_size_mm, main="m2", ylim= c(-100,220), ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(weights_clean$av_size_mm)/3, max(weights_clean$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$m2, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,4, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$m2, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")

# szyszko

plot(data$syszko_off ~ data$av_size_mm, main="Model sensu Szyszko (1983)",ylim= c(-100,220), ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(data$av_size_mm)/3, max(data$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$syszko, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,4, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$syszko, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")

# Booij et al.

plot(data$booij_off ~ data$av_size_mm, main="Model sensu Booij et al (1994)",ylim= c(-100,220),  ylab= "percent of error", xlab="carabid body length (mm)",pch=20)
abline(h=0, col="red")
abline(v=c(max(data$av_size_mm)/3, max(data$av_size_mm)/3*2), lty="dashed")
plot( model_performance$section,model_performance$booij, type="l",lwd=2, ylim= c(0,60),xlim= c(0.65,3.5), xaxt="n", main= "Mean of absolute error per range section", xlab="", ylab="percent of error"  )
axis(1, at = seq(1,4, by = 1), las=1 )
text(model_performance$section,10, round(model_performance$booij, 2), cex=1.5, col="red")
abline(v=c(1.5,2.5), lty="dashed")


# overall mean error

models <- c("Szyszko", "Booij et al", "m0",  "m1","m1r", "m2")
models <- ordered(models, levels = c("Szyszko", "Booij et al", "m0",  "m1","m1r", "m2" ))

balanced_mean <- c(
  mean(model_performance$syszko),
  mean(model_performance$booij),
  mean(model_performance$m0),
  mean(model_performance$m1),
  mean(model_performance$m1r),
  mean(model_performance$m2))

balanced_mean <- data.frame(mean = balanced_mean, model =models)

par(mfrow=c(1,1))
plot(balanced_mean$mean,balanced_mean$models, ylim= c(20,35), xaxt="n", xlab="",main=" ", ylab="percent of error", pch=18,cex=2, col="black")
axis(1, at= c(1:6), labels = c("Szyszko", "Booij et al", "m0",  "m1","m1r", "m2"  ), las=2, cex.axis=1.1 )
abline(v=c(2.5), lty="dashed")
text(balanced_mean$model,balanced_mean$mean-2, round(balanced_mean$mean, 2), cex=1.2, col="black")


## catch simulation ####

par(mfrow=c(1,3))

# for smaller beetles 

random_catch_results <- data.frame(true=double(), szyszko=double(), booij=double(), m0=double() , m1=double(), m1r=double(), m2=double())

data2 <- data[data$size_section == 1,]

its = c(1:1000)
species = c(10:20)
individuals = c(80:120)

set.seed(128)
for (i in its){
  no_species <- sample(species,1)
  no_individuals <- sample(individuals,1)
  
  random_species <- data2[sample(nrow(data2), no_species, replace = F),]
  
  random_catch <- random_species[sample(nrow(random_species), no_individuals, replace = T),]
  
  random_catch_results[i,]$true <- sum(random_catch$weight_mg)
  random_catch_results[i,]$m0 <- sum(random_catch$m0_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1 <- sum(random_catch$m1_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1r <- sum(random_catch$m1r_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$m2 <- sum(random_catch$m2_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$szyszko <- sum(random_catch$syszko_weight)- sum(random_catch$weight_mg)
  random_catch_results[i,]$booij <- sum(random_catch$booij_weight)- sum(random_catch$weight_mg)
}

random_catch_results$szyszko_rel <- random_catch_results$szyszko/random_catch_results$true*100
random_catch_results$booij_rel <- random_catch_results$booij/random_catch_results$true*100
random_catch_results$m0_rel <- random_catch_results$m0/random_catch_results$true*100
random_catch_results$m1_rel <- random_catch_results$m1/random_catch_results$true*100
random_catch_results$m1r_rel <- random_catch_results$m1r/random_catch_results$true*100
random_catch_results$m2_rel <- random_catch_results$m2/random_catch_results$true*100

random_catch_synth_small <- data.frame(
  
  model=c("true", "szyszko", "booij", "m0", "m1", "m1r", "m2"), 
  
  mean=c(mean(random_catch_results$true),NA,NA,NA,NA,NA,NA), 
  sd=c(sd(random_catch_results$true),NA,NA,NA,NA,NA,NA),
  max=c(max(random_catch_results$true), NA,NA,NA,NA,NA,NA), 
  min=c(min(random_catch_results$true) ,NA,NA,NA,NA,NA,NA))

random_catch_synth_small$rel_err_mean <- c(NA, mean(abs(random_catch_results$szyszko_rel)), mean(abs(random_catch_results$booij_rel)),mean(abs(random_catch_results$m0_rel)),mean(abs(random_catch_results$m1_rel)),mean(abs     (random_catch_results$m1r_rel)),mean(abs(random_catch_results$m2_rel))) 

random_catch_synth_small$rel_err_median <- c(NA, median(random_catch_results$szyszko_rel), median(random_catch_results$booij_rel), median(random_catch_results$m0_rel),median(random_catch_results$m1_rel),median(random_catch_results$m1r_rel),median(random_catch_results$m2_rel))


plot(10, 1, xlim=c(1,6), ylim=c(-100,100), xaxt='n', xlab='', ylab="biomass error (%)", main="a", cex.main= 3,cex.lab= 2,cex.axis= 2)
axis(1, labels=colnames(random_catch_results[,2:7]), at=1:6,cex.axis= 2)

abline(h= 0, col="red")

for(i in 1:ncol(random_catch_results[,8:13])) {
  p <- random_catch_results[,i+7]
  p <- p[! p %in% 0]
  boxplot(p, add=T, at=i,cex.axis= 2 )
}


# for larger beetles

random_catch_results <- data.frame(true=double(), szyszko=double(), booij=double(), m0=double() , m1=double(), m1r=double(), m2=double())

data2 <- data[data$size_section != 1,]

its = c(1:1000)
species = c(10:20)
individuals = c(80:120)

set.seed(56)
for (i in its){
  no_species <- sample(species,1)
  no_individuals <- sample(individuals,1)
  
  random_species <- data2[sample(nrow(data2), no_species, replace = F),]
  
  random_catch <- random_species[sample(nrow(random_species), no_individuals, replace = T),]
  
  random_catch_results[i,]$true <- sum(random_catch$weight_mg)
  random_catch_results[i,]$m0 <- sum(random_catch$m0_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1 <- sum(random_catch$m1_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1r <- sum(random_catch$m1r_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$m2 <- sum(random_catch$m2_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$szyszko <- sum(random_catch$syszko_weight)- sum(random_catch$weight_mg)
  random_catch_results[i,]$booij <- sum(random_catch$booij_weight)- sum(random_catch$weight_mg)
}

random_catch_results$szyszko_rel <- random_catch_results$szyszko/random_catch_results$true*100
random_catch_results$booij_rel <- random_catch_results$booij/random_catch_results$true*100
random_catch_results$m0_rel <- random_catch_results$m0/random_catch_results$true*100
random_catch_results$m1_rel <- random_catch_results$m1/random_catch_results$true*100
random_catch_results$m1r_rel <- random_catch_results$m1r/random_catch_results$true*100
random_catch_results$m2_rel <- random_catch_results$m2/random_catch_results$true*100

random_catch_synth_large <- data.frame(
  
  model=c("true", "szyszko", "booij", "m0", "m1", "m1r", "m2"), 
  
  mean=c(mean(random_catch_results$true),NA,NA,NA,NA,NA,NA), 
  sd=c(sd(random_catch_results$true),NA,NA,NA,NA,NA,NA),
  max=c(max(random_catch_results$true), NA,NA,NA,NA,NA,NA), 
  min=c(min(random_catch_results$true) ,NA,NA,NA,NA,NA,NA))

random_catch_synth_large$rel_err_mean <- c(NA, mean(abs(random_catch_results$szyszko_rel)), mean(abs(random_catch_results$booij_rel)),mean(abs(random_catch_results$m0_rel)),mean(abs(random_catch_results$m1_rel)),mean(abs     (random_catch_results$m1r_rel)),mean(abs(random_catch_results$m2_rel))) 

random_catch_synth_large$rel_err_median <- c(NA, median(random_catch_results$szyszko_rel), median(random_catch_results$booij_rel), median(random_catch_results$m0_rel),median(random_catch_results$m1_rel),median(random_catch_results$m1r_rel),median(random_catch_results$m2_rel))

plot(10, 1, xlim=c(1,6), ylim=c(-100,100), xaxt='n', xlab='', ylab="", main="b", cex.main= 3,cex.lab= 1.7,cex.axis= 2)
axis(1, labels=colnames(random_catch_results[,2:7]), at=1:6, cex.main= 2,cex.lab= 2,cex.axis= 2)

abline(h= 0, col="red")

for(i in 1:ncol(random_catch_results[,8:13])) {
  p <- random_catch_results[,i+7]
  p <- p[! p %in% 0]
  boxplot(p, add=T, at=i, cex.axis=2)
}


# mixed scenario 


random_catch_results <- data.frame(true=double(), szyszko=double(), booij=double(), m0=double() , m1=double(), m1r=double(), m2=double())

weights_clean_large <- data[data$size_section != 1,]
weights_clean_small <- data[data$size_section == 1,]

its = c(1:1000)
species1 = c(5:10)
species2 = c(5:10)
individuals = c(80:120)

set.seed(12)
for (i in its){
  no_species1 <- sample(species1,1)
  no_species2 <- sample(species2,1)
  
  no_individuals <- sample(individuals,1)
  
  random_species1 <- weights_clean_small[sample(nrow(data2),  no_species1, replace = F),]
  random_species2 <- weights_clean_large[sample(nrow(data2),  no_species2, replace = F),]
  random_species <- rbind(random_species1,random_species2)
  
  random_catch <- random_species[sample(nrow(random_species), no_individuals, replace = T),]
  
  random_catch_results[i,]$true <- sum(random_catch$weight_mg)
  random_catch_results[i,]$m0 <- sum(random_catch$m0_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1 <- sum(random_catch$m1_loocv) - sum(random_catch$weight_mg)
  random_catch_results[i,]$m1r <- sum(random_catch$m1r_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$m2 <- sum(random_catch$m2_loocv)- sum(random_catch$weight_mg)
  random_catch_results[i,]$szyszko <- sum(random_catch$syszko_weight)- sum(random_catch$weight_mg)
  random_catch_results[i,]$booij <- sum(random_catch$booij_weight)- sum(random_catch$weight_mg)
}

random_catch_results$szyszko_rel <- random_catch_results$szyszko/random_catch_results$true*100
random_catch_results$booij_rel <- random_catch_results$booij/random_catch_results$true*100
random_catch_results$m0_rel <- random_catch_results$m0/random_catch_results$true*100
random_catch_results$m1_rel <- random_catch_results$m1/random_catch_results$true*100
random_catch_results$m1r_rel <- random_catch_results$m1r/random_catch_results$true*100
random_catch_results$m2_rel <- random_catch_results$m2/random_catch_results$true*100


random_catch_synth_mixed <- data.frame(
  
  model=c("true", "szyszko", "booij", "m0", "m1", "m1r", "m2"), 
  
  mean=c(mean(random_catch_results$true),NA,NA,NA,NA,NA,NA), 
  sd=c(sd(random_catch_results$true),NA,NA,NA,NA,NA,NA),
  max=c(max(random_catch_results$true), NA,NA,NA,NA,NA,NA), 
  min=c(min(random_catch_results$true) ,NA,NA,NA,NA,NA,NA))

random_catch_synth_mixed$rel_err_mean <- c(NA, mean(abs(random_catch_results$szyszko_rel)), mean(abs(random_catch_results$booij_rel)),mean(abs(random_catch_results$m0_rel)),mean(abs(random_catch_results$m1_rel)),mean(abs(random_catch_results$m1r_rel)),mean(abs(random_catch_results$m2_rel))) 

random_catch_synth_mixed$rel_err_median <- c(NA, median(random_catch_results$szyszko_rel), median(random_catch_results$booij_rel), median(random_catch_results$m0_rel),median(random_catch_results$m1_rel),median(random_catch_results$m1r_rel),median(random_catch_results$m2_rel))


plot(10, 1, xlim=c(1,6), ylim=c(-100,100), xaxt='n', xlab='', ylab="", main="c", cex.main= 3,cex.lab= 1.7,cex.axis= 2)
axis(1, labels=colnames(random_catch_results[,2:7]), at=1:6, cex.axis=2)

abline(h= 0, col="red")

for(i in 1:ncol(random_catch_results[,8:13])) {
  p <- random_catch_results[,i+7]
  p <- p[! p %in% 0]
  boxplot(p, add=T, at=i, cex.axis=2)
}


## Further plots and tables ####

# data plot (introduction)

par(mfrow= c(1,2))
plot(weights_clean2$weight_mg ~ weights_clean2$av_size_mm, pch=8, main= "a", ylab="weight (mg)", xlab="body length (mm)", col="grey")

pred_size <- seq(min(weights_clean$av_size_mm), max(weights_clean$av_size_mm), length.out = 100)

lines(pred_size, syszko_weight(pred_size), lty="dotted", lwd=2)
lines(pred_size, booij_weight(pred_size), lty="twodash", lwd=2)

#on the log scale

plot(weights_clean$ln_weight ~ weights_clean$ln_size, pch=8, main= "b", ylab="ln(weight)", xlab="ln(body length)",col="grey")

pred_ln_size <- log(pred_size)

lines(pred_ln_size, log(syszko_weight(pred_size)),lty="dotted", lw=2)
lines(pred_ln_size, log(booij_weight(pred_size)),lty="twodash", lw=2)




# data table for supporting information

# predict beetle weights with different models (without CV)
pred<- predict(m0, newdata = weights_clean2, type="response")
weights_clean2$m0_pred <- exp(pred)

pred<- predict(m1, newdata = weights_clean2, type="response")
weights_clean2$m1_pred <- exp(pred)

pred<- predict(m1r, newdata = weights_clean2, re.form=NA)
weights_clean2$m1r_pred <- exp(pred)

pred<- predict(m2, newdata = weights_clean2, type="response")
weights_clean2$m2_pred <- exp(pred)


appendix_table <- weights_clean2[,c(6,1,13,12,22,23,24:27)]
appendix_table$syszko_weight <- round(appendix_table$syszko_weight, digits = 2)
appendix_table$booij_weight <- round(appendix_table$booij_weight, digits = 2)
appendix_table$m0_pred <- round(appendix_table$m0_pred, digits = 2)
appendix_table$m1_pred <- round(appendix_table$m1_pred, digits = 2)
appendix_table$m1r_pred <- round(appendix_table$m1r_pred, digits = 2)
appendix_table$m2_pred <- round(appendix_table$m2_pred, digits = 2)

appendix_table$source_weight[appendix_table$source_weight == "schultz"] <- "1)"
appendix_table$source_weight[appendix_table$source_weight == "booij"] <- "2)"
appendix_table$source_seize[appendix_table$source_seize == "booij"] <- "2)"
appendix_table$source_seize[appendix_table$source_seize == "freude"] <- "3)"

write_xlsx(appendix_table,"C:/Users/fweiss/Promotion/PhD - own research/projects/Biomasse/data/appendix_data.xlsx")





