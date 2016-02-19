# author - Devashish Das 
# heavily borrowed things from 
# http://faculty.bscb.cornell.edu/~hooker/ShortCourseLab.R
#------------------------------

library(lubridate)
library(weatherData)
# obtain data from weatherData
year = as.character(seq(2001, 2015))
temp = matrix(0,365,15)
precip =matrix(0,365,15)
counter = 0
for(i in seq(2001, 2015)) {
	counter = counter + 1
	print(i)
	t0 = paste(year[counter],"-01-01",sep="")
	t1 = paste(year[counter],"-12-31",sep="")
	df = getSummarizedWeather("KRST", t0, t1,
		opt_custom_columns=TRUE,
		custom_columns=c(3,20))
	temp[,counter] = df[1:365,2]
	a = df[1:365,3]
	a[a == "T"] = "0.0"
	precip[, counter] = a
}
precip = (type.convert(precip))
total.precip = apply(precip, 2, sum)


#points where temp, rain in calculated
daytime = (1:365)-0.5

#scatter plot to see it
par(mfrow=c(1,2))
matplot(daytime,temp,type='p',pch=20, col="blue", 
	main = "Rochester, MN temperature (2001 -2015)",
	xlab = "day", ylab="temperature (F)")
plot(daytime,temp[,15],type='p',pch=20, col="blue", 
	main = "Rochester, MN temperature (2015)",
	xlab = "day", ylab="temperature (F)")

#install.packages("fda")
library("fda")

fbasis = create.fourier.basis(c(0,365),31)

# we can plot to see how it looks
plot(create.fourier.basis(c(0,365,5))

#a penalty term to penalize second derivative
curv.Lfd = int2Lfd(2)

#choose penalty
lambda = 1
curv.fdPar = fdPar(fbasis,curv.Lfd,lambda)
tempSmooth1 = smooth.basis(daytime,temp,curv.fdPar)

#see the contests of fd objects
names(tempSmooth1)

#lets plot them out
par(mfrow=c(1,2))
matplot(daytime,temp,type='p',pch=20, col=scales::alpha("blue",.25), 
	main = "Rochester, MN temperature (2001 -2015)",
	xlab = "day", ylab="temperature (F)")
lines(tempSmooth1$fd, lwd = 2, lty = 1, col = "black")

plot(daytime,temp[,15],type='p',pch=20, col=scales::alpha("blue",.25), 
	main = "Rochester, MN temperature (2001 -2015)",
	xlab = "day", ylab="temperature (F)")
lines(tempSmooth1$fd[15], lwd = 2, lty = 1, col = "black")


# very important to choose the smoothing parameters
lambdas = 10^seq(-4,4,by=0.5)    # lambdas to look over

mean.gcv = rep(0,length(lambdas)) # store mean gcv


for(ilam in 1:length(lambdas)){
  # Set lambda
  curv.fdPari = curv.fdPar
  curv.fdPari$lambda = lambdas[ilam]

  # Smooth
  tempSmoothi = smooth.basis(daytime,temp,curv.fdPari)

  # Record average gcv
  mean.gcv[ilam] = mean(tempSmoothi$gcv)
}

# We can plot what we have

plot(lambdas,mean.gcv,type='b',log='x')

# Lets select the lowest of these and smooth

best = which.min(mean.gcv)
lambdabest = lambdas[best]
curv.fdPar = fdPar(fbasis,curv.Lfd,lambdabest)
tempSmooth = smooth.basis(daytime,temp,curv.fdPar)

##################################################
## Functional PCA

tempfd = tempSmooth$fd
tempPCA = pca.fd(tempfd,nharm=15)

# Here we can look at proportion of variance explained:

plot(tempPCA$varprop,type='b')

## Looking at the principal components:

plot(tempPCA$harmonics[1:3])

plot(tempPCA,harm=1:9)



##------- FDA linear models

#### 1. Setup Data

# First we'll obtain log annual precipitation

annualprec = log10(total.precip)

# Now we need to set up a list of covariates.

xlist = list(len=2)

# First co-variate is just the intercept: a vector of ones

xlist[[1]] = rep(1,15)

# Second covariate is temperature

xlist[[2]] = tempfd


#### 2. fdPar objects for coeffients

# We also need a list of functional parameter objects to define the coefficient
# functions. 

bwtlist = list(len=2)

# First is a constant basis and add it to the fdPar object

cbasis = create.constant.basis(c(0,365))
bwtlist[[1]] = fdPar(cbasis)

# Now we need the coefficient of temperature, we'll use the same basis for it
# and add it as the second element in the list. 

beta.fdPar = curv.fdPar
bwtlist[[2]] = beta.fdPar



prec.model = fRegress(annualprec,xlist,bwtlist)

plot(prec.model$betaestlist[[2]])
