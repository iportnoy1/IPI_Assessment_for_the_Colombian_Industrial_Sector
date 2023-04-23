#Setting working directory
setwd("C:/Users/idpdl/Downloads")

#Calling built functions
source("Functions.R")

#Loading dataset
Data <- read.csv("Datos Depurados Bianca Collapsed.csv", header = T)

#Pre-treating dataset
rn <- Data[, 1]; Data <- as.data.frame(t(Data)); Data <- Data[-1,]; colnames(Data) <- rn

#Converting to numerical
for (i in 1:ncol(Data)) {Data[,i] = as.numeric(Data[,i])}

#Testing normality
pvals <- rep(0, ncol(Data))
for (i in 1:ncol(Data)) {pvals[i] <- (shapiro.test(Data[,i]))$p.value}

#PCA Testing
X_test = Data #Data[(which(rownames(Data) == "X2019.M10")):(nrow(Data)),]
X_train = Data[1:(which(rownames(Data) == "X2019.M10")),]
output <-  PCA_test(X_train, X_test)
T2 <- output$T2; T2a <- output$T2a; Q <- output$Q; Qa <- output$Qa; a <- output$a
x <- rep(0, nrow(Data)); Years <- 2005:2021
for (i in 1:17) {x[(12*(i-1)+1):(12*i)] <- Years[i]}

#plotting as time series
T2_ts <- ts(T2, start = c(2005,1),  end = c(2021,12), frequency = 12)
T2a_ts <- ts(T2a, start = c(2005,1),  end = c(2021,12), frequency = 12)
Q_ts <- ts(Q, start = c(2005,1),  end = c(2021,12), frequency = 12)
Qa_ts <- ts(Qa, start = c(2005,1),  end = c(2021,12), frequency = 12)

par(mfrow=c(1,2))
ts.plot(T2_ts, T2a_ts, xlab = "Period", ylab="T^2", gpars = list(col = c("Blue", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("topleft", legend=c("T^2", "Threshold"), col=c("blue", "red"), lty = 1:1, cex=0.8, 
       bty = "n")

ts.plot(Q_ts, Qa_ts, xlab = "Period", ylab="Q", gpars = list(col = c("Blue", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("topleft", legend=c("Q", "Threshold"), col=c("blue", "red"), lty = 1:1, cex=0.8, 
       bty = "n")

#Identifying most affected variables
Cont <- as.data.frame(contribution(Data, train_range = 1:178, test_range = 179:204))
Cumulative.Cont <- as.data.frame(rowSums(Cont))
colnames(Cumulative.Cont) <- "Cumulative.Contribution"
rownames(Cumulative.Cont)[which.max(Cumulative.Cont$Cumulative.Contribution)]
threshold = 70
important.vars<-as.data.frame(Cumulative.Cont[Cumulative.Cont$Cumulative.Contribution>=threshold,1])
rownames(important.vars)<-row.names(Cumulative.Cont)[Cumulative.Cont$Cumulative.Contribution>=threshold]
colnames(important.vars) <- "Contribution"

## Conveting variables to time series
Paper.Products <- ts(Data$`17 Paper and paper products`, start = c(2005,1), 
                     end = c(2021,12), frequency = 12)
Furniture <- ts(Data$`31 Furniture`, start = c(2005,1), end = c(2021,12), frequency = 12)
Basic.Metals <- ts(Data$`24 Basic metals`, start = c(2005,1), end = c(2021,12), frequency = 12)
Leather.Products <- ts(Data$`15 Leather and related products`, start = c(2005,1), end = c(2021,12), 
                       frequency = 12)
Electrical.Equipment <- ts(Data$`27 Electrical equipment`, start = c(2005,1), 
                           end = c(2021,12), frequency = 12)
Metal.Products <- ts(Data$`25 Fabricated metal products, except machinery`[1:178], 
                     start = c(2005,1), end = c(2021,12), frequency = 12)

Paper.Products.trunc <- ts(Data$`17 Paper and paper products`[1:178], start = c(2005,1), 
                           end = c(2019,10), frequency = 12)
Furniture.trunc <- ts(Data$`31 Furniture`[1:178], start = c(2005,1), end = c(2019,10), frequency = 12)
Basic.Metals.trunc <- ts(Data$`24 Basic metals`[1:178], 
                         start = c(2005,1), end = c(2019,10), frequency = 12)
Leather.Products.trunc <- ts(Data$`15 Leather and related products`[1:178], start = c(2005,1), 
                             end = c(2019,10), frequency = 12)
Electrical.Equipment.trunc <- ts(Data$`27 Electrical equipment`[1:178], start = c(2005,1), 
                                 end = c(2019,10), frequency = 12)
Metal.Products.trunc <- ts(Data$`25 Fabricated metal products, except machinery`[1:178], 
                           start = c(2005,1), end = c(2019,10), frequency = 12)

## Forecasting with Multiplicative Model with the Holt-Winters Method
 #Paper Products
Paper.Products.hw = HoltWinters(Paper.Products.trunc, seasonal="mult")
Paper.Products.predict = predict(Paper.Products.hw, n.ahead=26)

 #Furniture
Furniture.hw = HoltWinters(Furniture.trunc, seasonal="mult")
Furniture.predict = predict(Furniture.hw, n.ahead=26)

 #Basic Metals
Basic.Metals.hw = HoltWinters(Basic.Metals.trunc, seasonal="mult")
Basic.Metals.predict = predict(Basic.Metals.hw, n.ahead=26)

 #Leather Prducts
Leather.Products.hw = HoltWinters(Leather.Products.trunc, seasonal="mult")
Leather.Products.predict = predict(Leather.Products.hw, n.ahead=26)

 #Electrical Equipment
Electrical.Equipment.hw = HoltWinters(Electrical.Equipment.trunc, seasonal="mult")
Electrical.Equipment.predict = predict(Electrical.Equipment.hw, n.ahead=26)

 #Metal Products
Metal.Products.hw = HoltWinters(Metal.Products.trunc, seasonal="mult")
Metal.Products.predict = predict(Metal.Products.hw, n.ahead=26)

#Plotting all forecasts together
par(mfrow=c(3,2))
ts.plot(Paper.Products, Paper.Products.hw$fitted[,"xhat"], Paper.Products.predict, lty=c(1,2,2), 
        xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Paper Products", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")
#Furniture
ts.plot(Furniture, Furniture.hw$fitted[,"xhat"], Furniture.predict, lty=c(1,2,2), 
        xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Furniture", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")
#Basic Metals
ts.plot(Basic.Metals, Basic.Metals.hw$fitted[,"xhat"], Basic.Metals.predict, lty=c(1,2,2), 
        xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Basic Metals", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")
#Leather Prducts
ts.plot(Leather.Products, Leather.Products.hw$fitted[,"xhat"], Leather.Products.predict, 
        lty=c(1,2,2), xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Leather Products", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")
#Electrical Equipment
ts.plot(Electrical.Equipment, Electrical.Equipment.hw$fitted[,"xhat"], Electrical.Equipment.predict, 
        lty=c(1,2,2), xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Electrical Equipment", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")
#Metal Products
ts.plot(Metal.Products, Metal.Products.hw$fitted[,"xhat"], Metal.Products.predict, 
        lty=c(1,2,2), xlab = "Period", ylab = "Industrial Production Index", 
        gpars = list(col = c("Blue", "Red", "Red")))
abline(v = 2019.8333, col = "Magenta", lty = 2)
legend("bottomleft", legend=c("Metal Products", "Holt-Winters-Predicted"),
       col=c("blue", "red", "red"), lty = c(1,2,2), cex=0.8, bty = "n")