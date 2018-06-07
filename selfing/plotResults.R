#Plot results

#Read data
data <- read.csv("output.csv")

#Calculate mean on second half
meanLoad <- mean(data$load[-(1:(nrow(data)/2))])

#Save plot
jpeg('loadPlot.jpg')
plot(data$load, type="l", lwd=2, xlab = "generation", ylab="load", main = paste("Mean = ", meanLoad))
abline(h = meanLoad, lty=2)
dev.off()

#Calculate mean on second half
meanSelfingFrequency <- mean(data$selfingLocusFrequency[-(1:(nrow(data)/2))])

#Save plot
jpeg('selfingPlot.jpg')
plot(data$selfingLocusFrequency, type="l", lwd=2, xlab = "generation", ylab="selfing frequency", main = paste("Mean = ", meanSelfingFrequency))
abline(h = meanSelfingFrequency, lty=2)
dev.off()