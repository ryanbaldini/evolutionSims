getFitnesses <- function(genePool, s)
{
	n <- ncol(genePool)/2
	fitness <- rep(NA, n)
	# chromosomeFitness <- colSums(genePool)
	# for(i in 1:n)
	# {
	# 	fitness[i] <- chromosomeFitness[2*i-1] + chromosomeFitness[2*i]
	# }
	for(i in 1:n)
	{
		fitness[i] <- sum(genePool[,2*i-1]*genePool[,2*i])
	}
	1-s*fitness
}
getLoad <- function(genePool, fitness)
{
	sum(genePool)*s/(ncol(genePool)/2)
}

n <- 10000
l <- 100
s <- 0.01
mu <- 1e-5

t <- 100

#Asexual
load <- rep(NA, t)
genePool <- matrix(rbinom(2*n*l, 1, (1e-3)/2), ncol=2*n, nrow=l)
newGenePool <- matrix(0, ncol=2*n, nrow=l)
count <- 1:l
for(i in 1:t)
{
	#reproduce
	fitness <- getFitnesses(genePool, s)
	load[i] <- 1-mean(fitness)
	reproducers <- sample(1:n, replace=TRUE, prob=fitness)
	for(j in 1:n)
	{
		copy <- as.numeric(genePool[,(reproducers[j]*2-1):(reproducers[j]*2)])
		newGenePool[,j*2-1] <- copy[count + l*rbinom(100,1,.5)]
		newGenePool[,j*2] <- copy[count + l*rbinom(100,1,.5)]
	}
	genePool <- newGenePool
	
	#mutate
	genePool <- (genePool + rbinom(n*2*l, 1, mu)) %% 2
}


#Sexual
T1 <- Sys.time()
load <- rep(NA, t)
genePool <- matrix(rbinom(2*n*l, 1, sqrt((1e-3))), ncol=2*n, nrow=l)
newGenePool <- matrix(0, ncol=2*n, nrow=l)
count <- 1:l
for(i in 1:t)
{
	#reproduce
	fitness <- getFitnesses(genePool, s)
	load[i] <- 1-mean(fitness)
	reproducers1 <- sample(1:n, replace=TRUE, prob=fitness)
	reproducers2 <- sample(1:n, replace=TRUE, prob=fitness)
	for(j in 1:n)
	{
		copy1 <- as.numeric(genePool[,(reproducers1[j]*2-1):(reproducers1[j]*2)])
		copy2 <- as.numeric(genePool[,(reproducers2[j]*2-1):(reproducers2[j]*2)])
		newGenePool[,j*2-1] <- copy1[count + l*rbinom(100,1,.5)]
		newGenePool[,j*2] <- copy2[count + l*rbinom(100,1,.5)]
	}
	genePool <- newGenePool
	
	#mutate
	genePool <- (genePool + rbinom(n*2*l, 1, mu)) %% 2
}
T2 <- Sys.time()
T2-T1