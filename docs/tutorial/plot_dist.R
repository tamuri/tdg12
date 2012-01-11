# A very small example of plotting the distributions in R
# Asif Tamuri (atamuri@nimr.mrc.ac.uk)

# Where are the distribution .csv files?
setwd('./results/exact')


# Read them in 
muts <- read.csv('distribution.mutations.csv', header=F)
subs <- read.csv('distribution.substitutions.csv', header=F)

# Proportion of deleterious, nearly neutral and advantageous mutations
proportions <- function(x, y) {
	result <- c(sum(y[x < -2]), sum(y[x >= -2 & x <= 2]), sum(y[x > 2])) 
	return(result)
}

# Adding proportions to the plots
plot_proportions <- function(x, y) {
	mtext(paste("p = (", paste(round(proportions(x, y), digits=6), collapse=", "), ")", collapse=""), side=3, line=0, adj=1.0, cex=1)
}


# Plot them!
par(mfrow=c(2,2))

# All mutations
plot(muts$V1, muts$V2, ty='h', lwd=2, main="Mutations", xlab="S", ylab="")
plot_proportions(muts$V1, muts$V2)

# Non-synonymous mutations, truncate the highly deleterious mutations
plot(muts$V1, muts$V3, ty='h', lwd=2, main="Non-synonymous mutations (truncated)", xlab="S", ylab="", 
	ylim=c(0, sort(muts$V3, dec=T)[2]))
plot_proportions(muts$V1, muts$V3)

# All substitutions
plot(subs$V1, subs$V2, ty='h', lwd=2, main="Substitutions", xlab="S", ylab="")
plot_proportions(subs$V1, subs$V2)

# Non-synonymous substitutions
plot(subs$V1, subs$V3, ty='h', lwd=2, main="Non-synonymous substitutions", xlab="S", ylab="")
plot_proportions(subs$V1, subs$V3)


