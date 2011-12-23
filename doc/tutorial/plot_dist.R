# A very small example of plotting the distributions in R
# Asif Tamuri (atamuri@nimr.mrc.ac.uk)

# Where are the distribution .csv files?
setwd('./results/exact')


# Read them in 
muts <- read.csv('distribution.mutations.csv', header=F)
subs <- read.csv('distribution.substitutions.csv', header=F)


# Plot them!
par(mfrow=c(2,2))

# All mutations
plot(muts$V1, muts$V2, ty='h', lwd=2, main="Mutations", xlab="S", ylab="")

# Non-synonymous mutations, truncate the highly deleterious mutations
plot(muts$V1, muts$V3, ty='h', lwd=2, main="Non-synonymous mutations\n(truncated)", xlab="S", ylab="", 
	ylim=c(0, sort(muts$V3, dec=T)[2]))

# All substitutions
plot(subs$V1, subs$V2, ty='h', lwd=2, main="Substitutions", xlab="S", ylab="")


# Non-synonymous substitutions
plot(subs$V1, subs$V3, ty='h', lwd=2, main="Non-synonymous substitutions", xlab="S", ylab="")






