# Load required libraries
library(MASS)
library(ggplot2)

# set sample size and no. of snps
n <- 10000   
p <- 100000  

# initialise seed and random-number generator
S <- 192398123
set.seed(S)

############ Genotype simulation:
# draw minor allele frequencies from a uniform distribution constrain to [0.05,0.5] interval
MAF <- runif(1, 0.05, 0.5)

# draw genotypes
# Function to generate genotypes for a single SNP
binom <- function(x, MAF) {
  rbinom(x, 2, MAF)  # Generates genotypes (0, 1, 2) for 'x' individuals
}

# Generate genotype matrix
geno <- replicate(p, binom(n, MAF))  # Generate p columns, each for n individuals


# Calculate observed MAF
calculate_maf <- function(genotypes) {
  # Count the minor allele (sum of all alleles divided by total alleles)
  allele_freq <- sum(genotypes) / (2 * length(genotypes))
  # MAF is the minimum of the allele frequency or 1 - allele frequency
  min(allele_freq, 1 - allele_freq)
}

# Calculate MAF for each SNP
(observed_mafs <- apply(geno, 2, calculate_maf))


######### Phenotype simulation:
# Simulate parameters first:
alpha1 = rnorm(p)
alpha2 = rnorm(p)
beta1 = rnorm(p)
beta2 = rnorm(p)
gamma = rnorm(p)

# Function to calculate mu:
calc_mu <- function(x, alpha1, alpha2) {
  c(mu1 = sum(x * alpha1), mu2 = sum(x * alpha2))
}

# Generate mu and scale:
mu <- t(sapply(1:n, function(i) {
  mu <- calc_mu(geno[i,], alpha1, alpha2) #mu(xi)
  print(i)
  scale(mu)
}))

# Function to calculate V
calc_V <- function(x, beta1, beta2, gamma) {
  # Measuring sigma1 and sigma2:
  sigma1 <- exp(sum(x * beta1) / p) / (1 + exp(sum(x * beta1) / p)) #make use of sigmoid function to bring sum(x * beta1) to the [0,1] domain
  sigma2 <- exp(sum(x * beta2) / p) / (1 + exp(sum(x * beta2) / p)) #divide by p to avoid te variances to explode
  # Measuring rho:
  numerator <- exp(sum(x * gamma)) - 1
  denominator <- exp(sum(x * gamma)) + 1
  # Handle cases of Inf or -Inf explicitly
  if (is.infinite(numerator) && numerator > 0) {
    rho=1                            # Assign 1 if numerator is +Inf
  } else if (is.infinite(numerator) && numerator < 0) {
    rho=-1                           # Assign -1 if numerator is -Inf
  } else {
    rho = (numerator / denominator)  # Regular calculation
  }  
  # Returning the variance-covariance matrix:
  matrix(c(sigma1^2, sigma1 * sigma2 * rho,
           sigma1 * sigma2 * rho, sigma2^2),
         nrow = 2)
}


# Generate phenotypes scaled
phenotypes_scaled <- t(sapply(1:n, function(i) {
  V <- calc_V(geno[i,], beta1, beta2, gamma) # V(xi)
  print(i)
  MASS::mvrnorm(1, mu[i,], V)
}))

# Results
colnames(phenotypes_scaled) <- c("y1", "y2")
head(phenotypes_scaled)

apply(phenotypes_scaled, 2, mean)
apply(phenotypes_scaled, 2, var) # variance is very high, how to control for that?

par(mfrow=c(1,1))
hist(phenotypes_scaled[,1], prob=TRUE)
hist(phenotypes_scaled[,2], add=TRUE, col=2, prob=TRUE)

cor(phenotypes_scaled) 

# Generate phenotypes (not scaling mu)
phenotypes <- t(sapply(1:n, function(i) {
  mu <- calc_mu(geno[i,], alpha1, alpha2) #mu(xi)
  V <- calc_V(geno[i,], beta1, beta2, gamma) # V(xi)
  print(V)
  MASS::mvrnorm(1, mu, V)
}))

# Results
colnames(phenotypes) <- c("y1", "y2")
head(phenotypes)

apply(phenotypes, 2, mean)
apply(phenotypes, 2, var) # variance is very high, how to control for that?

par(mfrow=c(1,1))
hist(phenotypes[,1], prob=TRUE)
hist(phenotypes[,2], add=TRUE, col=2, prob=TRUE)

cor(phenotypes) 
# CORRELATION: expected to be high, try to include in your simulation
#              study a way to account for sources of simulation
#              environment, genotype correlation, errors


# Simulate error terms: bivariate, according to the model
e <- rt(n, 6) # Chisq distribution with df = 6  (normqality for e to assume as a baseline)
e <- (e - 6)/sqrt(2*6) # standardized to have a mean of 0 and variance of 1
# Simulate Environment: according to the model
E <- rnorm(n)
i =0


# Simulate the effect
#  alpha(1/2) is the mean effect, 
#  beta_(1/2) is the variance effect, 
#  gamma is the covariance effect
###
# Draw betas, scaling, and genetic components
# If I want to introduce genetic correlation and 
# environmental covariance
# GENETIC CORRELATION: affects how alpha1 and alpha2 are related
# would expect that wiincreasing genetic correlation
# the number of significant SNPs relevant to the phenotypical variance
# would increase
# GENOTYPICAL CORRELATION: introduces shared env influence
rhoG <- 0.1
# ENVIRONMENTAL CORRELATION: introduces shared env influence
rhoE <- 0.3

# Set SNP-based heritability of y1 and y2: proportion of variance
# eplained by the genetic component
h2_1 <- 0.2
h2_2 <- 0.1

# Simulation of the MEAN EFFECTS:
alpha1_g = rnorm(p)
alpha2_g <- rhoG * alpha1 + sqrt(1 - rhoG^2) * rnorm(p) ## multiplies by (1 minus rhoG) so that var alpha2_g is 1
#s <- sqrt(2 * MAF * (1 - MAF))^(-1)         #1/√genotypic variance 
#y1_g <- rowSums(geno * (s * alpha1)) / sqrt(p)  #standardized 
#y2_g <- rowSums(geno * (s * alpha2)) / sqrt(p)
# Draw environmental components:
y1_e <- rnorm(n)
y2_e <- rhoE * y1_e + sqrt(1 - rhoE^2) * rnorm(n) #to simulate correlated environments

# Draw outcomes:
phenotypes_c <- t(sapply(1:n, function(i) {
  print(i)
  mu_g <- calc_mu(geno[i,], alpha1_g, alpha2_g) #1x2 vector
  V <- calc_V(geno[i,], beta1, beta2, gamma)
  Y <- MASS::mvrnorm(1, mu_g, V)
  y1 <- sqrt(h2_1) * Y[1] + sqrt(1 - h2_1) * y1_e[i]
  y2 <- sqrt(h2_2) * Y[2] + sqrt(1 - h2_2) * y2_e[i]
  return(c(y1,y2))
}))

# Results
colnames(phenotypes_c) <- c("y1", "y2")
apply(phenotypes_c, 2, mean)
apply(phenotypes_c, 2, var) # variance is very high, how to control for that?

par(mfrow=c(1,1))
hist(phenotypes_c[,1], prob=TRUE)
hist(phenotypes_c[,2], add=TRUE, col=2, prob=TRUE)

cor(phenotypes_c) 







##  Variance explianed by mean effects and variance effects
#rsq_g <- 0.005
#rsq_gxe <- 0.005
#rsq_e <- 1 - rsq_g * ifelse(alpha_X^2 >0, 1, 0) - rsq_gxe * ifelse(beta_X^2 >0, 1, 0)

# Simulate outcome
#g = geno * alpha_g
#gxe = geno * E * beta_g
#y = numscale(g)*sqrt(rsq_g) + numscale(gxe)*sqrt(rsq_gxe) + numscale(e)*sqrt(rsq_e)



# get centered variables and squared sums
#GC <- scale(G, scale = FALSE)
#SGC <- colSums(GC^2)
#XC <- X - mean(X)
#YC <- Y - mean(Y)

# get GWAS betas
#bX <- colSums(GC * XC) / SGC
#bY <- colSums(GC * YC) / SGC

# get GWAS residuals per SNP
#resX <- XC - sweep(GC, 2, bX, FUN = "*")
#resY <- YC - sweep(GC, 2, bX, FUN = "*")

# get GWAS est of sig2 per SNP
#sig2X <- colSums(resX^2) / (n - 2)
#sig2Y <- colSums(resY^2) / (n - 2)

# get GWAS standard errors
#SEbX <- sqrt(sig2X / SGC)
#SEbY <- sqrt(sig2Y / SGC)

# get GWAS z stats
#zX <- bX / SEbX
#zY <- bY / SEbY

# get p values
#pX <- 2 * pnorm(-abs(zX))
#pY <- 2 * pnorm(-abs(zY))

# get signed -10log(pval)
#ignedM10logPX <- sign(bX) * (-log10(pX))
#signedM10logPY <- sign(bY) * (-log10(pY))

# get most extreme p value
#tau <- ceiling(max(c(max(abs(signedM10logPX)), max(abs(signedM10logPY)))))

#GW_pval_threshold= c(-1, 1)*-log10(5*10^(-8))
#nominally_pval_threshold = c(-1,1)*-log10(5*10^(-6))

# create scatter plot
#df <- data.frame(signedM10logPX, signedM10logPY)
#gplot(df, aes(x = signedM10logPX, y = signedM10logPY, color = sqrt(signedM10logPX^2 + signedM10logPY^2))) +
#  geom_point(size = 0.5) +
#  scale_color_gradient(low = "blue", high = "red") +
#  xlim(-tau, tau) +
#  ylim(-tau, tau) +
#  xlab('Signed minus log(p-value) for X') +
#  ylab('Signed minus log(p-value) for Y') +
#  ggtitle(paste('LA plot\nrhoG =', rhoG, ', rhoE =', rhoE)) +
#  theme(aspect.ratio = 1) +
#  labs(color = "Magnitude") + 
# Euclidean distance from the origin (representing the magnitude of the 
# combined signed log p-values) 



