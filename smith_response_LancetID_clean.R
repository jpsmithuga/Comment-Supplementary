######################################
# Author: Jonathan Smith, PhD, MPH
# Email: jonathan.p.smith@yale.edu
######################################
# Clear global environment and graphics, if desired
  # rm(list = ls())
  # dev.off()

#' - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' Global Parameters 
#' - - - - - - - - - - - - - - - - - - - - - - - - - - 
#' Default figures are black and white (FALSE), but
#' setting global option for color figures

color_figures <- FALSE # If color figures are preferred, set to TRUE
colors <- ggsci::pal_lancet(palette = c("lanonc"), alpha = 0.9)(9)

#' R values from:
#'   3.0 from Schneider and Eichner
#'   0.36 from Blumberg and Lloyd-Smith
R <- c(3.0, 0.36)

#' k values are arbitrary
#'   Note: k = 1000000 approximates
#'   infinity (i.e. Poisson distribution/
#'   uniform transmission)
k <- c(0.10, 0.50, 2.00, 1000000)

###############################
# Figure 1A
###############################

#' Function to calculate proportion of infection (see supplementary document)
#'   @param R Negative Binomial R value
#'   @param k Negative Binomial k value
#'   @param p Proportion of infectious cases responsible for transmission

propinfection <- function(R, k, prop) {
  xp <- qgamma(1 - prop, shape = k, rate = k / R)
  tp <- 1 - pgamma(xp, shape = k + 1, rate = k / R)
  return(tp)
}

# Set proportion driving secondary transmission, from 0 to 1 by 0.001
xx <- seq(0, 1, 0.001) 

#' Calculate proportion of expected transmission 
#' for a given proportion of infectious cases
plotdata <- sapply(seq_along(k), 
                   function(x) propinfection(R = R[1], k = k[x], xx))

# Plot Figure, Panel A
plot(range(0, 1), range(0, 1), type = "n", bty = "n",
     xlab = "Proportion of Infectious Cases",
     ylab = "Expected Proportion of Secondary Cases")
if(color_figures) {
  sapply(seq_along(k), 
         function(x) lines(xx, plotdata[, x], lwd = 2, col = colors[x], lty = 1))
  legend(x = 0.6, y = 0.35, c("k = 0.1", "k = 0.5", "k = 2.0", "No Variation"), 
         lty = 1, bty = "n", lwd = 2, col = colors[1:4])
  } else {
  sapply(seq_along(k), 
         function(x) lines(xx, plotdata[, x], lwd = 2, col = 1, lty = c(4:1)[x]))
  legend(x = 0.6, y = 0.35, c("k = 0.1", "k = 0.5", "k = 2.0", "No Variation"),
         lty = 4:1, bty = "n", lwd = 2)
}
# mtext(bquote(R[0] ~ "= 3.0"), adj = 1)
# mtext("A)", adj = -0.1)


###############################
# Figure 1B
###############################

#' Function to calculate the probability of an outbreak of size Y
#'   @param y Final outbreak size
#'   @param n Number of index cases (default to 1)
#'   @param R Negative Binomial R value
#'   @param k Negative Binomial k value
outbreakprob <- function(y, n = 1, R, k) {
  logp_yn <- log(n) - log(y) + 
    lgamma(k * y + y - n) - (lgamma(k * y) + lgamma(y - n + 1)) + 
    (y - n) * log(R / k) - (k * y + y - n) * log(1 + R / k)
  return(exp(logp_yn))
}

# Set figure parameters
maxoutbreak <- 100 # Plot to size 100
pchh <- c(2, 3, 8, 19) # point selection if if B&W

# Calculate outbreak probabilities
outbrk_prob <- sapply(seq_along(k), 
                      function(x) outbreakprob(y = 1:maxoutbreak, R = R[2], k = k[x]))

# Setup plot to max cluster size from y = 1 to max size by 1
maxrange <- seq(1, maxoutbreak, 1)

## Create plot
plot(log(range(1, maxoutbreak)), range(-11,0), 
     type = 'n', xlab = '', ylab = '', axes = FALSE)
axis(side = 1, at = log(c(1, 5, seq(10, maxoutbreak + 10, b = 10))), 
     label = c(1, 5, seq(10, maxoutbreak + 10, b = 10)))
axis(side = 2, at = log(c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1)), 
     label = c(1/100000, 1/10000, 1/1000, 1/100, 1/10, 1/2, 1))
mtext(side = 1, 'Final Outbreak Size', padj = 4)
mtext(side = 2, 'Probability', padj = -4)
if(color_figures){
  sapply(seq_along(k), 
         function(x) points(log(maxrange), log(outbrk_prob[,x]), 
                            col = colors[x], pch = 19, cex = 1))
  legend('topright', c("k = 0.1", "k = 0.5", "k = 2.0", "No Variation"), 
         pch = 19, col = colors[1:4], cex = 1, bty = "n")
} else {
  sapply(seq_along(k), function(x) points(log(maxrange), log(outbrk_prob[,x]), 
                                          col = 1, pch = pchh[x], cex = 1))
  legend('topright', c("k = 0.1", "k = 0.5", "k = 2.0", "No Variation"), 
         pch = c(2, 3, 8, 19), cex = 1, bty = "n")
}
#mtext(bquote(R[0] ~ "= 0.36"), adj = 1)
#mtext("B)", adj = -0.1)
