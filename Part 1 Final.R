###Author - Nathan Wycoff


require(spBayes)
require(deSolve)

#Aphid data and model from https://cran.r-project.org/web/packages/deSolve/vignettes/deSolve.pdf
Aphid <- function(t, APHIDS, parameters) {
  deltax     <- c (0.5, rep(1, numboxes - 1), 0.5)
  Flux       <- -D * diff(c(0, APHIDS, 0)) / deltax
  dAPHIDS    <- -diff(Flux) / delx + APHIDS * r
  list(dAPHIDS )
} # end

D         <- 1    # m2/day  diffusion rate
r         <- 0.01   # /day    net growth rate
delx      <- 1      # m       thickness of boxes
numboxes  <- 60
# distance of boxes on plant, m, 1 m intervals
Distance  <- seq(from = 0.5, by = delx, length.out = numboxes)

# Initial conditions:  # ind/m2
APHIDS        <- rep(0, times = numboxes)
APHIDS[30:31] <- 100
state         <- c(APHIDS = APHIDS)      # initialise state variables

times <-seq(0, 300, by = 1)
print(system.time(
  out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")
  ))

head(out[,1:5])

image(out, method = "filled.contour", grid = Distance,
      xlab = "time, days", ylab = "Distance on plant, m",
      main = "Aphid density on a row of plants")