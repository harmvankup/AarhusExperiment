# Set the length and time parameters for the simulation
length <- 100 # length of the porous medium in millimeters
time <- 4 # time of the simulation in days
bulkconcentration <- 400 #in nmol/mL or umol/L 

# Set the diffusion coefficient and the porosity of the medium
D <- 5.27*10^(-6) # diffusion coefficient in cm^2/sec
porosity <- 0.9 # porosity of the medium
a <- 0.8
m <- 2.1
ffactor <- a*porosity^(-m) # Formation factor



# Set the time step and the number of iterations
dt <- 1 # time step in minutes
dx <- 0.1 # step in cm
iterations <- 24*60*time/dt # number of iterations

# Set the reaction rate at the left boundary
k <- 0.1 # reaction rate in mol/min^-1

# Set the initial and boundary conditions
concentration <- rep(bulkconcentration, length) # initial concentration of the substance is 1 at all points
concentration[1] <- 0 # concentration of the substance at the left boundary is 0
# Initialize the time and concentration arrays
times <- c() # time array
concentrations <- list() # concentration array
fluxes <- c() # the flux in nmol cm-2 s-1
accumulation <- c() # in nmol
theor_concentration <- c()
# Loop through the iterations and update the concentration at each point
for (i in 1:iterations){
  times <- c(times, (i-1)*(dt/60)) # add the current time to the time array in h
  concentrations[[i]] <- concentration # add the current concentration to the concentration array
  fluxes[[i]] <- concentration[2]*D/dx # in nmol cm-2 s-1
  accumulation[[1]] <- 0
  accumulation[[i+1]] <-accumulation[[i]] + fluxes[[i]]*dt*60*1.8*0.5
  theor_concentration[[i]] <-  accumulation[[i]]*0.094/(D*times[[i]]*60*60*1.8*0.5)
  #concentration[1] <- #concentration[1] + dt*k + D*dt*(concentration[2] - concentration[1])/((1/1000)^2*ffactor*porosity) # update concentration at the left boundary
  for (j in 2:(length-1)){
    concentration[j] <- concentration[j] + D*(dt*60)*(concentration[j+1] - 2*concentration[j] + concentration[j-1])/((dx)^2*ffactor*porosity)
  }
}
concentration[1] + D*dt*(concentration[j+1] - 2*concentration[j] + concentration[j-1])/((1/1000)^2*ffactor*porosity)

# Plot the concentration profiles at different time steps in the same figure
plot(concentrations[[2]], type="l", xlab="Position (mm)", ylab="Concentration", col="red", main="Concentration Profile")
legend_text <- c() # initialize the legend text array
for (i in seq(from=1, to=length(times), by=((iterations)/10))){
  lines(concentrations[[i]], col=rgb(0,0,i/length(times)))
  legend_text <- c(legend_text, paste(times[i], "h")) # add the time label to the legend text array
}
#legend("topright", legend=legend_text, col=rgb(0,0,seq(from=1, to=length(times), by=(iterations)/10)/length(times)), lty=1)

plot(times,fluxes)
plot(times,accumulation[-1])
plot(times,theor_concentration)

# Add the legend


testlength = seq(0.01,10,by = 0.01)
#test <-  testlength + k*dt*testlength + 
test <-   D*dt*(0 - testlength)/((1/1000)^2*ffactor)

plot(test, testlength, type="l")
