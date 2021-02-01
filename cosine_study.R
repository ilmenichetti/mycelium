w=(2*pi)/365
t=1:365
a=0.001

plot(t, 1+a*cos(w*t))



p=0.5
mu=0.08
t=seq(1,120)
t0=0
biomass= (p/mu)*(1-exp(-mu*(t-t0)))
  


calib_data$biom2

