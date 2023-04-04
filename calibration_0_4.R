
###loading the packages needed
library(rstan) #this is to call STAN from R
library(ggsci) #for color palettes
library(hydroGOF) #for RMSE
library(modeest) #
library(plyr) #
library(lubridate) #
library(shape) #


###defining some custom functions

# Add an alpha value to a colour for transparencies
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2,
        function(x)
          rgb(x[1], x[2], x[3], alpha=alpha))
}

#max and minimum in a table
colMax <- function (colData) {
  apply(colData, MARGIN=c(2), max)
}
colMin <- function (colData) {
  apply(colData, MARGIN=c(2), min)
}

#sd for a matrix by column
colSdApply <- function(x, ...)apply(X=x, MARGIN=2, FUN=sd, ...)



###creating some color palettes to get nice effects
chain_palette<-pal_uchicago("light", 1)(7)

parms_palette<-pal_npg("nrc", 1)(9)
parms_palette_alpha<-pal_npg("nrc", 0.5)(9)

times_palette<-pal_startrek("uniform", 1)(6)
times_palette_alpha<-pal_startrek("uniform", 0.5)(6)



###reading the data
dataset<-read.csv(file = "Incubation_2015.csv")
dataset$Biomass.g.m2<-as.numeric(as.character(dataset$Biomass.g.m2)) # in case the column is read as a factor
dataset$Incubation.time.days<-dataset$Incubation.time..months.*30 #to get the days
dataset$classes<-as.factor(dataset$classes)

dataset$Start<-as.Date(dataset$Start,"%d/%m/%y")
dataset$End<-as.Date(dataset$End,"%d/%m/%y")

dataset$Start_day<-yday(dataset$Start)
dataset$End_day<-yday(dataset$End)

#create a vector with the periods
min(dataset$Start_day)
max(dataset$Start_day)
periods=c(min(dataset$Start_day), min(dataset$Start_day)+floor((304-181)/4), min(dataset$Start_day)+floor((304-181)/4)*2, min(dataset$Start_day)+floor((304-181)/4)*3)

periods_vec<-dataset$Start_day
periods_vec<-NA
periods_vec[(dataset$Start_day>=periods[1] & dataset$Start_day<(periods[2]+1))]=1
periods_vec[(dataset$Start_day>=periods[2] & dataset$Start_day<(periods[3]+1))]=2
periods_vec[(dataset$Start_day>=periods[3] & dataset$Start_day<(periods[4]+1))]=3
periods_vec[(dataset$Start_day>=periods[4])]=4


#fill the zeroes
which_ones<-which(dataset$Biomass.g.m2==0)
for(i in 1:length(which_ones)){
  class<-dataset[which_ones[i],]$classes
  time_inc<-dataset[which_ones[i],]$Incubation.time.days
  treat<-dataset[which_ones[i],]$Treatment.name
  dataset$Biomass.g.m2[which_ones[i]]<-mean(dataset[dataset$classes==class & dataset$Incubation.time.days==time_inc & dataset$Treatment.name==treat,]$Biomass.g.m2)
  }


###selecting the four blocks as separate inputs into the Bayesian model
biom1<-(dataset[dataset$classes==levels(dataset$classes)[1],]$Biomass.g.m2)
start1<-dataset[dataset$classes==levels(dataset$classes)[1],]$Start_day-182
time1<-dataset[dataset$classes==levels(dataset$classes)[1],]$End_day-182
tree1<-dataset[dataset$classes==levels(dataset$classes)[1],]$Tree_vol_increment.m3.ha.Y.
periods_vec1<-periods_vec[dataset$classes==levels(dataset$classes)[1]]

biom2<-(dataset[dataset$classes==levels(dataset$classes)[2],]$Biomass.g.m2)
start2<-dataset[dataset$classes==levels(dataset$classes)[2],]$Start_day-182
time2<-dataset[dataset$classes==levels(dataset$classes)[2],]$End_day-182
tree2<-dataset[dataset$classes==levels(dataset$classes)[2],]$Tree_vol_increment.m3.ha.Y.
periods_vec2<-periods_vec[dataset$classes==levels(dataset$classes)[2]]

biom3<-(dataset[dataset$classes==levels(dataset$classes)[3],]$Biomass.g.m2)
start3<-dataset[dataset$classes==levels(dataset$classes)[3],]$Start_day-182
time3<-dataset[dataset$classes==levels(dataset$classes)[3],]$End_day-182
tree3<-dataset[dataset$classes==levels(dataset$classes)[3],]$Tree_vol_increment.m3.ha.Y.
periods_vec3<-periods_vec[dataset$classes==levels(dataset$classes)[3]]

biom4<-(dataset[dataset$classes==levels(dataset$classes)[4],]$Biomass.g.m2)
start4<-dataset[dataset$classes==levels(dataset$classes)[4],]$Start_day-182
time4<-dataset[dataset$classes==levels(dataset$classes)[4],]$End_day-182
tree4<-dataset[dataset$classes==levels(dataset$classes)[4],]$Tree_vol_increment.m3.ha.Y.
periods_vec4<-periods_vec[dataset$classes==levels(dataset$classes)[4]]

biom5<-(dataset[dataset$classes==levels(dataset$classes)[5],]$Biomass.g.m2)
start5<-dataset[dataset$classes==levels(dataset$classes)[5],]$Start_day-182
time5<-dataset[dataset$classes==levels(dataset$classes)[5],]$End_day-182
tree5<-dataset[dataset$classes==levels(dataset$classes)[5],]$Tree_vol_increment.m3.ha.Y.
periods_vec5<-periods_vec[dataset$classes==levels(dataset$classes)[5]]

biom6<-(dataset[dataset$classes==levels(dataset$classes)[6],]$Biomass.g.m2)
start6<-dataset[dataset$classes==levels(dataset$classes)[6],]$Start_day-182
time6<-dataset[dataset$classes==levels(dataset$classes)[6],]$End_day-182
tree6<-dataset[dataset$classes==levels(dataset$classes)[6],]$Tree_vol_increment.m3.ha.Y.
periods_vec6<-periods_vec[dataset$classes==levels(dataset$classes)[6]]


###assemble the calibration data for STAN
calib_data<-list(N=length(biom1),
                 max_sd1=sd(biom1),
                 max_sd2=sd(biom2),
                 max_sd3=sd(biom3),
                 max_sd4=sd(biom4),
                 max_sd5=sd(biom5),
                 max_sd6=sd(biom6),
                 biom1=biom1,
                 time1=time1,
                 start1=start1,
                 tree1=tree1,
                 period1=periods_vec1,
                 biom2=biom2,
                 time2=time2,
                 start2=start2,
                 tree2=tree2,
                 period2=periods_vec2,
                 biom3=biom3,
                 time3=time3,
                 start3=start3,
                 tree3=tree3,
                 period3=periods_vec3,
                 biom4=biom4,
                 time4=time4,
                 start4=start4,
                 tree4=tree4,
                 period4=periods_vec4,
                 biom5=biom5,
                 time5=time5,
                 start5=start5,
                 tree5=tree5,
                 period5=periods_vec5,
                 biom6=biom6,
                 time6=time6,
                 start6=start6,
                 tree6=tree6,
                 period6=periods_vec6)



########################################################################################################
########################################  MCMC CALIBRATION  ############################################
Sys.setenv(STAN_NUM_THREADS=4)
iterations=10000
chain_length=floor(iterations/2) #the algorithm will discard the initial iterations for warmup
fit <- stan(file = 'Myc_model_variablemu.stan', data = calib_data,  chains = 4, iter = iterations, cores=3, warmup = iterations-chain_length)
########################################################################################################


# extracting the posteriors in an array
# coordinates of the array are [iterations, chain, parameter]
posteriors_cp <- as.array(fit)
dim(posteriors_cp) # doublecheck the dimensions of the array

class_names<-c("control/apatite", "control/quartz", "control/urea", "phosphorous/apatite", "phosphorous/quartz", "phosphorous/urea") # this is for plotting the names in the graphs



biom_AOV<-rbind(as.data.frame(cbind(rep(class_names[1], length(calib_data$biom1)),calib_data$biom1, calib_data$time1)),
      as.data.frame(cbind(rep(class_names[2], length(calib_data$biom2)),calib_data$biom2, calib_data$time2)),
      as.data.frame(cbind(rep(class_names[3], length(calib_data$biom3)),calib_data$biom3, calib_data$time3)),
      as.data.frame(cbind(rep(class_names[4], length(calib_data$biom4)),calib_data$biom4, calib_data$time4)),
      as.data.frame(cbind(rep(class_names[5], length(calib_data$biom5)),calib_data$biom5, calib_data$time5)),
      as.data.frame(cbind(rep(class_names[6], length(calib_data$biom6)),calib_data$biom6, calib_data$time6)))
      
library(RColorBrewer)
palette_seq=c(brewer.pal(5, "Blues"),
              brewer.pal(5, "Reds"),
              brewer.pal(5, "Purples"),
              brewer.pal(5, "Greens"),
              brewer.pal(5, "Oranges"),
              brewer.pal(5, "Greys"))

par(mar=c(12,4,2,2))
    boxplot(as.numeric(biom_AOV$V2) ~ as.factor(biom_AOV$V3)*as.factor(biom_AOV$V1), las=2, col=palette_seq, xlab="", ylab="biomass")


########################################################################################################
##################################  PLOTTING AND DATA TREATMENT  #######################################

###parameter names vector
parms_names<-c(expression(paste(mu,"1")),
               expression(paste(mu,"2")),
               expression(paste(mu,"3")),
               expression(paste(mu,"4")),"a", "z", expression(paste(P[0]," control/apatite")),
                                         expression(paste(P[0]," control/quartz")),
                                         expression(paste(P[0]," control/urea")),
                                         expression(paste(P[0]," phosphorous/apatite")),
                                         expression(paste(P[0]," phosphorous/quartz")),
                                         expression(paste(P[0]," phosphorous/urea")),
               expression(paste(nu,"1")),
               expression(paste(nu,"2")),
               expression(paste(nu,"3")),
               expression(paste(nu,"4")),
               expression(paste(nu,"5")),
               expression(paste(nu,"6"))) # this is for plotting the names in the graphs
parms_list<-c("mu1","mu2","mu3","mu4", "a", "z", "p1", "p2", "p3", "p4", "p5", "p6", "nu1", "nu2", "nu3", "nu4", "nu5","nu6")


###plot the chains for diagnostic purposes (should be mixing well, if you see some straight line it's really bad)
png("random_walks.png", width=4000, height=2500, res=300)
par(mfrow=c(3,4))
plot(posteriors_cp[,1,"mu1"], type="l", col=chain_palette[1], ylab=parms_names[1], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"mu1"], col=chain_palette[i])}
plot(posteriors_cp[,1,"mu2"], type="l", col=chain_palette[1], ylab=parms_names[2], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"mu2"], col=chain_palette[i])}
plot(posteriors_cp[,1,"mu3"], type="l", col=chain_palette[1], ylab=parms_names[3], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"mu3"], col=chain_palette[i])}
plot(posteriors_cp[,1,"mu4"], type="l", col=chain_palette[1], ylab=parms_names[4], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"mu4"], col=chain_palette[i])}
# plot(posteriors_cp[,1,"a"], type="l", col=chain_palette[1], ylab=parms_names[2], xlab="MCMC index")
# for(i in 2:4){lines(posteriors_cp[,i,"a"], col=chain_palette[i])}
plot(posteriors_cp[,1,"z"], type="l", col=chain_palette[1], ylab=parms_names[6], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"z"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p1"], type="l", col=chain_palette[1], ylab=parms_names[7], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p1"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p2"], type="l", col=chain_palette[1], ylab=parms_names[8], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p2"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p3"], type="l", col=chain_palette[1], ylab=parms_names[9], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p3"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p4"], type="l", col=chain_palette[1], ylab=parms_names[10], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p4"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p5"], type="l", col=chain_palette[1], ylab=parms_names[11], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p5"], col=chain_palette[i])}
plot(posteriors_cp[,1,"p6"], type="l", col=chain_palette[1], ylab=parms_names[12], xlab="MCMC index")
for(i in 2:4){lines(posteriors_cp[,i,"p6"], col=chain_palette[i])}
dev.off()



###plot the parameter density distributions
mu_names<-c(paste("day ", periods[1], "-", periods[2]),
            paste("day ", periods[2], "-", periods[3]),
            paste("day ", periods[3], "-", periods[4]),
            paste("day ", periods[4], "-", 304))

png("parms_distributions.png", width=3500, height=2600, res=350)
par(mfrow=c(2,2))
#layout(matrix(c(1,3,3, 2,3,3), 2, 3, byrow = TRUE))
#layout(matrix(c(1,1,3,3,3, 2,2,3,3,3), 2, 5, byrow = TRUE))

d_mu1<-density(posteriors_cp[,,parms_list[1]]) #kernel density estimation
d_mu2<-density(posteriors_cp[,,parms_list[2]]) #kernel density estimation
d_mu3<-density(posteriors_cp[,,parms_list[3]]) #kernel density estimation
d_mu4<-density(posteriors_cp[,,parms_list[4]]) #kernel density estimation
plot(d_mu1, type="l", col=parms_palette[1], ylab="Probability density", xlab=expression(mu), main="", ylim=c(0, 1.08*max(c(d_mu1$y,d_mu2$y,d_mu3$y,d_mu4$y))), xlim=c(0, 1.08*max(c(d_mu1$x,d_mu2$x,d_mu3$x,d_mu4$x))))
polygon(d_mu1, col=parms_palette_alpha[1], border=parms_palette[1]) #fill the area under the curve with a semitransparent color
text(mlv(d_mu1$x, method = "Parzen"),max(d_mu1$y)*1.01, mu_names[1], col=parms_palette[1], cex=0.8) #add a text on top of the distribution
polygon(d_mu2, col=parms_palette_alpha[2], border=parms_palette[2]) #fill the area under the curve with a semitransparent color
text(mlv(d_mu2$x, method = "Parzen"), max(d_mu2$y)*1.01, mu_names[2], col=parms_palette[2], cex=0.8) #add a text on top of the distribution
polygon(d_mu3, col=parms_palette_alpha[3], border=parms_palette[3]) #fill the area under the curve with a semitransparent color
text(mlv(d_mu3$x, method = "Parzen")*0.8, max(d_mu3$y)*1.01, mu_names[3], col=parms_palette[3], cex=0.8) #add a text on top of the distribution
polygon(d_mu4, col=parms_palette_alpha[4], border=parms_palette[4]) #fill the area under the curve with a semitransparent color
text(mlv(d_mu4$x, method = "Parzen")*0.8, max(d_mu4$y)*1.01, mu_names[4], col=parms_palette[4], cex=0.8) #add a text on top of the distribution
legend("topright","(a)", bty="n")

d_nu1<-density(posteriors_cp[,,parms_list[13]]) #kernel density estimation
d_nu2<-density(posteriors_cp[,,parms_list[14]]) #kernel density estimation
d_nu3<-density(posteriors_cp[,,parms_list[15]]) #kernel density estimation
d_nu4<-density(posteriors_cp[,,parms_list[16]]) #kernel density estimation
d_nu5<-density(posteriors_cp[,,parms_list[17]]) #kernel density estimation
d_nu6<-density(posteriors_cp[,,parms_list[18]]) #kernel density estimation
plot(d_nu1, type="l", col=parms_palette[1], ylab="Probability density", xlab=expression(nu), main="", ylim=c(0, 1.08*max(c(d_nu1$y,d_nu2$y,d_nu3$y,d_nu4$y))), xlim=c(0, 1.08*max(c(d_nu1$x,d_nu2$x,d_nu3$x,d_nu4$x))))
polygon(d_nu1, col=times_palette_alpha[1], border=times_palette[1]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu1$x, method = "Parzen"),max(d_nu1$y)*1.01, nu_names[1], col=times_palette[1], cex=0.8) #add a text on top of the distribution
polygon(d_nu2, col=times_palette_alpha[2], border=class_names[2]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu2$x, method = "Parzen"), max(d_nu2$y)*1.01, class_names[2], col=times_palette[2], cex=0.8) #add a text on top of the distribution
polygon(d_nu3, col=times_palette_alpha[3], border=times_palette[3]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu3$x, method = "Parzen")*0.8, max(d_nu3$y)*1.01, class_names[3], col=times_palette[3], cex=0.8) #add a text on top of the distribution
polygon(d_nu4, col=times_palette_alpha[4], border=times_palette[4]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu4$x, method = "Parzen")*0.8, max(d_nu4$y)*1.01, class_names[4], col=times_palette[4], cex=0.8) #add a text on top of the distribution
polygon(d_nu5, col=times_palette_alpha[5], border=times_palette[4]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu5$x, method = "Parzen")*0.8, max(d_nu4$y)*1.01, class_names[5], col=times_palette[4], cex=0.8) #add a text on top of the distribution
polygon(d_nu6, col=times_palette_alpha[6], border=times_palette[4]) #fill the area under the curve with a semitransparent color
text(mlv(d_nu6$x, method = "Parzen")*0.8, max(d_nu4$y)*1.01, class_names[6], col=times_palette[4], cex=0.8) #add a text on top of the distribution
legend("topright","(a)", bty="n")

d_z<-density(posteriors_cp[,,parms_list[6]]) #kernel density estimation
plot(d_z, type="l", col=parms_palette[5], ylab="Probability density", xlab=parms_names[6], main="")
polygon(d_z, col=parms_palette_alpha[5], border=parms_palette[5]) #fill the area under the curve with a semitransparent color
legend("topright","(b)", bty="n")


#one plot for all the p values
d_p1<-density(posteriors_cp[,,parms_list[7]]) #kernel density estimation
d_p2<-density(posteriors_cp[,,parms_list[8]]) #kernel density estimation
d_p3<-density(posteriors_cp[,,parms_list[9]]) #kernel density estimation
d_p4<-density(posteriors_cp[,,parms_list[10]]) #kernel density estimation
d_p5<-density(posteriors_cp[,,parms_list[11]]) #kernel density estimation
d_p6<-density(posteriors_cp[,,parms_list[12]]) #kernel density estimation

plot(d_p1, type="l", col=times_palette[1], ylab="Probability density", xlab=expression(P[0]), main="", xlim=c(0, max(c(d_p1$x,d_p2$x,d_p3$x,d_p4$x,d_p5$x,d_p6$x))),
     ylim=c(0, max(c(d_p1$y,d_p2$y,d_p3$y,d_p4$y,d_p5$y,d_p6$y))*1.05))
polygon(d_p1, col=times_palette_alpha[1], border=times_palette[1]) #fill the area under the curve with a semitransparent color
text(mlv(d_p1$x, method = "Parzen"), max(d_p1$y)*1.05, class_names[1], col=times_palette[1], cex=0.8) #add a text on top of the distribution
polygon(d_p2, col=times_palette_alpha[2], border=times_palette[2]) #fill the area under the curve with a semitransparent color
text(mlv(d_p2$x, method="Parzen"), max(d_p2$y)*1.05, class_names[2], col=times_palette[2], cex=0.8) #add a text on top of the distribution
polygon(d_p3, col=times_palette_alpha[3], border=times_palette[3]) #fill the area under the curve with a semitransparent color
text(mlv(d_p3$x, method="Parzen"), max(d_p3$y)*1.05, class_names[3], col=times_palette[3], cex=0.8) #add a text on top of the distribution
polygon(d_p4, col=times_palette_alpha[4], border=times_palette[4]) #fill the area under the curve with a semitransparent color
text(mlv(d_p4$x, method="Parzen"), max(d_p4$y)*1.03, class_names[4], col=times_palette[4], cex=0.8) #add a text on top of the distribution
polygon(d_p5, col=times_palette_alpha[5], border=times_palette[5]) #fill the area under the curve with a semitransparent color
text(mlv(d_p5$x, method="Parzen"), max(d_p5$y)*1.03, class_names[5], col=times_palette[5], cex=0.8) #add a text on top of the distribution
polygon(d_p6, col=times_palette_alpha[6], border=times_palette[6]) #fill the area under the curve with a semitransparent color
text(mlv(d_p6$x, method="Parzen"), max(d_p6$y)*1.03, class_names[6], col=times_palette[6], cex=0.8) #add a text on top of the distribution
legend("topright","(c)", bty="n")

dev.off()






########################################################################################################
###########################################  SIMULATION ################################################


posteriors<-as.data.frame(fit) #extract the posterios in an array of two dimensions (collapsing the chains)

# thinning the posteriors to plot the simulations
thinning=10
resampling_vector<-seq(from=1, to=chain_length*4, by=thinning)
resampled_posteriors<-posteriors[resampling_vector,,]

dim(posteriors)
dim(resampled_posteriors)

# plot the simulations
time_sim=180
time_simulation=seq(from=0, to=time_sim)+153


#initial plotting of the data (to get the idea of what the data looks like)
biom_range<-range(c(calib_data$biom1, calib_data$biom2, calib_data$biom3, calib_data$biom4, calib_data$biom5, calib_data$biom6))

png("Data_plot.png", width=4000, height=2000, res=350)
par(mfrow=c(2,3), mar=c(5,5,2,2))
#print("Empty plot specify the axes limits of the graphic:")
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[1])
for(i in 1:length(calib_data$start1)){Arrows(calib_data$start1[i], 0, calib_data$time1[i], calib_data$biom1[i], col=times_palette_alpha[1], pch=16, arr.width=0.1, arr.length=0.1)}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[2])
for(i in 1:length(calib_data$start2)){Arrows(calib_data$start2[i], 0, calib_data$time2[i], calib_data$biom2[i], col=times_palette_alpha[2], pch=16, arr.width=0.1, arr.length=0.1)}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[3])
for(i in 1:length(calib_data$start3)){Arrows(calib_data$start3[i], 0, calib_data$time3[i], calib_data$biom3[i], col=times_palette_alpha[3], pch=16, arr.width=0.1, arr.length=0.1)}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[4])
for(i in 1:length(calib_data$start4)){Arrows(calib_data$start4[i], 0, calib_data$time4[i], calib_data$biom4[i], col=times_palette_alpha[4], pch=16, arr.width=0.1, arr.length=0.1)}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[5])
for(i in 1:length(calib_data$start5)){Arrows(calib_data$start5[i], 0, calib_data$time5[i], calib_data$biom5[i], col=times_palette_alpha[5], pch=16, arr.width=0.1, arr.length=0.1)}
plot(1, type="n", xlab="", ylab="", xlim=c(0, 180), ylim=c(0, biom_range[2]),main=class_names[6])
for(i in 1:length(calib_data$start6)){Arrows(calib_data$start6[i], 0, calib_data$time6[i], calib_data$biom6[i], col=times_palette_alpha[6], pch=16, arr.width=0.1, arr.length=0.1)}
dev.off()

#plotting the rates
png("Boxplot_rates.png", width=2500, height=2500, res=350)
par(mar=c(10,5,2,2))
boxplot(((dataset$Biomass.g.m2)/(dataset$End_day- dataset$Start_day))~dataset$classes, col=times_palette_alpha, pch=16, arr.width=0.1, arr.length=0.1, las=2, xlab="", ylab=expression(paste("Average biomass increase rate, g m"^"-2")))
dev.off()

png("Boxplot_rates_times.png", width=2500, height=2500, res=350)
par(mar=c(10,5,2,2))
boxplot(((dataset$Biomass.g.m2)/(dataset$End_day- dataset$Start_day))~dataset$classes*dataset$Incubation.time..months.,
        col=times_palette_alpha, pch=16, arr.width=0.1, arr.length=0.1, las=2, xlab="", ylab=expression(paste("Average biomass increase rate, g m"^"-2")))
dev.off()


#create the table to hold the posterior probability density for each of the data points

biomass_simulated_p1<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))
biomass_simulated_p2<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))
biomass_simulated_p3<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))
biomass_simulated_p4<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))
biomass_simulated_p5<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))
biomass_simulated_p6<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$biom1))

RMSE_simulated_p1<-mat.or.vec(dim(resampled_posteriors)[1], 1)
RMSE_simulated_p2<-mat.or.vec(dim(resampled_posteriors)[1], 1)
RMSE_simulated_p3<-mat.or.vec(dim(resampled_posteriors)[1], 1)
RMSE_simulated_p4<-mat.or.vec(dim(resampled_posteriors)[1], 1)
RMSE_simulated_p5<-mat.or.vec(dim(resampled_posteriors)[1], 1)
RMSE_simulated_p6<-mat.or.vec(dim(resampled_posteriors)[1], 1)

res_simulated_p1<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time1))
res_simulated_p2<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time2))
res_simulated_p3<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time3))
res_simulated_p4<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time4))
res_simulated_p5<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time5))
res_simulated_p6<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time6))

cv_simulated_p1<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time1))
cv_simulated_p2<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time2))
cv_simulated_p3<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time3))
cv_simulated_p4<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time4))
cv_simulated_p5<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time5))
cv_simulated_p6<-mat.or.vec(dim(resampled_posteriors)[1], length(calib_data$time6))


  #run the model for each parameter set
  for( i in 1: dim(resampled_posteriors)[1]){

   #calculate the biomass
    p1=resampled_posteriors$p1[i]
    p2=resampled_posteriors$p2[i]
    p3=resampled_posteriors$p3[i]
    p4=resampled_posteriors$p4[i]
    p5=resampled_posteriors$p5[i]
    p6=resampled_posteriors$p6[i]
    mu1=resampled_posteriors$mu1[i]
    mu2=resampled_posteriors$mu2[i]
    mu3=resampled_posteriors$mu3[i]
    mu4=resampled_posteriors$mu4[i]
    om=(2*pi)/365.25
    z=resampled_posteriors$z[i]
    nu1=resampled_posteriors$nu1[i]
    nu2=resampled_posteriors$nu2[i]
    nu3=resampled_posteriors$nu3[i]
    nu4=resampled_posteriors$nu4[i]
    nu5=resampled_posteriors$nu5[i]
    nu6=resampled_posteriors$nu6[i]

    a1=z*mean(tree1);
    for(k in 1:length(periods_vec1)){
        if(periods_vec1[k]==1){
        biomass_simulated_p1[i,k]<- p1/(mu1*nu1)*(1-exp(-(mu1*nu1)*(calib_data$time1[k]-calib_data$start1[k])))+
          ((a1*p1)/(mu1*nu1))*(cos(om)*calib_data$time1[k]-exp(-(mu1*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k])-
          ((a1*p1)/(om*(mu1*nu1)))*(1/(1+(1/(om^2*(mu1*nu1)^2))))*((1/(mu1*nu1))*(sin(om)*calib_data$time1[k]-exp(-(mu1*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*sin(om)*calib_data$start1[k])+
                                                      (1/(om*(mu1*nu1)^2))*(cos(om)*calib_data$time1[k]-exp(-(mu1*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k]));
        } else if (periods_vec1[k]==2){
          biomass_simulated_p1[i,k]<- p1/(mu2*nu1)*(1-exp(-(mu2*nu1)*(calib_data$time1[k]-calib_data$start1[k])))+
            ((a1*p1)/(mu2*nu1))*(cos(om)*calib_data$time1[k]-exp(-(mu2*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k])-
            ((a1*p1)/(om*(mu2*nu1)))*(1/(1+(1/(om^2*(mu2*nu1)^2))))*((1/(mu2*nu1))*(sin(om)*calib_data$time1[k]-exp(-(mu2*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*sin(om)*calib_data$start1[k])+
                                                           (1/(om*(mu2*nu1)^2))*(cos(om)*calib_data$time1[k]-exp(-(mu2*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k]));
        } else if (periods_vec1[k]==3){
          biomass_simulated_p1[i,k]<- p1/(mu3*nu1)*(1-exp(-(mu3*nu1)*(calib_data$time1[k]-calib_data$start1[k])))+
            ((a1*p1)/(mu3*nu1))*(cos(om)*calib_data$time1[k]-exp(-(mu3*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k])-
            ((a1*p1)/(om*(mu3*nu1)))*(1/(1+(1/(om^2*(mu3*nu1)^2))))*((1/(mu3*nu1))*(sin(om)*calib_data$time1[k]-exp(-(mu3*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*sin(om)*calib_data$start1[k])+
                                                           (1/(om*(mu3*nu1)^2))*(cos(om)*calib_data$time1[k]-exp(-(mu3*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k]));
        } else if (periods_vec1[k]==4){
          biomass_simulated_p1[i,k]<- p1/(mu4*nu1)*(1-exp(-(mu4*nu1)*(calib_data$time1[k]-calib_data$start1[k])))+
            ((a1*p1)/(mu4*nu1))*(cos(om)*calib_data$time1[k]-exp(-(mu4*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k])-
            ((a1*p1)/(om*(mu4*nu1)))*(1/(1+(1/(om^2*(mu4*nu1)^2))))*((1/(mu4*nu1))*(sin(om)*calib_data$time1[k]-exp(-(mu4*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*sin(om)*calib_data$start1[k])+
                                                           (1/(om*(mu4*nu1)^2))*(cos(om)*calib_data$time1[k]-exp(-(mu4*nu1)*(calib_data$time1[k]-calib_data$start1[k]))*cos(om)*calib_data$start1[k]));
        }
    }
    a2=z*mean(tree2);
      for(k in 1:length(periods_vec1)){
        if(periods_vec2[k]==1){
        biomass_simulated_p2[i,k]<- p2/(mu1*nu2)*(1-exp(-(mu1*nu2)*(calib_data$time2[k]-calib_data$start2[k])))+
        ((a2*p2)/(mu1*nu2))*(cos(om)*calib_data$time2[k]-exp(-(mu1*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k])-
        ((a2*p2)/(om*(mu1*nu2)))*(1/(1+(1/(om^2*(mu1*nu2)^2))))*((1/(mu1*nu2))*(sin(om)*calib_data$time2[k]-exp(-(mu1*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*sin(om)*calib_data$start2[k])+
                                                    (1/(om*(mu1*nu2)^2))*(cos(om)*calib_data$time2[k]-exp(-(mu1*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k]));
      } else if (periods_vec2[k]==2){
        biomass_simulated_p2[i,k]<- p2/(mu2*nu2)*(1-exp(-(mu2*nu2)*(calib_data$time2[k]-calib_data$start2[k])))+
          ((a2*p2)/(mu2*nu2))*(cos(om)*calib_data$time2[k]-exp(-(mu2*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k])-
          ((a2*p2)/(om*(mu2*nu2)))*(1/(1+(1/(om^2*(mu2*nu2)^2))))*((1/(mu2*nu2))*(sin(om)*calib_data$time2[k]-exp(-(mu2*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*sin(om)*calib_data$start2[k])+
                                                         (1/(om*(mu2*nu2)^2))*(cos(om)*calib_data$time2[k]-exp(-(mu2*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k]));
      } else if (periods_vec2[k]==3){
        biomass_simulated_p2[i,k]<- p2/(mu3*nu2)*(1-exp(-(mu3*nu2)*(calib_data$time2[k]-calib_data$start2[k])))+
          ((a2*p2)/(mu3*nu2))*(cos(om)*calib_data$time2[k]-exp(-(mu3*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k])-
          ((a2*p2)/(om*(mu3*nu2)))*(1/(1+(1/(om^2*(mu3*nu2)^2))))*((1/(mu3*nu2))*(sin(om)*calib_data$time2[k]-exp(-(mu3*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*sin(om)*calib_data$start2[k])+
                                                         (1/(om*(mu3*nu2)^2))*(cos(om)*calib_data$time2[k]-exp(-(mu3*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k]));
      } else if (periods_vec2[k]==4){
        biomass_simulated_p2[i,k]<- p2/(mu4*nu2)*(1-exp(-(mu4*nu2)*(calib_data$time2[k]-calib_data$start2[k])))+
          ((a2*p2)/(mu4*nu2))*(cos(om)*calib_data$time2[k]-exp(-(mu4*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k])-
          ((a2*p2)/(om*(mu4*nu2)))*(1/(1+(1/(om^2*(mu4*nu2)^2))))*((1/(mu4*nu2))*(sin(om)*calib_data$time2[k]-exp(-(mu4*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*sin(om)*calib_data$start2[k])+
                                                         (1/(om*(mu4*nu2)^2))*(cos(om)*calib_data$time2[k]-exp(-(mu4*nu2)*(calib_data$time2[k]-calib_data$start2[k]))*cos(om)*calib_data$start2[k]));
      }
    }

    a3=z*mean(tree3);
    for(k in 1:length(periods_vec1)){
        if(periods_vec3[k]==1){
        biomass_simulated_p3[i,k]<- p3/(mu1*nu3)*(1-exp(-(mu1*nu3)*(calib_data$time3[k]-calib_data$start3[k])))+
        ((a3*p3)/(mu1*nu3))*(cos(om)*calib_data$time3[k]-exp(-(mu1*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k])-
        ((a3*p3)/(om*(mu1*nu3)))*(1/(1+(1/(om^2*(mu1*nu3)^2))))*((1/(mu1*nu3))*(sin(om)*calib_data$time3[k]-exp(-(mu1*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*sin(om)*calib_data$start3[k])+
                                                    (1/(om*(mu1*nu3)^2))*(cos(om)*calib_data$time3[k]-exp(-(mu1*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k]));
      } else if (periods_vec3[k]==2){
        biomass_simulated_p3[i,k]<- p3/(mu2*nu3)*(1-exp(-(mu2*nu3)*(calib_data$time3[k]-calib_data$start3[k])))+
          ((a3*p3)/(mu2*nu3))*(cos(om)*calib_data$time3[k]-exp(-(mu2*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k])-
          ((a3*p3)/(om*(mu2*nu3)))*(1/(1+(1/(om^2*(mu2*nu3)^2))))*((1/(mu2*nu3))*(sin(om)*calib_data$time3[k]-exp(-(mu2*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*sin(om)*calib_data$start3[k])+
                                                         (1/(om*(mu2*nu3)^2))*(cos(om)*calib_data$time3[k]-exp(-(mu2*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k]));
      } else if (periods_vec3[k]==3){
        biomass_simulated_p3[i,k]<- p3/(mu3*nu3)*(1-exp(-(mu3*nu3)*(calib_data$time3[k]-calib_data$start3[k])))+
          ((a3*p3)/(mu3*nu3))*(cos(om)*calib_data$time3[k]-exp(-(mu3*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k])-
          ((a3*p3)/(om*(mu3*nu3)))*(1/(1+(1/(om^2*(mu3*nu3)^2))))*((1/(mu3*nu3))*(sin(om)*calib_data$time3[k]-exp(-(mu3*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*sin(om)*calib_data$start3[k])+
                                                         (1/(om*(mu3*nu3)^2))*(cos(om)*calib_data$time3[k]-exp(-(mu3*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k]));
      } else if (periods_vec3[k]==4){
        biomass_simulated_p3[i,k]<- p3/(mu4*nu3)*(1-exp(-(mu4*nu3)*(calib_data$time3[k]-calib_data$start3[k])))+
          ((a3*p3)/(mu4*nu3))*(cos(om)*calib_data$time3[k]-exp(-(mu4*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k])-
          ((a3*p3)/(om*(mu4*nu3)))*(1/(1+(1/(om^2*(mu4*nu3)^2))))*((1/(mu4*nu3))*(sin(om)*calib_data$time3[k]-exp(-(mu4*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*sin(om)*calib_data$start3[k])+
                                                         (1/(om*(mu4*nu3)^2))*(cos(om)*calib_data$time3[k]-exp(-(mu4*nu3)*(calib_data$time3[k]-calib_data$start3[k]))*cos(om)*calib_data$start3[k]));
      }
    }

    a4=z*mean(tree4);
    for(k in 1:length(periods_vec1)){
        if(periods_vec4[k]==1){
        biomass_simulated_p4[i,k]<- p4/(mu1*nu4)*(1-exp(-(mu1*nu4)*(calib_data$time4[k]-calib_data$start4[k])))+
        ((a4*p4)/(mu1*nu4))*(cos(om)*calib_data$time4[k]-exp(-(mu1*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k])-
        ((a4*p4)/(om*(mu1*nu4)))*(1/(1+(1/(om^2*(mu1*nu4)^2))))*((1/(mu1*nu4))*(sin(om)*calib_data$time4[k]-exp(-(mu1*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*sin(om)*calib_data$start4[k])+
                                                    (1/(om*(mu1*nu4)^2))*(cos(om)*calib_data$time4[k]-exp(-(mu1*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k]));
      } else if (periods_vec4[k]==2){
        biomass_simulated_p4[i,k]<- p4/(mu2*nu4)*(1-exp(-(mu2*nu4)*(calib_data$time4[k]-calib_data$start4[k])))+
          ((a4*p4)/(mu2*nu4))*(cos(om)*calib_data$time4[k]-exp(-(mu2*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k])-
          ((a4*p4)/(om*(mu2*nu4)))*(1/(1+(1/(om^2*(mu2*nu4)^2))))*((1/(mu2*nu4))*(sin(om)*calib_data$time4[k]-exp(-(mu2*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*sin(om)*calib_data$start4[k])+
                                                         (1/(om*(mu2*nu4)^2))*(cos(om)*calib_data$time4[k]-exp(-(mu2*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k]));
      } else if (periods_vec4[k]==3){
        biomass_simulated_p4[i,k]<- p4/(mu3*nu4)*(1-exp(-(mu3*nu4)*(calib_data$time4[k]-calib_data$start4[k])))+
          ((a4*p4)/(mu3*nu4))*(cos(om)*calib_data$time4[k]-exp(-(mu3*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k])-
          ((a4*p4)/(om*(mu3*nu4)))*(1/(1+(1/(om^2*(mu3*nu4)^2))))*((1/(mu3*nu4))*(sin(om)*calib_data$time4[k]-exp(-(mu3*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*sin(om)*calib_data$start4[k])+
                                                         (1/(om*(mu3*nu4)^2))*(cos(om)*calib_data$time4[k]-exp(-(mu3*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k]));
      } else if (periods_vec4[k]==4){
        biomass_simulated_p4[i,k]<- p4/(mu4*nu4)*(1-exp(-(mu4*nu4)*(calib_data$time4[k]-calib_data$start4[k])))+
          ((a4*p4)/(mu4*nu4))*(cos(om)*calib_data$time4[k]-exp(-(mu4*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k])-
          ((a4*p4)/(om*(mu4*nu4)))*(1/(1+(1/(om^2*(mu4*nu4)^2))))*((1/(mu4*nu4))*(sin(om)*calib_data$time4[k]-exp(-(mu4*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*sin(om)*calib_data$start4[k])+
                                                         (1/(om*(mu4*nu4)^2))*(cos(om)*calib_data$time4[k]-exp(-(mu4*nu4)*(calib_data$time4[k]-calib_data$start4[k]))*cos(om)*calib_data$start4[k]));
      }
    }

    a5=z*mean(tree5);
    for(k in 1:length(periods_vec1)){
      if(periods_vec5[k]==1){
        biomass_simulated_p5[i,k]<- p5/(mu1*nu5)*(1-exp(-(mu1*nu5)*(calib_data$time5[k]-calib_data$start5[k])))+
        ((a5*p5)/(mu1*nu5))*(cos(om)*calib_data$time5[k]-exp(-(mu1*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k])-
        ((a5*p5)/(om*(mu1*nu5)))*(1/(1+(1/(om^2*(mu1*nu5)^2))))*((1/(mu1*nu5))*(sin(om)*calib_data$time5[k]-exp(-(mu1*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*sin(om)*calib_data$start5[k])+
                                                     (1/(om*(mu1*nu5)^2))*(cos(om)*calib_data$time5[k]-exp(-(mu1*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k]));
      } else if (periods_vec5[k]==2){
        biomass_simulated_p5[i,k]<- p5/(mu2*nu5)*(1-exp(-(mu2*nu5)*(calib_data$time5[k]-calib_data$start5[k])))+
          ((a5*p5)/(mu2*nu5))*(cos(om)*calib_data$time5[k]-exp(-(mu2*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k])-
          ((a5*p5)/(om*(mu2*nu5)))*(1/(1+(1/(om^2*(mu2*nu5)^2))))*((1/(mu2*nu5))*(sin(om)*calib_data$time5[k]-exp(-(mu2*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*sin(om)*calib_data$start5[k])+
                                                         (1/(om*(mu2*nu5)^2))*(cos(om)*calib_data$time5[k]-exp(-(mu2*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k]));
      } else if (periods_vec5[k]==3){
        biomass_simulated_p5[i,k]<- p5/(mu3*nu5)*(1-exp(-(mu3*nu5)*(calib_data$time5[k]-calib_data$start5[k])))+
          ((a5*p5)/(mu3*nu5))*(cos(om)*calib_data$time5[k]-exp(-(mu3*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k])-
          ((a5*p5)/(om*(mu3*nu5)))*(1/(1+(1/(om^2*(mu3*nu5)^2))))*((1/(mu3*nu5))*(sin(om)*calib_data$time5[k]-exp(-(mu3*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*sin(om)*calib_data$start5[k])+
                                                         (1/(om*(mu3*nu5)^2))*(cos(om)*calib_data$time5[k]-exp(-(mu3*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k]));
      } else if (periods_vec5[k]==4){
        biomass_simulated_p5[i,k]<- p5/(mu4*nu5)*(1-exp(-(mu4*nu5)*(calib_data$time5[k]-calib_data$start5[k])))+
          ((a5*p5)/(mu4*nu5))*(cos(om)*calib_data$time5[k]-exp(-(mu4*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k])-
          ((a5*p5)/(om*(mu4*nu5)))*(1/(1+(1/(om^2*(mu4*nu5)^2))))*((1/(mu4*nu5))*(sin(om)*calib_data$time5[k]-exp(-(mu4*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*sin(om)*calib_data$start5[k])+
                                                         (1/(om*(mu4*nu5)^2))*(cos(om)*calib_data$time5[k]-exp(-(mu4*nu5)*(calib_data$time5[k]-calib_data$start5[k]))*cos(om)*calib_data$start5[k]));
      }
    }

    a6=z*mean(tree6);
    for(k in 1:length(periods_vec1)){
      if(periods_vec6[k]==1){
        biomass_simulated_p6[i,k]<- p6/(mu1*nu6)*(1-exp(-(mu1*nu6)*(calib_data$time6[k]-calib_data$start6[k])))+
        ((a6*p6)/(mu1*nu6))*(cos(om)*calib_data$time6[k]-exp(-(mu1*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k])-
        ((a6*p6)/(om*(mu1*nu6)))*(1/(1+(1/(om^2*(mu1*nu6)^2))))*((1/(mu1*nu6))*(sin(om)*calib_data$time6[k]-exp(-(mu1*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*sin(om)*calib_data$start6[k])+
                                                     (1/(om*(mu1*nu6)^2))*(cos(om)*calib_data$time6[k]-exp(-(mu1*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k]));
      } else if(periods_vec6[k]==2){
        biomass_simulated_p6[i,k]<- p6/(mu2*nu6)*(1-exp(-(mu2*nu6)*(calib_data$time6[k]-calib_data$start6[k])))+
          ((a6*p6)/(mu2*nu6))*(cos(om)*calib_data$time6[k]-exp(-(mu2*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k])-
          ((a6*p6)/(om*(mu2*nu6)))*(1/(1+(1/(om^2*(mu2*nu6)^2))))*((1/(mu2*nu6))*(sin(om)*calib_data$time6[k]-exp(-(mu2*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*sin(om)*calib_data$start6[k])+
                                                         (1/(om*(mu2*nu6)^2))*(cos(om)*calib_data$time6[k]-exp(-(mu2*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k]));
      } else if(periods_vec6[k]==3){
        biomass_simulated_p6[i,k]<- p6/(mu3*nu6)*(1-exp(-(mu3*nu6)*(calib_data$time6[k]-calib_data$start6[k])))+
          ((a6*p6)/(mu3*nu6))*(cos(om)*calib_data$time6[k]-exp(-(mu3*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k])-
          ((a6*p6)/(om*(mu3*nu6)))*(1/(1+(1/(om^2*(mu3*nu6)^2))))*((1/(mu3*nu6))*(sin(om)*calib_data$time6[k]-exp(-(mu3*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*sin(om)*calib_data$start6[k])+
                                                         (1/(om*(mu3*nu6)^2))*(cos(om)*calib_data$time6[k]-exp(-(mu3*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k]));
      } else if(periods_vec6[k]==4){
        biomass_simulated_p6[i,k]<- p6/(mu4*nu6)*(1-exp(-(mu4*nu6)*(calib_data$time6[k]-calib_data$start6[k])))+
          ((a6*p6)/(mu4*nu6))*(cos(om)*calib_data$time6[k]-exp(-(mu4*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k])-
          ((a6*p6)/(om*(mu4*nu6)))*(1/(1+(1/(om^2*(mu4*nu6)^2))))*((1/(mu4*nu6))*(sin(om)*calib_data$time6[k]-exp(-(mu4*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*sin(om)*calib_data$start6[k])+
                                                         (1/(om*(mu4*nu6)^2))*(cos(om)*calib_data$time6[k]-exp(-(mu4*nu6)*(calib_data$time6[k]-calib_data$start6[k]))*cos(om)*calib_data$start6[k]));
      }
    }

    # #calculate the RMSE for each parameter set
    RMSE_simulated_p1[i]<-rmse(biomass_simulated_p1[i,],calib_data$biom1)
    RMSE_simulated_p2[i]<-rmse(biomass_simulated_p2[i,],calib_data$biom2)
    RMSE_simulated_p3[i]<-rmse(biomass_simulated_p3[i,],calib_data$biom3)
    RMSE_simulated_p4[i]<-rmse(biomass_simulated_p4[i,],calib_data$biom4)
    RMSE_simulated_p5[i]<-rmse(biomass_simulated_p5[i,],calib_data$biom5)
    RMSE_simulated_p6[i]<-rmse(biomass_simulated_p6[i,],calib_data$biom6)

    res_simulated_p1[i,]<-abs(biomass_simulated_p1[i,] - calib_data$biom1)
    res_simulated_p2[i,]<-abs(biomass_simulated_p2[i,] - calib_data$biom2)
    res_simulated_p3[i,]<-abs(biomass_simulated_p3[i,] - calib_data$biom3)
    res_simulated_p4[i,]<-abs(biomass_simulated_p4[i,] - calib_data$biom4)
    res_simulated_p5[i,]<-abs(biomass_simulated_p5[i,] - calib_data$biom5)
    res_simulated_p6[i,]<-abs(biomass_simulated_p6[i,] - calib_data$biom6)


    cv_simulated_p1[i,]<-abs(biomass_simulated_p1[i,] - calib_data$biom1)/mean(calib_data$biom1)
    cv_simulated_p2[i,]<-abs(biomass_simulated_p2[i,] - calib_data$biom2)/mean(calib_data$biom2)
    cv_simulated_p3[i,]<-abs(biomass_simulated_p3[i,] - calib_data$biom3)/mean(calib_data$biom3)
    cv_simulated_p4[i,]<-abs(biomass_simulated_p4[i,] - calib_data$biom4)/mean(calib_data$biom4)
    cv_simulated_p5[i,]<-abs(biomass_simulated_p5[i,] - calib_data$biom5)/mean(calib_data$biom5)
    cv_simulated_p6[i,]<-abs(biomass_simulated_p6[i,] - calib_data$biom6)/mean(calib_data$biom6)

  }


dim(biomass_simulated_p2)

plot(calib_data$biom2, biomass_simulated_p2[1500,], ylim=c(0, biom_range[2]), col=times_palette[1], pch=15, main=class_names[1], ylab=expression(paste("Biomass, g m"^"-2")), xlab="")
plot(calib_data$biom2, biomass_simulated_p2[1701,], ylim=c(0, biom_range[2]), col=times_palette[1], pch=15, main=class_names[1], ylab=expression(paste("Biomass, g m"^"-2")), xlab="")

plot(calib_data$biom1, biomass_simulated_p1[1500,], ylim=c(0, biom_range[2]), col=times_palette[1], pch=15, main=class_names[1], ylab=expression(paste("Biomass, g m"^"-2")), xlab="")
plot(calib_data$biom1, biomass_simulated_p1[1001,], ylim=c(0, biom_range[2]), col=times_palette[1], pch=15, main=class_names[1], ylab=expression(paste("Biomass, g m"^"-2")), xlab="")


biom_range<-range(c(calib_data$biom1, calib_data$biom2, calib_data$biom3, calib_data$biom4, calib_data$biom5, calib_data$biom6),
                  c(biomass_simulated_p1, biomass_simulated_p2, biomass_simulated_p3, biomass_simulated_p4, biomass_simulated_p5, biomass_simulated_p6))



#ID for each point
colnames(dataset)
ID1<-(dataset[dataset$classes==levels(dataset$classes)[1],]$No.1)
ID2<-(dataset[dataset$classes==levels(dataset$classes)[2],]$No.1)
ID3<-(dataset[dataset$classes==levels(dataset$classes)[3],]$No.1)
ID4<-(dataset[dataset$classes==levels(dataset$classes)[4],]$No.1)
ID5<-(dataset[dataset$classes==levels(dataset$classes)[5],]$No.1)
ID6<-(dataset[dataset$classes==levels(dataset$classes)[6],]$No.1)



#plotting the biomass simulation together with the data
png("Scatter_plot.png", width=4000, height=2500, res=350)
par(mfrow=c(2,3), mar=c(5,5,2,2))
means<-colMeans(biomass_simulated_p1)
sds<-colSdApply(biomass_simulated_p1)
plot(calib_data$biom1, means, ylim=c(0, max(c(means,calib_data$biom1))), xlim=c(0, biom_range[2]), col=times_palette[1], main=class_names[1], pch=as.numeric(as.factor(round(calib_data$time1/10)*10)), ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom1, means-sds, calib_data$biom1, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[1])
legend("topright","(a)", bty="n")
legend("topleft", levels(as.factor(round(calib_data$time1/10)*10)), pch=1:5,cex=1.2, bty="n")


means<-colMeans(biomass_simulated_p2)
sds<-colSdApply(biomass_simulated_p2)
plot(calib_data$biom2, means, ylim=c(0, max(c(means,calib_data$biom2))), xlim=c(0, biom_range[2]), col=times_palette[2], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[2], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom2, means-sds, calib_data$biom2, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[2])
legend("topright","(b)", bty="n")

means<-colMeans(biomass_simulated_p3)
sds<-colSdApply(biomass_simulated_p3)
plot(calib_data$biom3, means, ylim=c(0, max(c(means,calib_data$biom3))), xlim=c(0, biom_range[2]), col=times_palette[3], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[3], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom3, means-sds, calib_data$biom3, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[3])
legend("topright","(c)", bty="n")

means<-colMeans(biomass_simulated_p4)
sds<-colSdApply(biomass_simulated_p4)
plot(calib_data$biom4, means, ylim=c(0, max(c(means,calib_data$biom4))), xlim=c(0, biom_range[2]), col=times_palette[4], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[4], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom4, means-sds, calib_data$biom4, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[4])
legend("topright","(d)", bty="n")

means<-colMeans(biomass_simulated_p5)
sds<-colSdApply(biomass_simulated_p5)
plot(calib_data$biom5, means, ylim=c(0, max(c(means,calib_data$biom5))), xlim=c(0, biom_range[2]), col=times_palette[5], pch=as.numeric(as.factor(round(calib_data$time3/10)*10)), main=class_names[5], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom5, means-sds, calib_data$biom5, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[5])
legend("topright","(e)", bty="n")

means<-colMeans(biomass_simulated_p6)
sds<-colSdApply(biomass_simulated_p6)
plot(calib_data$biom6, means, ylim=c(0, max(c(means,calib_data$biom6))), xlim=c(0, biom_range[2]), col=times_palette[6], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[6], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
arrows(calib_data$biom6, means-sds, calib_data$biom6, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[6])
legend("topright","(f)", bty="n")

dev.off()




#plotting the biomass simulation together with the data
png("Scatter_plot_with_IDS.png", width=4000, height=2500, res=390)
par(mfrow=c(2,3), mar=c(5,5,2,2))
means<-colMeans(biomass_simulated_p1)
sds<-colSdApply(biomass_simulated_p1)
plot(calib_data$biom1, means, ylim=c(0, max(c(means,calib_data$biom1))), xlim=c(0, biom_range[2]), col=times_palette[1], main=class_names[1], pch=as.numeric(as.factor(round(calib_data$time1/10)*10)), ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom1, means-sds, calib_data$biom1, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[1])
text(calib_data$biom1, means+.5, ID1, cex=0.5)
legend("topright","(a)", bty="n")
legend("topleft", levels(as.factor(round(calib_data$time1/10)*10)), pch=1:5,cex=1.2, bty="n")

means<-colMeans(biomass_simulated_p2)
sds<-colSdApply(biomass_simulated_p2)
plot(calib_data$biom2, means, ylim=c(0, max(c(means,calib_data$biom2))), xlim=c(0, biom_range[2]), col=times_palette[2], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[2], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom2, means-sds, calib_data$biom2, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[2])
legend("topright","(b)", bty="n")
text(calib_data$biom2, means+.5, ID2, cex=0.5)

means<-colMeans(biomass_simulated_p3)
sds<-colSdApply(biomass_simulated_p3)
plot(calib_data$biom3, means, ylim=c(0, max(c(means,calib_data$biom3))), xlim=c(0, biom_range[2]), col=times_palette[3], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[3], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom3, means-sds, calib_data$biom3, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[3])
legend("topright","(c)", bty="n")
text(calib_data$biom3, means+.5, ID3, cex=0.5)

means<-colMeans(biomass_simulated_p4)
sds<-colSdApply(biomass_simulated_p4)
plot(calib_data$biom4, means, ylim=c(0, max(c(means,calib_data$biom4))), xlim=c(0, biom_range[2]), col=times_palette[4], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[4], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom4, means-sds, calib_data$biom4, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[4])
legend("topright","(d)", bty="n")
text(calib_data$biom4, means+.5, ID4, cex=0.5)

means<-colMeans(biomass_simulated_p5)
sds<-colSdApply(biomass_simulated_p5)
plot(calib_data$biom5, means, ylim=c(0, max(c(means,calib_data$biom5))), xlim=c(0, biom_range[2]), col=times_palette[5], pch=as.numeric(as.factor(round(calib_data$time3/10)*10)), main=class_names[5], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom5, means-sds, calib_data$biom5, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[5])
legend("topright","(e)", bty="n")
text(calib_data$biom5, means+.5, ID5, cex=0.5)

means<-colMeans(biomass_simulated_p6)
sds<-colSdApply(biomass_simulated_p6)
plot(calib_data$biom6, means, ylim=c(0, max(c(means,calib_data$biom6))), xlim=c(0, biom_range[2]), col=times_palette[6], pch=as.numeric(as.factor(round(calib_data$time2/10)*10)), main=class_names[6], ylab=expression(paste("Predicted, Biomass, g m"^"-2")), xlab=expression(paste("Measured, Biomass, g m"^"-2")))
#arrows(calib_data$biom6, means-sds, calib_data$biom6, means+sds, length=0.03, angle=90, code=3, col=times_palette_alpha[6])
legend("topright","(f)", bty="n")
text(calib_data$biom6, means+.5, ID6, cex=0.5)

dev.off()

png("Scatter_plot_MAE_vs_meas_means.png", width=4000, height=2500, res=390)
par(mfrow=c(2,3), mar=c(5,5,2,2))
plot(colMeans(biomass_simulated_p1),colMeans(t(dataframe_MAE_p1[,-1])),  ylab="MAE (average)", xlab="Simulated biomass", main=class_names[1])
plot( colMeans(biomass_simulated_p2),colMeans(t(dataframe_MAE_p2[,-1])), ylab="MAE (average)", xlab="Simulated biomass", main=class_names[2])
plot( colMeans(biomass_simulated_p3),colMeans(t(dataframe_MAE_p3[,-1])), ylab="MAE (average)", xlab="Simulated biomass", main=class_names[3])
plot( colMeans(biomass_simulated_p4),colMeans(t(dataframe_MAE_p4[,-1])), ylab="MAE (average)", xlab="Simulated biomass", main=class_names[4])
plot( colMeans(biomass_simulated_p5),colMeans(t(dataframe_MAE_p5[,-1])), ylab="MAE (average)", xlab="Simulated biomass", main=class_names[5])
plot( colMeans(biomass_simulated_p6),colMeans(t(dataframe_MAE_p6[,-1])), ylab="MAE (average)", xlab="Simulated biomass", main=class_names[6])
dev.off()

plot((t(dataframe_MAE_p1[,-1])), (biomass_simulated_p1), xlab="MAE (average)", ylab="Simulated biomass", main=class_names[1])


png("Boxplot_RMSE.png", width=2000, height=2500, res=350)
par(mar=c(10,5,2,2))
RMSE_tab<-cbind(RMSE_simulated_p1, RMSE_simulated_p2, RMSE_simulated_p3, RMSE_simulated_p4, RMSE_simulated_p5, RMSE_simulated_p6)
boxplot(RMSE_tab, col=times_palette_alpha, names = class_names, las=2, ylab=expression(paste("RMSE, g m"^"-2")))
dev.off()


# plotting mean absolute error per time point
dataframe_MAE_p1<-as.data.frame(t(rbind(calib_data$time1,res_simulated_p1)))
dataframe_MAE_p2<-as.data.frame(t(rbind(calib_data$time2,res_simulated_p2)))
dataframe_MAE_p3<-as.data.frame(t(rbind(calib_data$time3,res_simulated_p3)))
dataframe_MAE_p4<-as.data.frame(t(rbind(calib_data$time4,res_simulated_p4)))
dataframe_MAE_p5<-as.data.frame(t(rbind(calib_data$time5,res_simulated_p5)))
dataframe_MAE_p6<-as.data.frame(t(rbind(calib_data$time6,res_simulated_p6)))

colnames(dataframe_MAE_p1)[1]<-"time"
colnames(dataframe_MAE_p2)[1]<-"time"
colnames(dataframe_MAE_p3)[1]<-"time"
colnames(dataframe_MAE_p4)[1]<-"time"
colnames(dataframe_MAE_p5)[1]<-"time"
colnames(dataframe_MAE_p6)[1]<-"time"

dataframe_MAE_p1_mean<-(aggregate(dataframe_MAE_p1, list(dataframe_MAE_p1$time), mean))[,-c(1,2)]
dataframe_MAE_p2_mean<-(aggregate(dataframe_MAE_p2, list(dataframe_MAE_p2$time), mean))[,-c(1,2)]
dataframe_MAE_p3_mean<-(aggregate(dataframe_MAE_p3, list(dataframe_MAE_p3$time), mean))[,-c(1,2)]
dataframe_MAE_p4_mean<-(aggregate(dataframe_MAE_p4, list(dataframe_MAE_p4$time), mean))[,-c(1,2)]
dataframe_MAE_p5_mean<-(aggregate(dataframe_MAE_p5, list(dataframe_MAE_p5$time), mean))[,-c(1,2)]
dataframe_MAE_p6_mean<-(aggregate(dataframe_MAE_p6, list(dataframe_MAE_p6$time), mean))[,-c(1,2)]

png("Boxplot_MAE.png", width=4000, height=2000, res=350)
MAE_range<-range(c(dataframe_MAE_p1_mean,dataframe_MAE_p2_mean, dataframe_MAE_p3_mean, dataframe_MAE_p4_mean, dataframe_MAE_p5_mean, dataframe_MAE_p6_mean))
par(mfrow=c(2,3), mar=c(5,5,2,2))
boxplot(t(dataframe_MAE_p1_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[1], main=class_names[1], ylim=MAE_range)
boxplot(t(dataframe_MAE_p2_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[2], main=class_names[2], ylim=MAE_range)
boxplot(t(dataframe_MAE_p3_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[3], main=class_names[3], ylim=MAE_range)
boxplot(t(dataframe_MAE_p4_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[4], main=class_names[4], ylim=MAE_range)
boxplot(t(dataframe_MAE_p5_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[5], main=class_names[5], ylim=MAE_range)
boxplot(t(dataframe_MAE_p6_mean), names=unique(calib_data$time1), ylab=expression(paste("MAE, g m"^"-2")), xlab="Day of the year", col=times_palette_alpha[6], main=class_names[6], ylim=MAE_range)
dev.off()



# plotting coefficient of variation % per time point
dataframe_CV_p1<-as.data.frame(t(rbind(calib_data$time1,cv_simulated_p1)))
dataframe_CV_p2<-as.data.frame(t(rbind(calib_data$time2,cv_simulated_p2)))
dataframe_CV_p3<-as.data.frame(t(rbind(calib_data$time3,cv_simulated_p3)))
dataframe_CV_p4<-as.data.frame(t(rbind(calib_data$time4,cv_simulated_p4)))
dataframe_CV_p5<-as.data.frame(t(rbind(calib_data$time5,cv_simulated_p5)))
dataframe_CV_p6<-as.data.frame(t(rbind(calib_data$time6,cv_simulated_p6)))

colnames(dataframe_CV_p1)[1]<-"time"
colnames(dataframe_CV_p2)[1]<-"time"
colnames(dataframe_CV_p3)[1]<-"time"
colnames(dataframe_CV_p4)[1]<-"time"
colnames(dataframe_CV_p5)[1]<-"time"
colnames(dataframe_CV_p6)[1]<-"time"

dataframe_CV_p1_mean<-(aggregate(dataframe_CV_p1, list(dataframe_CV_p1$time), mean))[,-c(1,2)]
dataframe_CV_p2_mean<-(aggregate(dataframe_CV_p2, list(dataframe_CV_p2$time), mean))[,-c(1,2)]
dataframe_CV_p3_mean<-(aggregate(dataframe_CV_p3, list(dataframe_CV_p3$time), mean))[,-c(1,2)]
dataframe_CV_p4_mean<-(aggregate(dataframe_CV_p4, list(dataframe_CV_p4$time), mean))[,-c(1,2)]
dataframe_CV_p5_mean<-(aggregate(dataframe_CV_p5, list(dataframe_CV_p5$time), mean))[,-c(1,2)]
dataframe_CV_p6_mean<-(aggregate(dataframe_CV_p6, list(dataframe_CV_p6$time), mean))[,-c(1,2)]

png("Boxplot_CV.png", width=4000, height=2000, res=350)
CV_range<-range(unlist(c(dataframe_CV_p1_mean,dataframe_CV_p2_mean, dataframe_CV_p3_mean, dataframe_CV_p4_mean, dataframe_CV_p5_mean, dataframe_CV_p6_mean))*100)
par(mfrow=c(2,3), mar=c(5,5,2,2))
boxplot(t(dataframe_CV_p1_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[1], main=class_names[1], ylim=CV_range)
boxplot(t(dataframe_CV_p2_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[2], main=class_names[2], ylim=CV_range)
boxplot(t(dataframe_CV_p3_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[3], main=class_names[3], ylim=CV_range)
boxplot(t(dataframe_CV_p4_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[4], main=class_names[4], ylim=CV_range)
boxplot(t(dataframe_CV_p5_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[5], main=class_names[5], ylim=CV_range)
boxplot(t(dataframe_CV_p6_mean)*100, names=unique(calib_data$time1), ylab="Coefficient of variation, %", xlab="Day of the year", col=times_palette_alpha[6], main=class_names[6], ylim=CV_range)
dev.off()


