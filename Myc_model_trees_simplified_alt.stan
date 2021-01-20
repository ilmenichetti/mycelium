data { // this block is to define the vectors coming in as data
  int N; // N is the number of lines in the input table.
  real max_sd1; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd2; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd3; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd4; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd5; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd6; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  vector[N] biom1; //
  vector[N] t_1; //
  vector[N] t0_1; //
  vector[N] biom2; //
  vector[N] t_2; //
  vector[N] t0_2; //
  vector[N] biom3; //
  vector[N] t_3; //
  vector[N] t0_3; //
  vector[N] biom4; //
  vector[N] t_4; //
  vector[N] t0_4; //
  vector[N] biom5; //
  vector[N] t_5; //
  vector[N] t0_5; //
  vector[N] biom6; //
  vector[N] t_6; //
  vector[N] t0_6; //

} // <<<<<<<<<<<<<<< end of data block


parameters { // this block defines just the parameters to calibrate
             // here we can specify also the boundaries of them

  real<lower=0> p1; // production group 1
  real<lower=0> p2; // production group 2
  real<lower=0> p3; // production group 3
  real<lower=0> p4; // production group 4
  real<lower=0> p5; // production group 5
  real<lower=0> p6; // production group 6
  real<lower=0> mu1; //  mortality time 1
  real<lower=0> mu2; //  mortality time 3
  real<lower=0> mu3; //  mortality time 4
  real<lower=0> mu4; //  mortality time 5
  real<lower=0> mu5; //  mortality time 5
  real<lower=0> mu6; //  mortality time 5
  real<lower=0> a; // a

} // <<<<<<<<<<<<<<<  end of parameters block



model { //this block runs the model and the MCMC sampling

  //declaring the vectors used in the simulation
  vector[N] biomass1; // vector to store biomass simulated results
  vector[N] biomass2; // vector to store biomass simulated results
  vector[N] biomass3; // vector to store biomass simulated results
  vector[N] biomass4; // vector to store biomass simulated results
  vector[N] biomass5; // vector to store biomass simulated results
  vector[N] biomass6; // vector to store biomass simulated results

  real om; // om is calculated below
  om=(2*pi())/365.25; //om calculation


  // defining the parameter priors
  // we can use many distributions, but usually I decide either normal or uniform,
  p1    ~ normal(0.5,0.25);
  p2    ~ normal(0.5,0.25);
  p3    ~ normal(0.5,0.25);
  p4    ~ normal(0.5,0.25);
  p5    ~ normal(0.5,0.25);
  p6    ~ normal(0.5,0.25);
  mu1   ~ normal(0.05,0.025);
  mu2   ~ normal(0.05,0.025);
  mu3   ~ normal(0.05,0.025);
  mu4   ~ normal(0.05,0.025);
  mu5   ~ normal(0.05,0.025);
  mu6   ~ normal(0.05,0.025);
  a     ~ normal(1,0.5);

  // running the model over time
  for (i in 1:N){ // loop for each class

        biomass1[i]= (p1/mu1)*(1-exp(-mu1*(t_1[i]-t0_1[i])))+
                    ((a*p1)/mu1)*(cos(om*t_1[i])-exp(-mu1*(t_1[i]-t0_1[i]))*cos(om*t0_1[i]))-
                     ((a*p1)/(om*mu1))*(1/(1+(1/(om^2*mu1^2))))*((1/mu1)*(sin(om*t_1[i])-exp(-mu1*(t_1[i]-t0_1[i]))*sin(om*t0_1[i]))+
                     (1/(om*mu1^2))*(cos(om*t_1[i])-exp(-mu1*(t_1[i]-t0_1[i]))*cos(om*t0_1[i])));

         biom1[i] ~ normal(biomass1[i], // this line compares the simulation with the results
                           max_sd1); //the sigma is the standard deviation in the data



        biomass2[i]= (p2/mu2)*(1-exp(-mu2*(t_2[i]-t0_2[i])))+
                    ((a*p2)/mu2)*(cos(om*t_2[i])-exp(-mu2*(t_2[i]-t0_2[i]))*cos(om*t0_2[i]))-
                     ((a*p2)/(om*mu2))*(1/(1+(1/(om^2*mu2^2))))*((1/mu2)*(sin(om*t_2[i])-exp(-mu2*(t_2[i]-t0_2[i]))*sin(om*t0_2[i]))+
                     (1/(om*mu2^2))*(cos(om*t_2[i])-exp(-mu2*(t_2[i]-t0_2[i]))*cos(om*t0_2[i])));

        biom2[i] ~ normal(biomass2[i], // this line compares the simulation with the results
                           max_sd2); //the sigma is the standard deviation in the data



        biomass3[i]= (p3/mu3)*(1-exp(-mu3*(t_3[i]-t0_3[i])))+
                    ((a*p3)/mu3)*(cos(om*t_3[i])-exp(-mu3*(t_3[i]-t0_3[i]))*cos(om*t0_3[i]))-
                     ((a*p3)/(om*mu3))*(1/(1+(1/(om^2*mu3^2))))*((1/mu3)*(sin(om*t_3[i])-exp(-mu3*(t_3[i]-t0_3[i]))*sin(om*t0_3[i]))+
                     (1/(om*mu3^2))*(cos(om*t_3[i])-exp(-mu3*(t_3[i]-t0_3[i]))*cos(om*t0_3[i])));

        biom3[i] ~ normal(biomass3[i], // this line compares the simulation with the results
                           max_sd3); //the sigma is the standard deviation in the data



        biomass4[i]= (p4/mu4)*(1-exp(-mu4*(t_4[i]-t0_4[i])))+
                    ((a*p4)/mu4)*(cos(om*t_4[i])-exp(-mu4*(t_4[i]-t0_4[i]))*cos(om*t0_4[i]))-
                     ((a*p4)/(om*mu4))*(1/(1+(1/(om^2*mu4^2))))*((1/mu4)*(sin(om*t_4[i])-exp(-mu4*(t_4[i]-t0_4[i]))*sin(om*t0_4[i]))+
                     (1/(om*mu4^2))*(cos(om*t_4[i])-exp(-mu4*(t_4[i]-t0_4[i]))*cos(om*t0_4[i])));


        biom4[i] ~ normal(biomass4[i], // this line compares the simulation with the results
                           max_sd4); //the sigma is the standard deviation in the data



        biomass5[i]= (p5/mu5)*(1-exp(-mu5*(t_5[i]-t0_5[i])))+
                    ((a*p5)/mu5)*(cos(om*t_5[i])-exp(-mu5*(t_5[i]-t0_5[i]))*cos(om*t0_5[i]))-
                     ((a*p5)/(om*mu5))*(1/(1+(1/(om^2*mu5^2))))*((1/mu5)*(sin(om*t_5[i])-exp(-mu5*(t_4[i]-t0_5[i]))*sin(om*t0_5[i]))+
                     (1/(om*mu5^2))*(cos(om*t_5[i])-exp(-mu5*(t_5[i]-t0_5[i]))*cos(om*t0_5[i])));

        biom5[i] ~ normal(biomass5[i], // this line compares the simulation with the results
                           max_sd5); //the sigma is the standard deviation in the data

        biomass6[i]= (p6/mu6)*(1-exp(-mu6*(t_6[i]-t0_6[i])))+
                    ((a*p6)/mu6)*(cos(om*t_6[i])-exp(-mu6*(t_6[i]-t0_6[i]))*cos(om*t0_6[i]))-
                     ((a*p6)/(om*mu6))*(1/(1+(1/(om^2*mu6^2))))*((1/mu6)*(sin(om*t_6[i])-exp(-mu6*(t_6[i]-t0_6[i]))*sin(om*t0_6[i]))+
                     (1/(om*mu6^2))*(cos(om*t_6[i])-exp(-mu6*(t_6[i]-t0_6[i]))*cos(om*t0_6[i])));

        biom6[i] ~ normal(biomass6[i], // this line compares the simulation with the results
                           max_sd6); //the sigma is the standard deviation in the data

    } // ----------------  end of model iteration over time

} // <<<<<<<<<<<<<<<  end of model block



