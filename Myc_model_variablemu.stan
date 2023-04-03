data { // this block is to define the vectors coming in as data
  int N; // N is the number of lines in the input table.
  real max_sd1; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd2; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd3; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd4; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd5; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  real max_sd6; // max standard deviation of the subset of the dataset (assumed to be the error of all measurements)
  vector[N] biom1; //
  vector[N] tree1;//
  vector[N] time1; //
  vector[N] start1; //
  vector[N] period1; //
  vector[N] biom2; //
  vector[N] tree2;//
  vector[N] time2; //
  vector[N] start2; //
  vector[N] period2; //
  vector[N] biom3; //
  vector[N] tree3;//
  vector[N] time3; //
  vector[N] period3; //
  vector[N] start3; //
  vector[N] biom4; //
  vector[N] tree4;//
  vector[N] time4; //
  vector[N] start4; //
  vector[N] period4; //
  vector[N] biom5; //
  vector[N] tree5;//
  vector[N] time5; //
  vector[N] start5; //
  vector[N] period5; //
  vector[N] biom6; //
  vector[N] tree6;//
  vector[N] time6; //
  vector[N] start6; //
  vector[N] period6; //

} // <<<<<<<<<<<<<<< end of data block


parameters { // this block defines just the parameters to calibrate
             // here we can specify also the boundaries of them

  real<lower=0> p1; // production group 1
  real<lower=0> p2; // production group 2
  real<lower=0> p3; // production group 3
  real<lower=0> p4; // production group 4
  real<lower=0> p5; // production group 4
  real<lower=0> p6; // production group 4
  real<lower=0, upper=1> mu1; //  mortality
  real<lower=0, upper=1> mu2; //  mortality
  real<lower=0, upper=1> mu3; //  mortality
  real<lower=0, upper=1> mu4; //  mortality
  real<lower=0> z; // proportionality coefficient of tree growth
  real<lower=0, upper=1> nu1; //  mortality
  real<lower=0, upper=1> nu2; //  mortality
  real<lower=0, upper=1> nu3; //  mortality
  real<lower=0, upper=1> nu4; //  mortality
  real<lower=0, upper=1> nu5; //  mortality
  real<lower=0, upper=1> nu6; //  mortality

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
  real a1; // this creates the real to store the result of a calculation
  real a2; // this creates the real to store the result of a calculation
  real a3; // this creates the real to store the result of a calculation
  real a4; // this creates the real to store the result of a calculation
  real a5; // this creates the real to store the result of a calculation
  real a6; // this creates the real to store the result of a calculation
  om=(2*pi())/365.25; //om calculation


  // defining the parameter priors
  // we can use many distributions, but usually I decide either normal or uniform,
  p1    ~ normal(0.25,0.5);
  p2    ~ normal(0.25,0.5);
  p3    ~ normal(0.25,0.5);
  p4    ~ normal(0.25,0.5);
  p5    ~ normal(0.25,0.5);
  p6    ~ normal(0.25,0.5);
  mu1    ~ uniform(0,1);
  mu2    ~ uniform(0,1);
  mu3    ~ uniform(0,1);
  mu4    ~ uniform(0,1);
  z     ~ uniform(0,100);
  nu1    ~ normal(0.8,1.2);
  nu2    ~ normal(0.8,1.2);
  nu3    ~ normal(0.8,1.2);
  nu4    ~ normal(0.8,1.2);
  nu5    ~ normal(0.8,1.2);
  nu6    ~ normal(0.8,1.2);

  // running the model over time
  for (i in 1:N){ // loop for each class

        a1= z*tree1[i]; // calculating a
        if(period1[i]==1)
        biomass1[i]= p1/(mu1*nu1)*(1-exp(-(mu1*nu1)*(time1[i]-start1[i])))+
                    ((a1*p1)/(mu1*nu1))*(cos(om)*time1[i]-exp(-(mu1*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i])-
                     ((a1*p1)/(om*(mu1*nu1)))*(1/(1+(1/(om^2*(mu1*nu1)^2))))*((1/(mu1*nu1))*(sin(om)*time1[i]-exp(-(mu1*nu1)*(time1[i]-start1[i]))*sin(om)*start1[i])+
                     (1/(om*(mu1*nu1)^2))*(cos(om)*time1[i]-exp(-(mu1*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i]));
        else if(period1[i]==2)
        biomass1[i]= p1/(mu2*nu1)*(1-exp(-(mu2*nu1)*(time1[i]-start1[i])))+
                    ((a1*p1)/(mu2*nu1))*(cos(om)*time1[i]-exp(-(mu2*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i])-
                     ((a1*p1)/(om*(mu2*nu1)))*(1/(1+(1/(om^2*(mu2*nu1)^2))))*((1/(mu2*nu1))*(sin(om)*time1[i]-exp(-(mu2*nu1)*(time1[i]-start1[i]))*sin(om)*start1[i])+
                     (1/(om*(mu2*nu1)^2))*(cos(om)*time1[i]-exp(-(mu2*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i]));
        else if(period1[i]==3)
            biomass1[i]= p1/(mu3*nu1)*(1-exp(-(mu3*nu1)*(time1[i]-start1[i])))+
                        ((a1*p1)/(mu3*nu1))*(cos(om)*time1[i]-exp(-(mu3*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i])-
                         ((a1*p1)/(om*(mu3*nu1)))*(1/(1+(1/(om^2*(mu3*nu1)^2))))*((1/(mu3*nu1))*(sin(om)*time1[i]-exp(-(mu3*nu1)*(time1[i]-start1[i]))*sin(om)*start1[i])+
                         (1/(om*(mu3*nu1)^2))*(cos(om)*time1[i]-exp(-(mu3*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i]));
        else if(period1[i]==4)
            biomass1[i]= p1/(mu4*nu1)*(1-exp(-(mu4*nu1)*(time1[i]-start1[i])))+
                        ((a1*p1)/(mu4*nu1))*(cos(om)*time1[i]-exp(-(mu4*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i])-
                         ((a1*p1)/(om*(mu4*nu1)))*(1/(1+(1/(om^2*(mu4*nu1)^2))))*((1/(mu4*nu1))*(sin(om)*time1[i]-exp(-(mu4*nu1)*(time1[i]-start1[i]))*sin(om)*start1[i])+
                         (1/(om*(mu4*nu1)^2))*(cos(om)*time1[i]-exp(-(mu4*nu1)*(time1[i]-start1[i]))*cos(om)*start1[i]));
         biom1[i] ~ normal(biomass1[i], // this line compares the simulation with the results
                           max_sd1); //the sigma is the standard deviation in the data



        a2= z*tree2[i]; // calculating a
        if(period2[i]==1)
        biomass2[i]= p2/(mu1*nu2)*(1-exp(-(mu1*nu2)*(time2[i]-start2[i])))+
                    ((a2*p2)/(mu1*nu2))*(cos(om)*time2[i]-exp(-(mu1*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i])-
                     ((a2*p2)/(om*(mu1*nu2)))*(1/(1+(1/(om^2*(mu1*nu2)^2))))*((1/(mu1*nu2))*(sin(om)*time2[i]-exp(-(mu1*nu2)*(time2[i]-start2[i]))*sin(om)*start2[i])+
                     (1/(om*(mu1*nu2)^2))*(cos(om)*time2[i]-exp(-(mu1*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i]));
        else if(period2[i]==2)
        biomass2[i]= p2/(mu2*nu2)*(1-exp(-(mu2*nu2)*(time2[i]-start2[i])))+
                    ((a2*p2)/(mu2*nu2))*(cos(om)*time2[i]-exp(-(mu2*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i])-
                     ((a2*p2)/(om*(mu2*nu2)))*(1/(1+(1/(om^2*(mu2*nu2)^2))))*((1/(mu2*nu2))*(sin(om)*time2[i]-exp(-(mu2*nu2)*(time2[i]-start2[i]))*sin(om)*start2[i])+
                     (1/(om*(mu2*nu2)^2))*(cos(om)*time2[i]-exp(-(mu2*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i]));
        else if(period2[i]==3)
        biomass2[i]= p2/(mu3*nu2)*(1-exp(-(mu3*nu2)*(time2[i]-start2[i])))+
                    ((a2*p2)/(mu3*nu2))*(cos(om)*time2[i]-exp(-(mu3*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i])-
                     ((a2*p2)/(om*(mu3*nu2)))*(1/(1+(1/(om^2*(mu3*nu2)^2))))*((1/(mu3*nu2))*(sin(om)*time2[i]-exp(-(mu3*nu2)*(time2[i]-start2[i]))*sin(om)*start2[i])+
                     (1/(om*(mu3*nu2)^2))*(cos(om)*time2[i]-exp(-(mu3*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i]));
        else if(period2[i]==4)
        biomass2[i]= p2/(mu4*nu2)*(1-exp(-(mu4*nu2)*(time2[i]-start2[i])))+
                    ((a2*p2)/(mu4*nu2))*(cos(om)*time2[i]-exp(-(mu4*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i])-
                     ((a2*p2)/(om*(mu4*nu2)))*(1/(1+(1/(om^2*(mu4*nu2)^2))))*((1/(mu4*nu2))*(sin(om)*time2[i]-exp(-(mu4*nu2)*(time2[i]-start2[i]))*sin(om)*start2[i])+
                     (1/(om*(mu4*nu2)^2))*(cos(om)*time2[i]-exp(-(mu4*nu2)*(time2[i]-start2[i]))*cos(om)*start2[i]));

        biom2[i] ~ normal(biomass2[i], // this line compares the simulation with the results
                           max_sd2); //the sigma is the standard deviation in the data



        a3= z*tree3[i]; // calculating a
        if(period3[i]==1)
        biomass3[i]= p3/(mu1*nu3)*(1-exp(-(mu1*nu3)*(time3[i]-start3[i])))+
                    ((a3*p3)/(mu1*nu3))*(cos(om)*time3[i]-exp(-(mu1*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i])-
                     ((a3*p3)/(om*(mu1*nu3)))*(1/(1+(1/(om^2*(mu1*nu3)^2))))*((1/(mu1*nu3))*(sin(om)*time3[i]-exp(-(mu1*nu3)*(time3[i]-start3[i]))*sin(om)*start3[i])+
                     (1/(om*(mu1*nu3)^2))*(cos(om)*time3[i]-exp(-(mu1*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i]));
        else if(period3[i]==2)
        biomass3[i]= p3/(mu2*nu3)*(1-exp(-(mu2*nu3)*(time3[i]-start3[i])))+
                    ((a3*p3)/(mu2*nu3))*(cos(om)*time3[i]-exp(-(mu2*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i])-
                     ((a3*p3)/(om*(mu2*nu3)))*(1/(1+(1/(om^2*(mu2*nu3)^2))))*((1/(mu2*nu3))*(sin(om)*time3[i]-exp(-(mu2*nu3)*(time3[i]-start3[i]))*sin(om)*start3[i])+
                     (1/(om*(mu2*nu3)^2))*(cos(om)*time3[i]-exp(-(mu2*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i]));
        else if(period3[i]==3)
        biomass3[i]= p3/(mu3*nu3)*(1-exp(-(mu3*nu3)*(time3[i]-start3[i])))+
                    ((a3*p3)/(mu3*nu3))*(cos(om)*time3[i]-exp(-(mu3*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i])-
                     ((a3*p3)/(om*(mu3*nu3)))*(1/(1+(1/(om^2*(mu3*nu3)^2))))*((1/(mu3*nu3))*(sin(om)*time3[i]-exp(-(mu3*nu3)*(time3[i]-start3[i]))*sin(om)*start3[i])+
                     (1/(om*(mu3*nu3)^2))*(cos(om)*time3[i]-exp(-(mu3*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i]));
        else if(period3[i]==4)
        biomass3[i]= p3/(mu4*nu3)*(1-exp(-(mu4*nu3)*(time3[i]-start3[i])))+
                    ((a3*p3)/(mu4*nu3))*(cos(om)*time3[i]-exp(-(mu4*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i])-
                     ((a3*p3)/(om*(mu4*nu3)))*(1/(1+(1/(om^2*(mu4*nu3)^2))))*((1/(mu4*nu3))*(sin(om)*time3[i]-exp(-(mu4*nu3)*(time3[i]-start3[i]))*sin(om)*start3[i])+
                     (1/(om*(mu4*nu3)^2))*(cos(om)*time3[i]-exp(-(mu4*nu3)*(time3[i]-start3[i]))*cos(om)*start3[i]));

        biom3[i] ~ normal(biomass3[i], // this line compares the simulation with the results
                           max_sd3); //the sigma is the standard deviation in the data



        a4= z*tree4[i]; // calculating a
        if(period4[i]==1)
        biomass4[i]= p4/(mu1*nu4)*(1-exp(-(mu1*nu4)*(time4[i]-start4[i])))+
                    ((a4*p4)/(mu1*nu4))*(cos(om)*time4[i]-exp(-(mu1*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i])-
                     ((a4*p4)/(om*(mu1*nu4)))*(1/(1+(1/(om^2*(mu1*nu4)^2))))*((1/(mu1*nu4))*(sin(om)*time4[i]-exp(-(mu1*nu4)*(time4[i]-start4[i]))*sin(om)*start4[i])+
                     (1/(om*(mu1*nu4)^2))*(cos(om)*time4[i]-exp(-(mu1*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i]));
        else if(period4[i]==2)
        biomass4[i]= p4/(mu2*nu4)*(1-exp(-(mu2*nu4)*(time4[i]-start4[i])))+
                    ((a4*p4)/(mu2*nu4))*(cos(om)*time4[i]-exp(-(mu2*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i])-
                     ((a4*p4)/(om*(mu2*nu4)))*(1/(1+(1/(om^2*(mu2*nu4)^2))))*((1/(mu2*nu4))*(sin(om)*time4[i]-exp(-(mu2*nu4)*(time4[i]-start4[i]))*sin(om)*start4[i])+
                     (1/(om*(mu2*nu4)^2))*(cos(om)*time4[i]-exp(-(mu2*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i]));
        else if(period4[i]==3)
        biomass4[i]= p4/(mu3*nu4)*(1-exp(-(mu3*nu4)*(time4[i]-start4[i])))+
                    ((a4*p4)/(mu3*nu4))*(cos(om)*time4[i]-exp(-(mu3*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i])-
                     ((a4*p4)/(om*(mu3*nu4)))*(1/(1+(1/(om^2*(mu3*nu4)^2))))*((1/(mu3*nu4))*(sin(om)*time4[i]-exp(-(mu3*nu4)*(time4[i]-start4[i]))*sin(om)*start4[i])+
                     (1/(om*(mu3*nu4)^2))*(cos(om)*time4[i]-exp(-(mu3*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i]));
        else if(period4[i]==4)
        biomass4[i]= p4/(mu4*nu4)*(1-exp(-(mu4*nu4)*(time4[i]-start4[i])))+
                    ((a4*p4)/(mu4*nu4))*(cos(om)*time4[i]-exp(-(mu4*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i])-
                     ((a4*p4)/(om*(mu4*nu4)))*(1/(1+(1/(om^2*(mu4*nu4)^2))))*((1/(mu4*nu4))*(sin(om)*time4[i]-exp(-(mu4*nu4)*(time4[i]-start4[i]))*sin(om)*start4[i])+
                     (1/(om*(mu4*nu4)^2))*(cos(om)*time4[i]-exp(-(mu4*nu4)*(time4[i]-start4[i]))*cos(om)*start4[i]));

        biom4[i] ~ normal(biomass4[i], // this line compares the simulation with the results
                           max_sd4); //the sigma is the standard deviation in the data



        a5= z*tree5[i]; // calculating a
        if(period5[i]==1)
        biomass5[i]= p5/(mu1*nu5)*(1-exp(-(mu1*nu5)*(time5[i]-start5[i])))+
                    ((a5*p5)/(mu1*nu5))*(cos(om)*time5[i]-exp(-(mu1*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i])-
                     ((a5*p5)/(om*(mu1*nu5)))*(1/(1+(1/(om^2*(mu1*nu5)^2))))*((1/(mu1*nu5))*(sin(om)*time5[i]-exp(-(mu1*nu5)*(time4[i]-start5[i]))*sin(om)*start5[i])+
                     (1/(om*(mu1*nu5)^2))*(cos(om)*time5[i]-exp(-(mu1*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i]));
        else if(period5[i]==2)
        biomass5[i]= p5/(mu2*nu5)*(1-exp(-(mu2*nu5)*(time5[i]-start5[i])))+
                    ((a5*p5)/(mu2*nu5))*(cos(om)*time5[i]-exp(-(mu2*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i])-
                     ((a5*p5)/(om*(mu2*nu5)))*(1/(1+(1/(om^2*(mu2*nu5)^2))))*((1/(mu2*nu5))*(sin(om)*time5[i]-exp(-(mu2*nu5)*(time4[i]-start5[i]))*sin(om)*start5[i])+
                     (1/(om*(mu2*nu5)^2))*(cos(om)*time5[i]-exp(-(mu2*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i]));
        else if(period5[i]==3)
        biomass5[i]= p5/(mu3*nu5)*(1-exp(-(mu3*nu5)*(time5[i]-start5[i])))+
                    ((a5*p5)/(mu3*nu5))*(cos(om)*time5[i]-exp(-(mu3*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i])-
                     ((a5*p5)/(om*(mu3*nu5)))*(1/(1+(1/(om^2*(mu3*nu5)^2))))*((1/(mu3*nu5))*(sin(om)*time5[i]-exp(-(mu3*nu5)*(time4[i]-start5[i]))*sin(om)*start5[i])+
                     (1/(om*(mu3*nu5)^2))*(cos(om)*time5[i]-exp(-(mu3*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i]));
        else if(period5[i]==4)
        biomass5[i]= p5/(mu4*nu5)*(1-exp(-(mu4*nu5)*(time5[i]-start5[i])))+
                    ((a5*p5)/(mu4*nu5))*(cos(om)*time5[i]-exp(-(mu4*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i])-
                     ((a5*p5)/(om*(mu4*nu5)))*(1/(1+(1/(om^2*(mu4*nu5)^2))))*((1/(mu4*nu5))*(sin(om)*time5[i]-exp(-(mu4*nu5)*(time4[i]-start5[i]))*sin(om)*start5[i])+
                     (1/(om*(mu4*nu5)^2))*(cos(om)*time5[i]-exp(-(mu4*nu5)*(time5[i]-start5[i]))*cos(om)*start5[i]));

        biom5[i] ~ normal(biomass5[i], // this line compares the simulation with the results
                           max_sd5); //the sigma is the standard deviation in the data

        a6= z*tree6[i]; // calculating a
        if(period6[i]==1)
        biomass6[i]= p6/(mu1*nu6)*(1-exp(-(mu1*nu6)*(time6[i]-start6[i])))+
                    ((a6*p6)/(mu1*nu6))*(cos(om)*time6[i]-exp(-(mu1*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i])-
                     ((a6*p6)/(om*(mu1*nu6)))*(1/(1+(1/(om^2*(mu1*nu6)^2))))*((1/(mu1*nu6))*(sin(om)*time6[i]-exp(-(mu1*nu6)*(time6[i]-start6[i]))*sin(om)*start6[i])+
                     (1/(om*(mu1*nu6)^2))*(cos(om)*time6[i]-exp(-(mu1*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i]));
        else if(period6[i]==2)
        biomass6[i]= p6/(mu2*nu6)*(1-exp(-(mu2*nu6)*(time6[i]-start6[i])))+
                    ((a6*p6)/(mu2*nu6))*(cos(om)*time6[i]-exp(-(mu2*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i])-
                     ((a6*p6)/(om*(mu2*nu6)))*(1/(1+(1/(om^2*(mu2*nu6)^2))))*((1/(mu2*nu6))*(sin(om)*time6[i]-exp(-(mu2*nu6)*(time6[i]-start6[i]))*sin(om)*start6[i])+
                     (1/(om*(mu2*nu6)^2))*(cos(om)*time6[i]-exp(-(mu2*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i]));
        else if(period6[i]==3)
        biomass6[i]= p6/(mu3*nu6)*(1-exp(-(mu3*nu6)*(time6[i]-start6[i])))+
                    ((a6*p6)/(mu3*nu6))*(cos(om)*time6[i]-exp(-(mu3*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i])-
                     ((a6*p6)/(om*(mu3*nu6)))*(1/(1+(1/(om^2*(mu3*nu6)^2))))*((1/(mu3*nu6))*(sin(om)*time6[i]-exp(-(mu3*nu6)*(time6[i]-start6[i]))*sin(om)*start6[i])+
                     (1/(om*(mu3*nu6)^2))*(cos(om)*time6[i]-exp(-(mu3*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i]));
        else if(period6[i]==4)
        biomass6[i]= p6/(mu4*nu6)*(1-exp(-(mu4*nu6)*(time6[i]-start6[i])))+
                    ((a6*p6)/(mu4*nu6))*(cos(om)*time6[i]-exp(-(mu4*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i])-
                     ((a6*p6)/(om*(mu4*nu6)))*(1/(1+(1/(om^2*(mu4*nu6)^2))))*((1/(mu4*nu6))*(sin(om)*time6[i]-exp(-(mu4*nu6)*(time6[i]-start6[i]))*sin(om)*start6[i])+
                     (1/(om*(mu4*nu6)^2))*(cos(om)*time6[i]-exp(-(mu4*nu6)*(time6[i]-start6[i]))*cos(om)*start6[i]));
        biom6[i] ~ normal(biomass6[i], // this line compares the simulation with the results
                           max_sd6); //the sigma is the standard deviation in the data

    } // ----------------  end of model iteration over time

} // <<<<<<<<<<<<<<<  end of model block



