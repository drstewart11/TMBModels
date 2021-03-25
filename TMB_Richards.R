#Load and install necessary packages
install.packages("TMB")
library(TMB)
library(devtools)
install_github("kaskr/TMB_contrib_R/TMBhelper")
library(TMBhelper)

tmb_model="
// Richards growth function
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

// DATA_SECTION
  DATA_VECTOR(Age); //Vector of age data
  DATA_VECTOR(Length); //Vector of length data
  int n = Age.size(); //Number of individuals

// PARAMETER_SECTION
  PARAMETER(Linf); // Asymptoptic length
  PARAMETER(k); // Slope of the inflection point (maximum relative growth rate)
  PARAMETER(b); // Mean length or vertical position of the inflection point
  PARAMETER(t0); // Inflection point of the curve
  PARAMETER(logSigma); // Fit sigma on a log scale to keep it > 0

// PRELIMINARY_CALCS_SECTION: (Transformed parameters)
  Type sigma = exp(logSigma);
  ADREPORT(sigma)
  

// PROCEDURE_SECTION
  vector<Type> La_i(n);
  La_i = Linf*pow((1-(1/b)*exp(-k*(Age-t0))),b); //Richards growth model

  Type nll = 0.0; // Initialize negative log-likelihood

  for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // negative log likelihood 
  nll -= dnorm(Length[i], La_i[i], sigma, TRUE); 
}

return nll;
}"
#Write model to C++
write(tmb_model, file = "NLL_Richards.cpp")

#Load model
compile("NLL_Richards.cpp") #Compile into a shared object file
dyn.load(dynlib("NLL_Richards"))

#Run the model
parameters<-list(Linf=650,k=0.4,b=1,t0=0,logSigma=5) #Starting values
obj = MakeADFun(data, parameters, DLL = "NLL_Richards",silent=TRUE) #Prepare model objects for fitting
#Define lower bounds to be positve for all parameters escept t0
Rich_opt = nlminb(start=obj$par, obj=obj$fn, gr=obj$gr,
                  lower=c(0,0,0,-10,1)) #Fit model using nlminb optimizer (quasi-Newton method)

#Likelihood profiles
Linf_CI_R=confint(tmbprofile(obj,"Linf",trace=FALSE))
k_CI_R=confint(tmbprofile(obj,"k",trace=FALSE))
b_CI_R=confint(tmbprofile(obj,"b",trace=FALSE))
t0_CI_R=confint(tmbprofile(obj,"t0",trace=FALSE))
logSigma_CI_R=confint(tmbprofile(obj,"logSigma",trace=FALSE))

(fit_Rich = summary(sdreport(obj))) #Summary of parameter estimates, Std. Errors, and 95% Confidence Intervals

