#Load and install necessary packages
install.packages("TMB")
library(TMB)
library(devtools)
install_github("kaskr/TMB_contrib_R/TMBhelper")
library(TMBhelper)



tmb_model="
// Gompertz growth function
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

// DATA_SECTION
  DATA_VECTOR(Age); //Vector of age data
  DATA_VECTOR(Length); //Vector of length data
  int n = Age.size(); //Number of individuals

// PARAMETER_SECTION
  PARAMETER(Linf); // Asymptoptic length
  PARAMETER(G); // Instantaneous rate of growth
  PARAMETER(t0); // Inflection point of the curve
  PARAMETER(logSigma); // Fit sigma on a log scale to keep it > 0

// PRELIMINARY_CALCS_SECTION: (Transformed parameters)
  Type sigma = exp(logSigma);
  ADREPORT(sigma)

// PROCEDURE_SECTION
  vector<Type> La_i(n);
  La_i = Linf*exp(-exp(-G*(Age-t0))); //Gompertz growth model

  Type nll = 0.0; // Initialize negative log-likelihood

  for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // negative log likelihood
  nll -= dnorm(Length[i], La_i[i], sigma, TRUE); 
}

return nll;
}"
#Write model to C++
write(tmb_model, file = "NLL_Gompertz.cpp")


dyn.load(dynlib("NLL_Gompertz"))#Load object file

#Run the model
parameters<-list(Linf=650,G=0.4,t0=0,logSigma=5) #Starting values
obj = MakeADFun(data, parameters, DLL = "NLL_Gompertz",silent=TRUE) #Prepare model objects for fitting
Gomp_opt = nlminb(start=obj$par, obj=obj$fn, gr=obj$gr) #Fit model using nlminb optimizer (quasi-Newton method)

#Likelihood profiles
Linf_CI_G=confint(tmbprofile(obj,"Linf",trace=FALSE))
G_CI_G=confint(tmbprofile(obj,"G",trace=FALSE))
t0_CI_G=confint(tmbprofile(obj,"t0",trace=FALSE))
logSigma_CI_G=confint(tmbprofile(obj,"logSigma",trace=FALSE))

(fit_Gomp = summary(sdreport(obj))) #Summary of parameter estimates, Std. Errors, and 95% Confidence Intervals



