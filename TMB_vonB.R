#Load and install necessary packages
install.packages("TMB")
library(TMB)
library(devtools)
install_github("kaskr/TMB_contrib_R/TMBhelper")
library(TMBhelper)


#Write C++ model template
#Note: Comments follow //
tmb_model="
// von Bertalanffy growth function
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {

// DATA_SECTION
  DATA_VECTOR(Age); //Vector of age data
  DATA_VECTOR(Length); //Vector of length data
  int n = Age.size(); //Number of individuals

// PARAMETER_SECTION
  PARAMETER(Linf); // Asymptoptic length
  PARAMETER(k); // Curvature parameter
  PARAMETER(t0); // Location parameter 
  PARAMETER(logSigma); // Fit sigma on a log scale to keep it > 0

// PRELIMINARY_CALCS_SECTION: (Transformed parameters)
  Type sigma = exp(logSigma);
  ADREPORT(sigma)

// PROCEDURE_SECTION
  vector<Type> La_i(n);
  La_i = Linf*(1-exp(-k*(Age-t0))); //von Bertalanffy growth equation

  Type nll = 0.0; // Initialize negative log-likelihood

  for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // negative log likelihood
  nll -= dnorm(Length[i], La_i[i], sigma, TRUE); 
}

return nll;
}"
#Write model to C++
write(tmb_model, file = "NLL_vbgf.cpp")

#Load model
compile("NLL_vbgf.cpp") #Compile into a shared object file

dyn.load(dynlib("NLL_vbgf"))


#Parameterize and run the model
data<-list(Length=la.new$Length, Age=la.new$Age) #Define data to be tested
parameters<-list(Linf=650,k=0.4,t0=0,logSigma=5) #Starting values
obj = MakeADFun(data, parameters, DLL = "NLL_vbgf",silent=TRUE) #Prepare model objects for fitting
vbgf_opt = nlminb(start=obj$par, obj=obj$fn, gr=obj$gr) #Fit model using nlminb optimizer (quasi-Newton method)
(fit_vonB = summary(sdreport(obj))) #Summary of optimal parameter estimates and their Std. Errors

#Plot observed and predicted data
{with(la.new, plot(Length~Age, data=la.new,xlab= "Age (yr)",ylab="Total Length (mm)",las=1,bty="n",col="gray",xlim=c(0,max(la.new$Age)+1),ylim=c(0,max(la.new$Length)+100)))
  curve(fit_vonB[1,1]*(1-exp(-fit_vonB[2,1]*(x-fit_vonB[3,1]))),from=1,to=max(la.new$Age),lwd=3,add=T)}


#Likelihood profiles
Linf_prof<-tmbprofile(obj,"Linf",trace=FALSE)
k_prof<-tmbprofile(obj,"k",trace=FALSE)
t0_prof<-tmbprofile(obj,"t0",trace=FALSE)
logSigma_prof<-tmbprofile(obj,"logSigma",trace=FALSE)

#95% confidence intervals
(CILinf=round(confint(Linf_prof)))
(CIk=round(confint(k_prof),2))
(CIt0=round(confint(t0_prof),2))
(CIlogsigma=round(confint(logSigma_prof),2))

par(mfrow=c(2,2)) #2x2 plot
par(mar=c(4.2,4.1,1.1,1)) #Change margins of plot: bottom, left, top, and right margins
plot(Linf_prof,xlab=expression("L"[infinity]), ylab="Value",cex=1.75,cex.axis=1.75,cex.lab=1.75)
plot(k_prof,xlab="k",ylab="Value",cex=1.75,cex.axis=1.75,cex.lab=1.75)
plot(t0_prof,xlab=expression("t"[0]), ylab="Value",cex=1.75,cex.axis=1.75,cex.lab=1.75)
plot(logSigma_prof,ylab="Value",cex=1.75,cex.axis=1.75,cex.lab=1.75)
