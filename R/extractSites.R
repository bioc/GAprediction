# Extract CpGs needed for prediction of GA (GAprediction package)
# type="se" (default) designates model with penalty term
# lambda within less one standard error of minimum
# (se="min", gives CpGs for minimum lambda model, se="all" gives all CpGs
# require by prediction model, even those not estimated)
# Returns a vector with CpG sites
extractSites<-function( type="se" ){
#  data(glmnetPredictor)
  tempMat<-NULL
  allCpGs<-NULL
  all<-FALSE
  if (type=="se"){
    # lambda's within one std. error of minimum retains fewer components
    # and model performance is comparable to minimum lambda, therefore it is default
    tempMat<-as.matrix(coef(UL.mod.cv, s="lambda.1se"))
  }
  # Minimum lambdas, slightly better performance, but more CpGs retained
  else if (type=="min") {
    tempMat<-as.matrix(coef(UL.mod.cv, s="lambda.min"))
  }
  else if (type=="all") {
    tempMat=as.matrix(coef(UL.mod.cv))
    all<-TRUE
  }
  else{
    message("Unknown type, please choose \"se\", \"min\" or \"all\" for the appropriate set of CpGs\n")
    stop("Exiting...")
  }
  # Extract all CpGs from trained model
  tempMat=data.frame(tempMat)
  names(tempMat)="est"
  # Extract only coeffcients that are used
  if(all==TRUE){
    allCpGs=rownames(tempMat)
    return(allCpGs[-1])
  }
  else{
    predictorSites=rownames(tempMat)[which(tempMat$est!=0)]
    return(predictorSites[-1])
  }
}
