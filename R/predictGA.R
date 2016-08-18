# Gestational age (GA) predictor
# Takes a methylation matrix with CpG names and delivers it to the pre-computed
# prediction model and returns a data frame with predicted GA and sample IDs. The model is based
# on ultrasound GA estimaes which where found to produce less noisy statistical models.

predictGA<-function( mldat, transp=TRUE, se=TRUE ){
# Will be needing the glmnet package
		temp<-dim(mldat)
# Make sure mldat a matrix
		mldat<-as.matrix(mldat)
# If the there are more rows than columns take the transpose
		if( temp[1]>temp[2] && transp==TRUE ){
			message("Transpose...\n")
			mldat<-t(mldat)
		}
# Get glmnet predictor
#		data("glmnetPredictor")
# Extract sites from predictor
		if (se==TRUE){
# lambda's within one std. error of minimum retains fewer components
# and model performance is comparable to minimum lambda, therfore it is default
		  predictorSites<-extractSites(type="se")
		}
# Minimum lambdas, slightly better performance, but more CpGs retained
		else{
		  predictorSites<-extractSites(type="min")
		}
# Extract all CpGs from trained model
		allCpGs<-extractSites(type="all")
		temp<-which(!predictorSites %in% colnames(mldat))
		if(length(temp)>0){
		  stop("Missing essential CpGs needed for prediction\n", predictorSites[temp],"\n")
		}
		mldat<-mldat[,predictorSites]
# See if there are any missing values in predictor CpGs
		for(i in predictorSites){
		  temp<-length(which(!is.finite(mldat[,i])))
		  if(temp>0){
		    stop("Missing values detected in essential CpG\n",i,"\n")
		  }
		}
		numSamples<-dim(mldat)[1]
		numCpGs<-length(allCpGs)
# Create a full matrix that correspond to the needs of the prediction model
# This might require quite a bit of memory...
		predMatr<-matrix(NA, ncol=numCpGs, nrow=numSamples)
		colnames(predMatr)<-allCpGs
# Place predictor CpGs from original matrix into new matrix
		predMatr[,predictorSites]<-mldat[,predictorSites]
		message("Predicting GA...\n")
		if( se==TRUE ){
		  predGA<-as.vector(predict(UL.mod.cv, newx=predMatr, s="lambda.1se"))
		}
		else{
		  predGA<-as.vector(predict(UL.mod.cv, newx=predMatr, s="lambda.min"))
		  message("Using minimum lambda model\n")
		}
# Make a dataframe with row-names in provided methylation matrix as sample IDs
# and a column GA with GA predictions
		predData<-data.frame(GA=predGA)
		rownames(predData)<-rownames(as.data.frame(mldat))
# Clean up the big stuff, leave the rest for the garbage collector
# Returns a data frame with GA with corresponding rownames taken from the
# provided matrix
		return(predData)
}
