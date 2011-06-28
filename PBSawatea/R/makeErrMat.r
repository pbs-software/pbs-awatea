#makeErrMat-----------------------------2011-05-05
# Make a simple ageing error matrix for Awatea.
#-----------------------------------------------RH
makeErrMat = function(N=60, ondiag=0.8, offdiag=0.1, corner=0.9) {
	errMat = diag(ondiag,N,N)
	for (i in 1:(N-1))
		errMat[i,i+1] = offdiag
	for (j in 1:(N-1))
		errMat[j+1,j] = offdiag
	errMat[1,1] = errMat[N,N] = corner
	write.table(errMat,file="errmat.dat",sep="\t",row.names=FALSE,col.names=FALSE)
	return(errMat)
}
#out=makeErrMat()