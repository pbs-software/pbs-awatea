#tabSAR---------------------------------2011-07-18
# Generate comma-delimited, two-dimensional output tables from reference point objects.
#  models - names of binary system files that store the decision tables.
#  pnam   - name of list object containing matrices of reference probabilities.
#  tnam   - names of matrices reporting times to reach reference points/criteria.
#  cats   - catch strategies (subset) to report in output tables.
#  digits - number of digits to retain after the decimal.
#-----------------------------------------------RH
tabSAR = function(models=paste("input-ymr",pad0(c(29,30),2),pad0(1,2),sep="."),
    #prefix="input-ymr", run=c(29,30), rwt=1,
    pnam = "refProbs3Gen90", tnam=c("Ttab0.5", "Ttab0.8", "Ttab0.95"),
    cats = seq(0,2500,500), digits=2 ) {

	#models = paste(prefix,pad0(run,2),pad0(rwt,2),sep=".")
	files  = paste(models,"Tables.RData",sep="")
	nfiles = length(files)
	for (i in 1:nfiles) {
		ifile = files[i]
		load(ifile)

		pcsv = gsub("Tables\\.RData","_prob.csv",ifile)  # output CSV name for probabilities
		cat(models[i],"\n",file=pcsv)
		cat("Annual catch,,,,Projection Year,,,\n",file=pcsv,append=TRUE)
		probs = get(pnam)
		cat(paste(c("strategy",dimnames(probs[[1]])[[2]]),collapse=","),"\n",file=pcsv,append=TRUE)
		for (j in names(probs)) {
			cat(paste("P(Bt > ",j,")",sep=""),"\n",file=pcsv,append=TRUE)
			ptab = probs[[j]]
			ptab = ptab[as.character(cats),]
			mess = paste(paste(dimnames(ptab)[[1]],apply(ptab,1,function(x){
				paste(show0(round(x,digits),digits,add2int=TRUE),collapse=",")}),sep=","),collapse="\n") # flatten table
			cat(mess,"\n",file=pcsv,append=TRUE)
		}

		tcsv = gsub("Tables\\.RData","_targ.csv",ifile)  # output CSV name for years to target
		cat(models[i],"\n",file=tcsv)
		cat("Annual catch,,,,Target Reference,,\n",file=tcsv,append=TRUE)
		for (k in tnam) {
			ttab = get(k)
			if (k==tnam[1])
				cat(paste(c("strategy",dimnames(ttab)[[2]]),collapse=","),"\n",file=tcsv,append=TRUE)
#browser();return()
			cat(paste(as.numeric(substring(k,5))*100,"% confidence",sep=""),"\n",file=tcsv,append=TRUE)
			ttab = ttab[as.character(cats),]
			mess = paste(paste(dimnames(ttab)[[1]],apply(ttab,1,function(x){
				paste(show0(round(x,digits),digits,add2int=TRUE),collapse=",")}),sep=","),collapse="\n") # flatten table
			cat(mess,"\n",file=tcsv,append=TRUE)
		}
		
	}
}
#-------------------------------------------tabSAR
#tabSAR(digits=2)
