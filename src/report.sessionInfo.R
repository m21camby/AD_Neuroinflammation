# Export the table of the sessionInfo() to Word or R markdown documents
# https://github.com/jfrancoiscollin/ClinReport/blob/master/R/report.sessionInfo.R



report.sessionInfo=function()
{
	

mkLabel <- function(L, n) {
	vers <- sapply(L[[n]], function(x) x[["Version"]])
	pkg <- sapply(L[[n]], function(x) x[["Package"]])
	paste(pkg, vers, sep = "_")
}

s=sessionInfo()



version=s$R.version$version.string
platform=s$platform
osversion=s$running
locale=strsplit(Sys.getlocale(),";")[[1]]
matrix=paste("Matrix products: ", s$matprod, "\n", sep = "")
basepkg=s$basePkgs
otherpkg=mkLabel(s, "otherPkgs")
loaded=mkLabel(s, "loadedOnly")

lab_r="R version:"
lab_platform="Platform:"
lab_osversion="Running under:"
lab_loc=rep("Locale:",length(locale))
lab_matrix="Matrix:"
lab_attach=rep("Attached base packages:",length(basepkg))
lab_other=rep("Other attached packages:",length(otherpkg))
lab_loaded=rep("Loaded via a namespace (and not attached):",length(loaded))



val=c(version,platform,
		osversion,locale,matrix,
		basepkg,
		otherpkg,
		loaded)

var=c(lab_r,
		lab_platform,
		lab_osversion,
		lab_loc,
		lab_matrix,
		lab_attach,
		lab_other,
		lab_loaded)



d=data.frame(Label=var,
		Information=val)



return(d)


}


