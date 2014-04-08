#######################################################################################
#should_diagnostic_script.R is licensed under a Creative Commons Attribution - Non commercial 3.0 Unported License. see full license at the end of this file.
#######################################################################################
#TEMPLATE_FOR _META_ANALYSIS_OF_DIAGNOSTIC_ACCURACY#
#this script follows a combination of guidelines proposed by Doebler and Holling, according to (http://cran.r-project.org/web/packages/mada/vignettes/mada.pdf)#
#
#
#####################################################################################
#SETTING ENVIRONMENT
#####################################################################################
#Installing packages needed for the analysis
#command below will install individual and is only run once. 
#Remove the hash tag if this is the first time you are running the code, and then you can add the hash tag again
#install.packages("mada",repos="http://cran.r-project.org")
#install.packages("car", repos="http://cran.r-project.org")
#install.packages("ggplot2", repos="http://cran.r-project.org")
#install.packages("RCurl", repos="http://cran.r-project.org")
#install.packages("gdata", repos="http://cran.r-project.org")
#install.packages("meta", repos="http://cran.r-project.org")

#Loading packages
lapply(c("metafor","ggplot2","car" ,"mada","RCurl","gdata","meta"), library, character.only=T)

#Load packages (after installed) with the library function
#lapply(c("metafor","ggplot2","gridExtra" ,"psych", "RCurl", "irr","pgirmess", "nortest", "moments",
#         "GPArotation","nFactors","gdata"), library, character.only=T)
########################################################################################################
#IMPORTING DATA AND RECODING
######################################################################################################
#Data for shoulder subscapular diagnostic   
#Importing data set from the Spredsheet in google docs (Insert link)
webdata <- getURL("https://docs.google.com/spreadsheet/pub?key=0AgQJ5wHo6dDwdERNX2c2RGVtWjBiQzd0NW5teEZoMWc&single=true&gid=9&output=csv"
                  ,ssl.verifypeer = FALSE)
data<-read.csv(textConnection(webdata)) #Readin .xml data into .csv
#Recoding MRI subgroups categories
data$MRI_Arthro<-car::recode(data$MRI_Arthro,"'arthro'='Arthro';'MRI'='MRI';else='Not Specified'")
#Recoding Tesla subgroups categories
data$Tesla<-car::recode(data$Tesla,"'3T'='3T';'under 1.5T'='1.5T';else='Not Specified'")
#Organizing dataset with general data
metanalise<-with(data,data.frame(TP,FP,TN,FN,names)) #Group variables needed for the analysis
#summary(metanalise)
#Organizing dataset with MRI data Specificity
specDataMRI<-with(data,data.frame(names,Positive,TP,MRI_Arthro))
#Organizing dataset with MRI data Sensitivity
sensDataMRI<-with(data,data.frame(names,Negative,TN,MRI_Arthro))
#Organizing dataset with Tesla data Specificity
specTesla<-with(data,data.frame(names,Positive,TP,Tesla))
#Organizing dataset with Tesla data Sensitivity
sensTesla<-with(data,data.frame(names,Negative,TN,Tesla))

#Data for subscapular diagnostic SRMA - QUADAS 
webqualitydata <- getURL("https://docs.google.com/spreadsheet/pub?key=0AgQJ5wHo6dDwdERNX2c2RGVtWjBiQzd0NW5teEZoMWc&single=true&gid=16&output=csv",ssl.verifypeer = FALSE)
qualitydata<-read.csv(textConnection(webqualitydata)) 

############################################################################
#Figure 1. Quality Assessment
############################################################################
summary(qualitydata)
#attach(qualitydata)
#Generate ggplot graph fro quality data information
ggplot(qualitydata, aes(Item, Author)) + geom_tile(aes(fill = Value),
colour = "white") + scale_fill_gradient(low = "white",
high = "steelblue", name="", breaks=c(0,5,10), labels=c("No","Not Clear","Yes")) +
 theme(axis.text.x = element_text(angle 
 = 330, hjust = 0, colour = "black",size=14),
        axis.text.y = element_text(colour = "black",size=14),
        axis.title.x = element_text(face="bold",size=14),
        axis.title.y = element_text(face="bold", size=14))
###########################################################################################
#General Diagnostic Metanalysis Calculations
###########################################################################################
#Get general metanalysis estimates for sensitivity and specificity
madad(metanalise)
#Forest plots for True Negative/Total Negative
mada::forest(madad(metanalise), type = "spec",plotci=TRUE,
	snames=metanalise$names)
#Forest plots for True Positives/Total Positives
mada::forest(madad(metanalise), type = "sens",plotci=TRUE,
	snames=metanalise$names)
###########################################################################################
#Figure 2. Diagnostic metanalysis model for MRI
###########################################################################################
#Forest Plot for Sensitivity subgrouped by MRI
metaSpecMRI<-metaprop(TP,Positive,names, sm="PLN",data=specDataMRI,byvar=MRI_Arthro)
meta::forest(metaSpecMRI)
#Forest Plot for Specificity subgrouped by MRI
metaSensMRI<-metaprop(TN,Negative,names, sm="PLN",data=sensDataMRI,byvar=MRI_Arthro)
meta::forest(metaSensMRI)

###########################################################################################
#Figure 3. Metanalysis model for Sensitivity
###########################################################################################
#Get general metanalysis estimates for sensitivity and specificity
#Forest Plot for Sensitivity subgrouped by Tesla
metaSpecTesla<-metaprop(TP,Positive,names, sm="PLN",data=specTesla,byvar=Tesla)
meta::forest(metaSpecTesla)
#Forest Plot for Specificity subgrouped by Tesla
mesaSensTesla<-metaprop(TN,Negative,names, sm="PLN",data=sensTesla,byvar=Tesla)
meta::forest(mesaSensTesla)

###########################################################################################
#Figure 4: ROC Curve
###########################################################################################

ROCCurve <- function (x, extrapolate = FALSE, plotsumm = TRUE, level = 0.95, 
											ylim = c(0, 1), xlim = c(0, 1), pch = 1, sroclty = 1, sroclwd = 1, 
											predict = FALSE, predlty = 3, predlwd = 1, type = "ruttergatsonis", 
											...) 
{
	plot(c(2, 2), ylim = ylim, xlim = xlim, xlab = "Specificity", 
			 ylab = "Sensitivity", ...)
	if (length(coef(x)) == 2) {
		FP <- x$freqdata$FP
		negatives <- FP + x$freqdata$TN
		FPR <- FP/negatives
		spec <- 1-FPR
		if (extrapolate) {
			bound = c(0, 1)
		}
		if (!extrapolate) {
			bound = c(min(spec), max(spec))
		}
		srocmat <- sroc(x, type = type)
		srocmat[,1] <- 1-srocmat[,1]
		lines(srocmat[cut(srocmat[, 1], bound, "withinbound") == 
										"withinbound", ], lty = sroclty, lwd = sroclwd)
	}
	else {
		warning("Not plotting any SROC for meta-regression")
	}
	if (plotsumm) {
		ROCellipse <- ROCellipse(x, level = level, add = FALSE, pch = pch, ...)
		ROCellipse$ROCellipse[,1] <- 1-ROCellipse$ROCellipse[,1]
		
		lines(ROCellipse$ROCellipse, ...)
		sumpt <- ROCellipse$fprsens
		points(1-sumpt[1], sumpt[2], pch =pch, ...)
		
	}
	if (predict) {
		alpha.sens <- x$alphasens
		alpha.fpr <- x$alphafpr
		alpha.spec <- 1- alpha.fpr
		mu <- x$coefficients["(Intercept)", ]
		Sigma <- x$Psi + vcov(x)
		talphaellipse <- ellipse(Sigma, centre = mu, level = level)
		predellipse <- matrix(0, ncol = 2, nrow = nrow(talphaellipse))
		predellipse[, 1] <- 1 - mada:::inv.trafo(alpha.fpr, talphaellipse[, 2])
		predellipse[, 2] <- mada:::inv.trafo(alpha.sens, talphaellipse[, 1])
		lines(predellipse, lty = predlty, lwd = predlwd)
	}
	return(invisible(NULL))
}


## Example
(fit <- reitsma(metanalise))
summary(fit)
ROCCurve(fit,xlim=c(1,0.5),predict=TRUE)
points(spec(metanalise), sens(metanalise), pch = 2, col=c("darkgrey"))
legend ("bottomright", c("Raw Data","Summary Estimate"),pch = c(2,1))
#legend ("bottomleft", c("SROC", "Conf.Region"),lwd = c(2,1))

#######################################################################################
#should_diagnostic_script.R is licensed under a Creative Commons Attribution - Non commercial 3.0 Unported License. You are free: to Share — to copy, distribute and transmit the work to Remix — to adapt the work, under the following conditions: Attribution — You must attribute the work in the manner specified by the author or licensor (but not in any way that suggests that they endorse you or your use of the work). Noncommercial — You may not use this work for commercial purposes. With the understanding that: Waiver — Any of the above conditions can be waived if you get permission from the copyright holder. Public Domain — Where the work or any of its elements is in the public domain under applicable law, that status is in no way affected by the license. Other Rights — In no way are any of the following rights affected by the license: Your fair dealing or fair use rights, or other applicable copyright exceptions and limitations; The author's moral rights; Rights other persons may have either in the work itself or in how the work is used, such as publicity or privacy rights. Notice — For any reuse or distribution, you must make clear to others the license terms of this work. The best way to do this is with a link to this web page. For more details see http://creativecommons.org/licenses/by-nc/3.0/
#######################################################################################

##TALITHA, VAMOS USAR ISSO?!?!?!! ############################
#ROCellipse (primeiros 3 comandos nao estao rodando)
metanalise2<-remove.vars(metanalise,"names")
rs <- rowSums(metanalise2)
weights <- 4 * rs / max(rs)
crosshair(metanalise2, xlim = c(0,0.6), ylim = c(0.4,1),
          col = 1:11, lwd = weights)
ROCellipse(fit$ROCellipse, phc = "")
points(spec(metanalise), sens(metanalise))
#diagnostic odds ratio (DOR)
(fit.DOR.DSL<-madauni(metanalise)) ## RANDOM EFFECT DIAG METANALYSIS
summary(fit.DOR.DSL)
(fit.DOR.MH<-madauni(metanalise, method = "MH"))
summary(fit.DOR.MH)
mada::forest(fit.DOR.MH)
#proportional hazards model approach
(fit.phm.homo <- phm(metanalise, hetero =  FALSE))
(fit.phm.het <- phm(metanalise))
#summary(fit.phm.homo)
summary(fit.phm.het)
#SROC curve
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(metanalise, add = TRUE)
#bivariate approach
(fit.reitsma <- reitsma(metanalise))
summary(fit.reitsma)
plot(fit.reitsma, sroclwd = 2,
    main = "SROC curve (bivariate model) for metanalise data")
points(fpr(metanalise), sens(metanalise), pch = 2)
legend ("bottomright", c("data","summary estimate"),pch = c(2,1))
legend ("bottomleft", c("SROC", "conf.region"),lwd = c(2,1))

