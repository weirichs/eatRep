### R code from vignette source 'eatRep.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(eatRep)


###################################################
### code chunk number 2: reading_writingdef
###################################################
library(eatRep)
data(lsa)
str(lsa, give.attr = FALSE)


###################################################
### code chunk number 3: readingMeansByCountrydef
###################################################
read  <- subset(lsa, domain == "reading")
readN1<- subset(read, nest == 1 )
read10<- subset(readN1, year == 2010 )
means <- jk2.mean(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = "country",
         dependent = "score")


###################################################
### code chunk number 4: reporting01
###################################################
res01 <- report(jk2.out = means, add = list(domain = "reading"))
print(res01, digits = 4)


###################################################
### code chunk number 5: readingMeansByCountry2def
###################################################
means <- jk2.mean(datL = read10, ID = "idstud", wgt = "wgt",
         imp = "imp", groups = "country", dependent = "score")
res02 <- report(jk2.out = means, add = list(domain = "reading"))
print(res02, digits = 4)


###################################################
### code chunk number 6: readingMeansByCountry3def
###################################################
means <- jk2.mean(datL = read10, ID = "idstud", imp = "imp", groups = "country",
         dependent = "score")
res03 <- report(jk2.out = means, add = list(domain = "reading"))
print(res03, digits = 4)


###################################################
### code chunk number 7: readingMeansByCountry4def
###################################################
means <- jk2.mean(datL = subset(read10,imp==1),  ID = "idstud",
         imp = "imp", groups = "country", dependent = "score")
res04 <- report(jk2.out = means, add = list(domain = "reading"))
print(res04, digits = 4)


###################################################
### code chunk number 8: readingMeansByCountry4def
###################################################
means <- jk2.mean(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = c("sex","country"),
         group.splits = c(0,2), group.differences.by = "sex", dependent = "score")
res05 <- report(jk2.out = means, add = list(domain = "reading"))
print(res05, digits = 4)


###################################################
### code chunk number 9: read5def
###################################################
means <- jk2.mean(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = c("sex","country"),
         group.splits = 0:2, group.differences.by = "sex", cross.differences = TRUE,
         dependent = "score")
res06 <- report(jk2.out = means, add = list(domain = "reading"))


###################################################
### code chunk number 10: hiseiFreqsByCountryGenderdef
###################################################
freqs <- jk2.table( datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = c("country", "sex"),
	       group.differences.by = "sex", chiSquare = TRUE, dependent = "comp")
res07 <- report(jk2.out = freqs, add = list(domain = "reading"))
print(res07, digits = 4)


###################################################
### code chunk number 11: hiseiFreqsExtraction
###################################################
options(scipen=4)
cols <- c("group", "coefficient", "value")
frqs <- freqs[["resT"]][["noTrend"]]
res  <- frqs[which(frqs[,"parameter"] == "chiSquareTest"), cols ]
wide <- reshape2::dcast(res, group~coefficient, value.var = "value")
wide <- wide[,-grep("Approx", colnames(wide))]
print(wide, digits = 1)


###################################################
### code chunk number 12: hiseiFreqsByCountryGenderdef2
###################################################
freqs <- jk2.table( datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = c("country", "sex"),
	       group.differences.by = "sex", chiSquare = FALSE, dependent = "passReg")


###################################################
### code chunk number 13: redefineValues2def
###################################################
read10[,"passedNA"] <- read10[,"passReg"]
read10[ sample(nrow(read10), 100, FALSE) ,"passedNA"]   <- NA
freqs2<- jk2.table( datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = c("country", "sex"),
	       dependent = "passedNA", separate.missing.indicator = TRUE)


###################################################
### code chunk number 14: regression1def
###################################################
mod1  <- jk2.glm(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp", groups = "country",
         formula = score~sex*ses, family=gaussian(link="identity"), poolMethod = "scalar")


###################################################
### code chunk number 15: regression1defResults1
###################################################
res   <- report(mod1, printGlm = TRUE)


###################################################
### code chunk number 16: regression2def
###################################################
mod1  <- jk2.glm(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp",
         formula = passReg~country*sex, family=binomial(link="logit") )
res   <- report(mod1, printGlm = TRUE)


###################################################
### code chunk number 17: transformdef
###################################################
#exp(mod1[c(1,3,5,7,9),"value"])


###################################################
### code chunk number 18: recodeSexdef
###################################################
read10[,"sexRecoded"] <- factor(read10[,"sex"], levels = c("male", "female") )


###################################################
### code chunk number 19: regressionReplicationdef
###################################################
mod1  <- jk2.glm(datL = read10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", imp = "imp",
         formula = passReg~country*sexRecoded, family=binomial(link="logit") )
res   <- report(mod1, printGlm = TRUE)


###################################################
### code chunk number 20: nestEx2def
###################################################
read      <- subset(lsa, domain == "reading")
readN1.10 <- subset(read, year == 2010 )
means <- jk2.mean(datL = readN1.10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", nest="nest", imp = "imp",
         groups = "country", dependent = "score")
res   <- report(means)


###################################################
### code chunk number 21: nestEx3def
###################################################
mod1  <- jk2.glm(datL = readN1.10, ID = "idstud", wgt = "wgt", type = "JK2",
         PSU = "jkzone", repInd = "jkrep", nest="nest", imp = "imp",
         groups = "country", formula = score~sex+ses, family=gaussian(link="identity") )
res   <- report(mod1)


