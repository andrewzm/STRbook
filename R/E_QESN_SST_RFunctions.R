#######################################################
################### File of functions
#######################################################
#### File reads all in R functions
#### needed for the E-QESN SST example
#######################################################
#######################################################


############################################################
############### Create Embed Data ##########################
############################################################
#' @export
createEmbedRNNData=function(m,tauEmb,tau,yTrain,rawDataInput,curXTestIndex){

  curTrainLen <- nrow(yTrain)
  numLocs <- ncol(yTrain)
  lenInSampleEmb=curTrainLen-(m*tauEmb)
  curTestLen <- length(curXTestIndex)

  altInSampleXRaw=array(NA,c(lenInSampleEmb,m+1,numLocs))

  for(i in 1:lenInSampleEmb){
    #####using laplacian eigenmaps
    altInSampleXRaw[i,,]=rawDataInput[seq(i,(m*tauEmb+i),by=tauEmb),]
  }

  #######Scale in-sample x and y
  altYInSampRaw=yTrain[(m*tauEmb+1):curTrainLen,]
  altYMean=mean(altYInSampRaw)
  altScaleFactor=sd(altYInSampRaw)

  altYInSamp=(altYInSampRaw-altYMean)/altScaleFactor


  meanXTrainMatrix=mean(rawDataInput[1:curTrainLen,])
  sdXTrainMatrix=sd(rawDataInput[1:curTrainLen,])


  altInSampleX=(altInSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix


  designMatrix=matrix(1,lenInSampleEmb,(m+1)*numLocs+1)
  for(i in 1:lenInSampleEmb){
    designMatrix[i,2:((m+1)*numLocs+1)]=as.vector(altInSampleX[i,,])
  }


  ####### Out Sample
  allOutSampIndexes=(curXTestIndex[1]-tau+1):curXTestIndex[length(curXTestIndex)]

  xMatchIndexes=match(xTestIndex,allOutSampIndexes)

  altOutSampleXRaw=array(NA,c(length(allOutSampIndexes),m+1,numLocs))
  for(i in 1:length(allOutSampIndexes)){
    altOutSampleXRaw[i,,]=rawDataInput[seq(allOutSampIndexes[i]-(m*tauEmb),allOutSampIndexes[i],by=tauEmb),]
  }

  #######Scale in-sample x and y
  altOutSampleX=(altOutSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix

  designMatrixOutSample=matrix(1,length(allOutSampIndexes),(m+1)*numLocs+1)
  for(i in 1:length(allOutSampIndexes)){
    designMatrixOutSample[i,2:((m+1)*numLocs+1)]=as.vector(altOutSampleX[i,,])
  }

  return(list(curInSampMean=altYMean,curInSampSD=altScaleFactor,curYInSamp=altYInSamp,curInSampX=altInSampleX,curOutSampX=altOutSampleX,lenInSampleEmb=lenInSampleEmb,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,curTestLen=curTestLen,xMatchIndexes=xMatchIndexes ))


}


################################################################
############### Ensembel ESN Set Parameters ####################
################################################################
#' @export
setParsEESN=function(regPar,nh,numLocs,m, quadInd){

  inSize = numLocs*(m+1)
  nColsU=inSize+1

  #########sampleVec
  sampVecESN=0:(nh-1)

  strValuesXtemp=rep(0,nh)

  #####create Ridge Matrix
  if(quadInd){
    ridgeMat=regPar*diag(2*nh )
  }else{
    ridgeMat=regPar*diag(nh)
  }
  return(list(strValuesXtemp=strValuesXtemp,sampVecESN=sampVecESN,ridgeMat=ridgeMat,nColsU=nColsU))
}



################################################################
################################################################
############### genResR ########################################
################################################################
################################################################
###### Generates a single reservoir

#' @export
genResR=function(nh,wWidth,uWidth,piW,piU,nuESN,quadInd,DataObj,setParObj,testLen){

  echoTrainLen=DataObj$lenInSampleEmb
  designMat=DataObj$designMatrix
  designMatrixOutSample=DataObj$designMatrixOutSample
  startValuesHMat=setParObj$strValuesXtemp
  ridgeMat=setParObj$ridgeMat
  scaledInSampY=DataObj$curYInSamp
  nColsU=ncol(designMat)
  inSampMean=DataObj$curInSampMean
  inSampStDev=DataObj$curInSampSD
  xMatchIndexes=DataObj$xMatchIndexes
  # DataObj,setParObj

  U=matrix(runif(nh*nColsU,-uWidth,uWidth),nh)
  W=matrix(runif(nh*nh,-wWidth,wWidth),nh)


  ######## Make unscaled W matrix sparse
  for(i in 1:nh){
    numNonZeroW=rbinom(1,nh,piW)
    numZeroW=nh-numNonZeroW

    curWIndexZero=sample(1:nh,numZeroW)
    W[curWIndexZero,i]=0
  }


  ######## Make U matrix sparse
  for(i in 1:nColsU){
    numNonZeroU=rbinom(1,nh,piU)
    numZeroU=nh-numNonZeroU

    curUIndexZero=sample(1:nh,numZeroU)
    U[curUIndexZero,i]=0
  }

  ######calculate the spectral radius
  spectralRadius = abs(eigen(W,only.values=TRUE)$values[1])

  ##### Scale the W matrix
  wMatScaled = W * nuESN / spectralRadius

  ##### Calcaulte u Mat times design Matrix
  uProdMat=U%*%t(designMat)

  ##### Calculate in-sample hMat
  srtValsInSamp=startValuesHMat


  inSampleHMatObj=createHMatR(nh,echoTrainLen,wMatScaled,uProdMat,srtValsInSamp,quadInd)

  inSampleHMat=inSampleHMatObj$hMat
  ##### Save last in-sample h Value
  hLastInSamp=inSampleHMatObj$xTemp

  ######Calculate vMat
  vMatESN=t(scaledInSampY)%*%t(inSampleHMat)%*%solve(inSampleHMat%*%t(inSampleHMat)+ridgeMat)


  ###### Calculate Out-of-Sample Forecasts
  uProdMatOutSamp=U%*%t(designMatrixOutSample)
  ##### Calculate Out-of-Sample hMat
  outSampleHMat=createHMatR(nh,nrow(designMatrixOutSample),wMatScaled,uProdMatOutSamp,hLastInSamp,quadInd)$hMat

  #####Out-Sample Forecasts
  scaledForecast=vMatESN%*%outSampleHMat
  unScaledForecasts=(scaledForecast*inSampStDev+inSampMean)[,xMatchIndexes]


  return(list(U=U,wMatScaled=wMatScaled,vMatESN=vMatESN,inSampleHMat=inSampleHMat,unScaledForecasts=unScaledForecasts))
}

################################################################
############### createHMatR ####################################
################################################################
##### Calculates hidden unit matrix hMat
#' @export
createHMatR=function(nh,numTimePerds,wMatScaled,uProdMat,startValues,quadInd){
  xTemp=startValues



  if(quadInd){
    hMat = matrix(0,2*nh ,numTimePerds)
  }else{
    hMat = matrix(0,nh ,numTimePerds)
  }



  for (t in 1:numTimePerds){

    xTemp =tanh( wMatScaled %*% xTemp+uProdMat[,t] )

    ######allows for quadratic terms
    if(quadInd){
      hMat[,t]=c(xTemp,xTemp^2)
    }else{
      hMat[,t] = c(xTemp)
    }

  }
  return(list(hMat=hMat,xTemp=xTemp))
}


#######################################################
####### Find Month for Plotting  ######################
#######################################################
#' @export
findMonth=function(curPeriod){
  numMonth=curPeriod%%12

  if(numMonth==0){
    curMonth="Dec."
  }

  if(numMonth==1){
    curMonth="Jan."
  }

  if(numMonth==2){
    curMonth="Feb."
  }

  if(numMonth==3){
    curMonth="March"
  }
  if(numMonth==4){
    curMonth="April"
  }
  if(numMonth==5){
    curMonth="May"
  }
  if(numMonth==6){
    curMonth="June"
  }
  if(numMonth==7){
    curMonth="July"
  }
  if(numMonth==8){
    curMonth="Aug."
  }
  if(numMonth==9){
    curMonth="Sept."
  }

  if(numMonth==10){
    curMonth="Oct."
  }
  if(numMonth==11){
    curMonth="Nov."
  }

  return(curMonth)
}


#######################################################
####### Create month/year label  ######################
#######################################################
#' @export
yrMonthLabel=function(curPeriod,curYr){
  numMonth=curPeriod%%12

  curMontNum=ifelse(numMonth==0,12,numMonth)

  tempLastDigit=round((curYr*.01-as.integer(curYr*.01))*100)

  returnLabel=paste(curMontNum,"/",tempLastDigit, sep = "")

  return(returnLabel)
}


#######################################################
####### Find Year for Plotting  #######################
#######################################################
#' @export
findYear=function(strYear,curPerd){
  if((curPerd%%12)==0){
    returnYr=strYear+as.integer(curPerd/12)-1
  }else{
    returnYr=strYear+as.integer(curPerd/12)
  }
  return(returnYr)
}
