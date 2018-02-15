#' @export
Medwind_BHM <- function(GibbsInput,Mpre){

  #*********************************
  #  Gibbs Sampler Parameters
  #*********************************
  ngibbs <- GibbsInput$ngibbs
  nburn <- GibbsInput$nburn
  nreal <- GibbsInput$nreal

  #*********************************
  #  assign data from preprocessor
  #*********************************
  EUdat <- Mpre$Mdata$EUdat
  EVdat <- Mpre$Mdata$EVdat
  SUdat <- Mpre$Mdata$SUdat
  SVdat <- Mpre$Mdata$SVdat
  DxP <- Mpre$Mdata$DxP
  DyP <- Mpre$Mdata$DyP
  T <- Mpre$Mdata$T

  #***************************************************
  #  assign data incidence matrices from preprocessor
  #***************************************************
  He <-   Mpre$Mgrid$He
  HeHe <- Mpre$Mgrid$HeHe
  Hs <-   Mpre$Mgrid$Hs
  HsHs <- Mpre$Mgrid$HsHs
  ng <- Mpre$Mgrid$ng

  #***************************************************
  #  assign prior hyperparameters from preprocessor
  #***************************************************
  mu_uu <- Mpre$Mpriors$mu_uu
  mu_vv <- Mpre$Mpriors$mu_vv
  mu_uv <- Mpre$Mpriors$mu_uv
  mu_vu <- Mpre$Mpriors$mu_vu
  mu_up <- Mpre$Mpriors$mu_up
  mu_vp <- Mpre$Mpriors$mu_vp
  s2uu <-  Mpre$Mpriors$s2uu
  s2vv <-  Mpre$Mpriors$s2vv
  s2uv <-  Mpre$Mpriors$s2uv
  s2vu <-  Mpre$Mpriors$s2vu
  s2up <-  Mpre$Mpriors$s2up
  s2vp <-  Mpre$Mpriors$s2vp

  qu <- Mpre$Mpriors$qu
  ru <- Mpre$Mpriors$ru
  qv <- Mpre$Mpriors$qv
  rv <- Mpre$Mpriors$rv

  mu_u1 <- Mpre$Mpriors$mu_u1
  mu_v1 <- Mpre$Mpriors$mu_v1
  s2u1 <-  Mpre$Mpriors$s2u1
  s2v1 <-  Mpre$Mpriors$s2v1

  s2e <-  Mpre$Mpriors$s2e
  s2s <-  Mpre$Mpriors$s2s

  #***************************************************
  #  assign starting values from preprocessor
  #***************************************************
  u = Mpre$Mstrt$u
  v = Mpre$Mstrt$v
  theta_uu = Mpre$Mstrt$theta_uu
  theta_vv = Mpre$Mstrt$theta_vv
  theta_uv = Mpre$Mstrt$theta_uv
  theta_vu = Mpre$Mstrt$theta_vu
  theta_up = Mpre$Mstrt$theta_up
  theta_vp = Mpre$Mstrt$theta_vp
  s2u = Mpre$Mstrt$s2u
  s2v = Mpre$Mstrt$s2v

  #***********************************
  #** Parameters for Gibbs Sampler
  #***********************************
  real_int <- (ngibbs-nburn)/nreal   #interval between realizations to save (keeping nreal)
  rcnt <- 0                         #counter for realizations


  #******************************************************
  #  make sparse identity matrix of size ng (pred grid)
  #******************************************************
  Isp <- Matrix(0,ng,ng)
  diag(Isp) <- 1

#*************************************************
# initialize variables to save realizations
#*************************************************
uS <- matrix(0,ng,T)     #for mean, u
vS <- matrix(0,ng,T)     #for mean, v
u2S <- matrix(0,ng,T)    #for std dev, u
v2S <- matrix(0,ng,T)    #for std dev, v
uSreal <- list()         #save nreal u realizations at locsave
vSreal <- list()         #save nreal v realizations at locsave

#******************************************************
#   save theta and variance parameters, all iterations
#******************************************************
theta_uuS <- matrix(0,ngibbs,1)
theta_vvS <- matrix(0,ngibbs,1)
theta_uvS <- matrix(0,ngibbs,1)
theta_vuS <- matrix(0,ngibbs,1)
theta_upS <- matrix(0,ngibbs,1)
theta_vpS <- matrix(0,ngibbs,1)

s2uS <- matrix(0,ngibbs,1)
s2vS <- matrix(0,ngibbs,1)

#********************************************************************
#***************************************
#** main Gibbs sampler loop
#***************************************
#********************************************************************
for (k in 1:ngibbs){

  if (k %% 10 == 0)  print(k)

  #************************************
  #   sample u(t): t=1
  #************************************
  t <- 1
  Su <- SUdat[[t]]
  Eu <- EUdat[,t]
  Hs1 <- Hs[[t]]
  Hs2 <- HsHs[[t]]
  dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

  ct_v <-  theta_vv * v[,t] + theta_vp*DyP[,t]
  ct_u <- theta_uv * v[,t] + theta_up*DxP[,t]

  if (length(Su)== 0) {
#     At <- solve(HeHe/s2e +  (theta_vu^2/s2v)*Isp
#        + (theta_uu^2/s2u)*Isp + (1/s2u1)*Isp)
     At <- 1/(1/s2e +  (theta_vu^2/s2v) + (theta_uu^2/s2u) + (1/s2u1))
     bt <- t(t(Eu) %*% He/s2e +
           t(v[,t+1] - ct_v) * theta_vu/s2v +
           t(u[,t+1] - ct_u) * theta_uu/s2u +
           t(mu_u1)/s2u1)
     u[,t] <- At*bt + sqrt(At)*rnorm(ng)
    } else {
#     At <- solve(HeHe/s2e + Hs2/s2s + (theta_vu^2/s2v)*Isp
#               + (theta_uu^2/s2u)*Isp + (1/s2u1)*Isp)
     At <- 1/(1/s2e + dHs2/s2s + (theta_vu^2/s2v) + (theta_uu^2/s2u) + (1/s2u1))
     bt <- t(t(Eu) %*% He/s2e + t(Su) %*% Hs1/s2s +
           t(v[,t+1] - ct_v) * theta_vu/s2v +
           t(u[,t+1] - ct_u) * theta_uu/s2u +
           t(mu_u1)/s2u1)
     u[,t] <- At * bt + (At^.5)*rnorm(ng)
     }

  #************************************
  #   sample u(t): t=2,...,T-1
  #************************************
  for (t in 2:(T-1)) {
    Su <- SUdat[[t]]
    Eu <- EUdat[,t]
    Hs1 <- Hs[[t]]
    Hs2 <- HsHs[[t]]
    dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

    ct_v <- theta_vv * v[,t] + theta_vp*DyP[,t]
    ct_u <- theta_uv * v[,t] + theta_up*DxP[,t]
    ctm1_u <-  theta_uv * v[,t-1] + theta_up*DxP[,t-1]

    if (length(Su)== 0) {
#       At <- solve(HeHe/s2e + (theta_vu^2/s2v)*Isp
#                  + (theta_uu^2/s2u)*Isp + (1/s2u)*Isp)
       At <- 1/(1/s2e + (theta_vu^2/s2v) + (theta_uu^2/s2u) + (1/s2u))
       bt = t(t(Eu) %*% He/s2e +
               t(v[,t+1] - ct_v)*theta_vu/s2v +
               t(u[,t+1] - ct_u)*theta_uu/s2u +
               t(ctm1_u + theta_uu*u[,t-1])/s2u)
       u[,t] <- At*bt + sqrt(At)*rnorm(ng)
    }else {
#       At <- solve(HeHe/s2e + Hs2/s2s + (theta_vu^2/s2v)*Isp
#                  + (theta_uu^2/s2u)*Isp + (1/s2u)*Isp)
       At <- 1/(1/s2e + dHs2/s2s + (theta_vu^2/s2v)+ (theta_uu^2/s2u) + (1/s2u))
       bt = t(t(Eu)%*%He/s2e + t(Su)%*%Hs1/s2s +
               t(v[,t+1] - ct_v)*theta_vu/s2v +
               t(u[,t+1] - ct_u)*theta_uu/s2u +
               t(ctm1_u + theta_uu*u[,t-1])/s2u)
       u[,t] <- At * bt + (At^.5)*rnorm(ng)}

#   u[,t] <- At %*% bt + t(chol(At))%*%rnorm(ng)
  }

  #************************************
  #   sample u(t): t=T
  #************************************
  t <-T
  Su <- SUdat[[t]]
  Eu <- EUdat[,t]
  Hs1 <- Hs[[t]]
  Hs2 <- HsHs[[t]]
  dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

  ctm1_u <-  theta_uv * v[,t-1] + theta_up*DxP[,t-1]

  if (length(Su)== 0) {
#     At <- solve(HeHe/s2e + (1/s2u)*Isp)
     At <- 1/(1/s2e + 1/s2u)
     bt = t(t(Eu)%*%He/s2e +
              t(ctm1_u + theta_uu*u[,t-1])/s2u)
     u[,t] <- At*bt + sqrt(At)*rnorm(ng)
  }else {
#     At <- solve(HeHe/s2e + Hs2/s2s + (1/s2u)*Isp)
     At <- 1/(1/s2e + dHs2/s2s + 1/s2u)
     bt = t(t(Eu)%*%He/s2e + t(Su)%*%Hs1/s2s +
             t(ctm1_u + theta_uu*u[,t-1])/s2u)
     u[,t] <- At * bt + (At^.5)*rnorm(ng)
     }

#  u[,t] <- At %*% bt + t(chol(At))%*%rnorm(ng)


  #************************************
  #   sample v(t): t=1
  #************************************
  t <- 1
  Sv <- SVdat[[t]]
  Ev <- EVdat[,t]
  Hs1 <- Hs[[t]]
  Hs2 <- HsHs[[t]]
  dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

  ct_u <-  theta_uu * u[,t] + theta_up*DxP[,t]
  ct_v <- theta_vu * u[,t] + theta_vp*DyP[,t]

  if (length(Sv)== 0) {
#    At <- solve(HeHe/s2e +  (theta_uv^2/s2u)*Isp
#                + (theta_vv^2/s2v)*Isp + (1/s2v1)*Isp)
    At <- 1/(1/s2e +  (theta_uv^2/s2u) + (theta_vv^2/s2v) + 1/s2v1)
    bt <- t(t(Ev) %*% He/s2e +
              t(u[,t+1] - ct_u) * theta_uv/s2u +
              t(v[,t+1] - ct_v) * theta_vv/s2v +
              t(mu_v1)/s2v1)
     v[,t] <- At * bt + sqrt(At)*rnorm(ng)
  } else {
#    At <- solve(HeHe/s2e + Hs2/s2s + (theta_uv^2/s2u)*Isp
#                + (theta_vv^2/s2v)*Isp + (1/s2v1)*Isp)
    At <- 1/(1/s2e + dHs2/s2s + (theta_uv^2/s2u) + (theta_vv^2/s2v) + 1/s2v1)
    bt <- t(t(Ev) %*% He/s2e + t(Sv) %*% Hs1/s2s +
              t(u[,t+1] - ct_u) * theta_uv/s2u +
              t(v[,t+1] - ct_v) * theta_vv/s2v +
              t(mu_v1)/s2v1)
    v[,t] <- At * bt + (At^.5)*rnorm(ng)
    }

 # v[,t] <- At %*% bt + t(chol(At))%*%rnorm(ng)

  #************************************
  #   sample v(t): t=2,...,T-1
  #************************************
  for (t in 2:(T-1)) {
    Sv <- SVdat[[t]]
    Ev <- EVdat[,t]
    Hs1 <- Hs[[t]]
    Hs2 <- HsHs[[t]]
    dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

    ct_u <- theta_uu * u[,t] + theta_up*DxP[,t]
    ct_v <- theta_vu * u[,t] + theta_vp*DyP[,t]
    ctm1_v <-  theta_vu * u[,t-1] + theta_vp*DyP[,t-1]

    if (length(Sv)== 0) {
#      At <- solve(HeHe/s2e +  (theta_uv^2/s2u)*Isp
#                  + (theta_vv^2/s2v)*Isp + (1/s2v)*Isp)
      At <- 1/(1/s2e +  (theta_uv^2/s2u) + (theta_vv^2/s2v)+ 1/s2v)
      bt = t(t(Ev)%*%He/s2e +
               t(u[,t+1] - ct_u)*theta_uv/s2u +
               t(v[,t+1] - ct_v)*theta_vv/s2v +
               t(ctm1_v + theta_vv*v[,t-1])/s2v)
       v[,t] <- At * bt + sqrt(At)*rnorm(ng)
    }else {
#      At <- solve(HeHe/s2e + Hs2/s2s + (theta_uv^2/s2u)*Isp
#                  + (theta_vv^2/s2v)*Isp + (1/s2v)*Isp)
      At <- 1/(1/s2e + dHs2/s2s + (theta_uv^2/s2u) + (theta_vv^2/s2v) + 1/s2v)
      bt = t(t(Ev)%*%He/s2e + t(Sv)%*%Hs1/s2s +
               t(u[,t+1] - ct_u)*theta_uv/s2u +
               t(v[,t+1] - ct_v)*theta_vv/s2v +
               t(ctm1_v + theta_vv*v[,t-1])/s2v)
      v[,t] <- At * bt + (At^.5)*rnorm(ng)
      }

#    v[,t] <- At %*% bt + t(chol(At))%*%rnorm(ng)
  }

  #************************************
  #   sample v(t): t=T
  #************************************
  t <- T
  Sv <- SVdat[[t]]
  Ev <- EVdat[,t]
  Hs1 <- Hs[[t]]
  Hs2 <- HsHs[[t]]
  dHs2 <- data.matrix(diag(Hs2),rownames.force = NA)

  ctm1_v <-  theta_vu * u[,t-1] + theta_vp*DyP[,t-1]

  if (length(Sv)== 0) {
#    At <- solve(HeHe/s2e + (1/s2v)*Isp)
    At <- 1/(1/s2e + 1/s2v)
    bt = t(t(Ev)%*%He/s2e +
             t(ctm1_v + theta_vv*v[,t-1])/s2v)
    v[,t] <- At * bt + sqrt(At)*rnorm(ng)
  }else {
#    At <- solve(HeHe/s2e + Hs2/s2s  + (1/s2v)*Isp)
    At <- 1/(1/s2e + dHs2/s2s  + 1/s2v)
    bt = t(t(Ev)%*%He/s2e + t(Sv)%*%Hs1/s2s +
             t(ctm1_v + theta_vv*v[,t-1])/s2v)
    v[,t] <- At * bt + (At^.5)*rnorm(ng)}

 # v[,t] <- At %*% bt + t(chol(At))%*%rnorm(ng)

  rm(At,bt,Sv,Ev,Su,Eu,Hs1,Hs2,t)

  #************************************
  #   sample: thetas
  #************************************
  #**** theta_vv **********
  #tmpA <- 0
  #tmpb <- 0
  #for (t in 1:(T-1)){
  #   tmpA <- tmpA + t(v[,t]) %*% v[,t]/s2v
  #   Kvt  <- theta_vu*u[,t] + theta_vp*DyP[,t]
  #   tmpb <- tmpb + t(v[,t+1] - Kvt) %*% v[,t]/s2v
  #   }
  #A <- 1/(tmpA + 1/s2vv)
  #b <- tmpb + mu_vv/s2vv

  vtmp <- data.matrix(v,rownames.force = NA)
  utmp <- data.matrix(u,rownames.force = NA)

  tmpvv <- sum(diag(crossprod(vtmp[,1:T-1],vtmp[,1:T-1])))
  Kvt  <- theta_vu*utmp + theta_vp*DyP
  tdiff <- vtmp[,2:T] - Kvt[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,vtmp[,1:T-1])))/s2v
  A <- 1/((tmpvv/s2v) + 1/s2vv)
  b <- tmpb + mu_vv/s2vv
  theta_vv <- A*b + sqrt(A)*rnorm(1)

  #**** theta_uu **********
# tmpA <- 0
#  tmpb <- 0
#  for (t in 1:(T-1)){
#     tmpA <- tmpA + t(u[,t])%*%u[,t]/s2u
#     Kut <- theta_uv*v[,t] + theta_up*DxP[,t]
#     tmpb <- tmpb + t(u[,t+1] - Kut)%*%u[,t]/s2u
#  }
#  A <- 1/(tmpA + 1/s2uu)
#  b <- tmpb + mu_uu/s2uu

  tmpuu <- sum(diag(crossprod(utmp[,1:T-1],utmp[,1:T-1])))
  Kut  <- theta_uv*vtmp + theta_up*DxP
  tdiff <- utmp[,2:T] - Kut[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,utmp[,1:T-1])))/s2u
  A <- 1/((tmpuu/s2u) + 1/s2uu)
  b <- tmpb + mu_uu/s2uu
  theta_uu <- A*b + sqrt(A)*rnorm(1)

  #**** theta_vu **********
# tmpA <- 0
#  tmpb <- 0
#  for (t in 1:(T-1)){
#    tmpA <- tmpA + t(u[,t])%*%u[,t]/s2v
#    Kvt <- theta_vv*v[,t] + theta_vp*DyP[,t]
#    tmpb <- tmpb + t(v[,t+1] - Kvt)%*%u[,t]/s2v
#  }
#  A <- 1/(tmpA + 1/s2vu)
#  b <- tmpb + mu_vu/s2vu

  Kvt  <- theta_vv*vtmp + theta_vp*DyP
  tdiff <- vtmp[,2:T] - Kvt[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,utmp[,1:T-1])))/s2v
  A <- 1/((tmpuu/s2v) + 1/s2vu)
  b <- tmpb + mu_vu/s2vu
  theta_vu <- A*b + sqrt(A)*rnorm(1)

  #**** theta_uv **********
# tmpA <- 0
#  tmpb <- 0
#  for (t in 1:(T-1)){
#    tmpA <- tmpA + t(v[,t])%*%v[,t]/s2u
#    Kut <- theta_uu*u[,t] + theta_up*DxP[,t]
#    tmpb <- tmpb + t(u[,t+1] - Kut)%*%v[,t]/s2u
#  }
#  A <- 1/(tmpA + 1/s2uv)
#  b <- tmpb + mu_uv/s2uv

  Kut  <- theta_uu*utmp + theta_up*DxP
  tdiff <- utmp[,2:T] - Kut[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,vtmp[,1:T-1])))/s2u
  A <- 1/((tmpvv/s2u) + 1/s2uv)
  b <- tmpb + mu_uv/s2uv
  theta_uv <- A*b + sqrt(A)*rnorm(1)


  #**** theta_vp **********
# tmpA <- 0
#  tmpb <- 0
#  for (t in 1:(T-1)){
#    tmpA <- tmpA + t(DyP[,t])%*%DyP[,t]/s2v
#    Kvt <- theta_vv*v[,t] + theta_vu*u[,t]
#    tmpb <- tmpb + t(v[,t+1] - Kvt)%*%DyP[,t]/s2v
#  }
#  A <- 1/(tmpA + 1/s2vp)
#  b <- tmpb + mu_vp/s2vp

  tmpPP <- sum(diag(crossprod(DyP[,1:T-1],DyP[,1:T-1])))
  Kvt  <- theta_vv*vtmp + theta_vu*utmp
  tdiff <- vtmp[,2:T] - Kvt[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,DyP[,1:T-1])))/s2v
  A <- 1/((tmpPP/s2v) + 1/s2vp)
  b <- tmpb + mu_vp/s2vp
  theta_vp <- A*b + sqrt(A)*rnorm(1)

  #**** theta_up **********
#  tmpA <- 0
#  tmpb <- 0
#  for (t in 1:(T-1)){
#    tmpA <- tmpA + t(DxP[,t])%*%DxP[,t]/s2u
#    Kut <- theta_uu*u[,t] + theta_uv*v[,t]
#    tmpb <- tmpb + t(u[,t+1] - Kut)%*%DxP[,t]/s2u
#  }
#  A <- 1/(tmpA + 1/s2up)
#  b <- tmpb + mu_up/s2up
#   theta_up <- as.numeric(A*b + sqrt(A)*rnorm(1))

  tmpPP <- sum(diag(crossprod(DxP[,1:T-1],DxP[,1:T-1])))
  Kut  <- theta_uu*utmp + theta_uv*vtmp
  tdiff <- utmp[,2:T] - Kut[,1:T-1]
  tmpb <- sum(diag(crossprod(tdiff,DxP[,1:T-1])))/s2u
  A <- 1/((tmpPP/s2u) + 1/s2up)
  b <- tmpb + mu_up/s2up
  theta_up <- A*b + sqrt(A)*rnorm(1)

  rm(A,b,tmpb,tdiff,Kut,Kvt,tmpPP,tmpvv,tmpuu,utmp,vtmp)

  #************************************
  #   sample: s2u, s2v
  #************************************
  #********** s2u *********
  utmp <- data.matrix(u,rownames.force = NA)
  vtmp <- data.matrix(v,rownames.force = NA)
  Kut <- theta_uu*utmp+ theta_uv*vtmp + theta_up*DxP
  sumDiff=sum((utmp[,2:T] - Kut[,1:T-1])^2)
  qnew = qu + (T-1)*ng/2
  rnew = .5*(sumDiff)+ru
  s2u <- 1/rgamma(1,qnew,rnew)

  #********** s2v *********
  Kvt <- theta_vv*vtmp + theta_vu*utmp + theta_vp*DyP
  sumDiff=sum((vtmp[,2:T] - Kvt[,1:T-1])^2)
  qnew = qv + (T-1)*ng/2
  rnew = .5*(sumDiff)+rv
  s2v <- 1/rgamma(1,qnew,rnew)

  rm(utmp,vtmp,Kut,sumDiff,qnew,rnew)

  #**************************
  #  save variables
  #**************************
  theta_uuS[k] <- theta_uu
  theta_vvS[k] <- theta_vv
  theta_uvS[k] <- theta_uv
  theta_vuS[k] <- theta_vu
  theta_upS[k] <- theta_up
  theta_vpS[k] <- theta_vp
  s2uS[k] <- s2u
  s2vS[k] <- s2v

  if (k > nburn){
     uS <- uS + data.matrix(u, rownames.force = NA)
     vS <- vS + data.matrix(v, rownames.force = NA)
     u2S <- u2S + data.matrix(u, rownames.force = NA)^2
     v2S <- v2S + data.matrix(v, rownames.force = NA)^2
     if  ((k-nburn) %% real_int == 0){
         rcnt = rcnt+1;
         uSreal[[rcnt]] = data.matrix(u, rownames.force = NA)
         vSreal[[rcnt]] = data.matrix(v, rownames.force = NA)}
  }

}

uS <- uS/(ngibbs-nburn)
vS <- vS/(ngibbs-nburn)
u2S <- u2S/(ngibbs-nburn)
v2S <- v2S/(ngibbs-nburn)

#*** std dev
nn <- ngibbs-nburn
uSTD <- ((u2S*nn - ((uS*nn)^2)/nn)/(nn-1))^.5;
vSTD <- ((v2S*nn - ((vS*nn)^2)/nn)/(nn-1))^.5;


return(list(uS=uS,vS=vS,uSTD=uSTD,vSTD=vSTD,uSreal=uSreal,vSreal=vSreal,
            theta_uuS=theta_uuS,theta_vvS=theta_vvS,theta_uvS=theta_uvS,
            theta_vuS=theta_vuS,theta_upS=theta_upS,theta_vpS=theta_vpS,
            s2uS=s2uS,s2vS=s2vS))

}
