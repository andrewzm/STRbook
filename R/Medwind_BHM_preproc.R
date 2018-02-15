#' @export
Medwind_BHM_preproc <- function(Edat, Sdat, Predlocs, Inparm){

    EUdat <- Edat$EUdat
    EVdat <- Edat$EVdat
    EPdat <- Edat$EPdat
    ECMWFxylocs <- Edat$ECMWFxylocs

    Sxylocs <- Sdat$Sxylocs
    SUdat <- Sdat$SUdat
    SVdat <- Sdat$SVdat

    Mgridxylocs <- Predlocs

    gspdeg <- Inparm$gspdeg      # prediction grid spacing (deg)
    srad <- Inparm$srad          # search for ECMW data locs within 1/2 srad (deg) of pred grid loc
    erad <- Inparm$erad          # search for scat data locs within 1/2 erad (deg) of pred grid loc
    hx <- Inparm$hx              #x-grid (pred and ECMWF) spacing in meters
    hy <- Inparm$hy              #y-grid (pred and ECMWF) spacing in meters
    s2s <- Inparm$s2s            #scatterometer measurement error variance
    s2e <- Inparm$s2e            #ecmwf analysis measurement error variance
    mu_pri <- Inparm$mu_pri      #prior mean for theta parameters and initial wind components
    s2_pri <- Inparm$s2_pri      #prior variance for theta parms and intial wind components
    IGshape <- Inparm$IGshape    #inverse gamma prior shape parameter for process variances
    IGrate <- Inparm$IGrate      #inverse gamma prior rate parameter for process variances

    #***********************************
    #** Parameters for Data and Grid
    #***********************************
    T <- dim(EUdat)[2]            #number of time periods
    ng <- dim(Mgridxylocs)[1]     #number of pred grid locations
    ne <- dim(EUdat)[1]           #number of ECMWF grid locations
    mnLon <- min(Mgridxylocs[,1]) #min longitude pred grid
    mxLon <- max(Mgridxylocs[,1]) #max longitude pred grid
    mnLat <- min(Mgridxylocs[,2]) #min latitude pred grid
    mxLat <- max(Mgridxylocs[,2]) #max latitude pred grid
    nx <- length(seq(mnLon,mxLon,gspdeg))   #number of x (lon) pred grid locations
    ny <- length(seq(mnLat,mxLat,gspdeg))   #number of y (lat) pred grid locations

    ecmwf_search_r <- erad         #deg to search for ECMWF data centered on grid loc
    scat_search_r <- srad          #deg to search for scat data centered on grid loc

    #***********************************
    #** Preprocess Sparse Matrices
    #***********************************
    #  identity matrix
    #Isp <- Matrix(0,ng,ng)
    #diag(Isp) <- 1

    #  make incidence matrix for ECMWF analysis data
    #
    He <- make_H(Mgridxylocs,ECMWFxylocs,ecmwf_search_r);
    HeHe <- t(He)%*%He   #Kaw

    # make incidence matrices for scatterometer data
    #
    Hs = list();
    HsHs = list();
    for (t in 1:T){
        if (dim(Sxylocs[[t]])[1] == 0){
            Hs[[t]] <- Matrix(0,0,0)
            HsHs[[t]] <- Matrix(0,0,0)}
        else{
            Hs[[t]] <- make_H(Mgridxylocs,Sxylocs[[t]],scat_search_r)
            HsHs[[t]] <- t(Hs[[t]])%*%Hs[[t]]}
    }

    #******************************************
    #  make difference operators
    #   NOTE: assumes following vector order from grid; e.g.,
    #      for nx=3 and ny=2, this requires the following order:
    #            1  4  7  10
    #            2  5  8  11
    #            3  6  9  12
    #*****************************************

    Dx <- make_Dx_sparse(nx,ny)
    Dy <- make_Dy_sparse(nx,ny)
    Dx <- (1/(2*hx))*Dx;            #scale by grid spacing
    Dy <- (1/(2*hy))*Dy;            #scale by grid spacing

    #*****************************************************
    #  apply difference operators to ECMWF pressure
    #*****************************************************
    DxP <- Dx%*%data.matrix(EPdat, rownames.force = NA)*100
    DyP <- Dy%*%data.matrix(EPdat, rownames.force = NA)*100

    #***************************************************
    # prior hyperparameters
    #***************************************************
    # s2s <- 1   #scatterometer measurement error variance
    #s2e <- 10  #ecmwf analysis measurement error variance

    #********************************************
    # parameter prior means and variances
    #********************************************
    mu_uu <- mu_pri
    mu_vv <- mu_pri
    mu_uv <- mu_pri
    mu_vu <- mu_pri
    mu_up <- mu_pri
    mu_vp <- mu_pri
    s2uu <- s2_pri
    s2vv <- s2_pri
    s2uv <- s2_pri
    s2vu <- s2_pri
    s2up <- s2_pri
    s2vp <- s2_pri

    qu <- IGshape
    ru <- IGrate
    qv <- IGshape
    rv <- IGrate

    mu_u1 <- matrix(0,ng,1)
    mu_v1 <- matrix(0,ng,1)
    s2u1 <- s2_pri
    s2v1 <- s2_pri

    #**********************************************
    # starting values for u and v
    #**********************************************
    fnd <- which((ECMWFxylocs[,1] >= mnLon) & (ECMWFxylocs[,1] <= mxLon)
                 & (ECMWFxylocs[,2] >= mnLat) & (ECMWFxylocs[,2] <= mxLat))
    u <- EUdat[fnd,]
    v <- EVdat[fnd,]

    #**alternative way to do it when grids don't match up
    #Htrans <- solve(HeHe + .01*diag(ng))%*%t(He)
    #u <- Htrans%*%data.matrix(EUdat, rownames.force = NA)
    #v <- Htrans%*%data.matrix(EVdat, rownames.force = NA)

    # other starting values
    theta_uu = .9
    theta_vv = .9
    theta_uv = 0
    theta_vu = 0
    theta_up = 0
    theta_vp = 0
    s2u = 4
    s2v = 4

    #**********************************************
    #  make lists to return
    #*********************************************
    #  data list
    Mdata <- list(EVdat = EVdat,
                  EUdat = EUdat,
                  DxP=DxP,
                  DyP=DyP,
                  ECMWFxylocs=ECMWFxylocs,
                  SUdat=SUdat,
                  SVdat=SVdat,
                  Sxylocs=Sxylocs,
                  necmwf = ne,
                  T = T)

    #  prediction grid and incidence matrices list
    Mgrid <- list(Mgridxylocs = Mgridxylocs,
                  He = He,
                  HeHe = HeHe,
                  Hs = Hs,
                  HsHs=HsHs,
                  ng=ng,
                  nx=nx,
                  ny=ny)

    #  prior hyperparameters list
    Mpriors <- list(
        mu_uu = mu_uu,
        mu_vv = mu_vv,
        mu_uv = mu_uv,
        mu_vu = mu_vu,
        mu_up = mu_up,
        mu_vp = mu_vp,
        s2uu  = s2uu,
        s2vv  = s2vv,
        s2uv  = s2uv,
        s2vu  = s2vu,
        s2up  = s2up,
        s2vp  = s2vp,
        qu  = qu,
        ru  = ru,
        qv  = qv,
        rv  = rv,
        mu_u1 = mu_u1,
        mu_v1 = mu_v1,
        s2u1 = s2u1,
        s2v1 = s2v1,
        s2e = s2e,
        s2s = s2s)

    #   starting values list
    Mstrt <- list(
        u = u,
        v = v,
        theta_uu = theta_uu,
        theta_vv = theta_vv,
        theta_uv = theta_uv,
        theta_vu = theta_vu,
        theta_up = theta_up,
        theta_vp = theta_vp,
        s2u = s2u,
        s2v = s2v)

    return(list(Mdata=Mdata,Mgrid=Mgrid,Mpriors=Mpriors,Mstrt=Mstrt))
}
