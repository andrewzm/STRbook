make_H <- function(grdxy,obsxy,int){
#
# Output:
# H - (nobs x ng) sparse matrix that maps grid locations to
#          ncep observation locations
# Input:
# grdxy - (ng x 2) lon,lat locations of prediction grid
# obsxy - (nobs x 2) lon,lat locations of nscat locations
# int   -  interval in x,y directions (in deg) around pred grid loc
#            to look for obs locs; search will look for obs within
#            int/2 of the grid location
#********************************************************

nobs <- dim(obsxy)[1]
ng <- dim(grdxy)[1]

H <- Matrix(0,nrow=nobs,ncol=ng,sparse=TRUE)
#H <- matrix(0,nobs,ng)

for (k in 1:nobs){
   olat <- obsxy[k,2]
   olon <- obsxy[k,1]

   fnd <- which((olat > grdxy[,2]-int/2) & (olat <= grdxy[,2]+int/2)
         & (olon > grdxy[,1] - int/2) & (olon <= grdxy[,1]+int/2))

   if (length(fnd) != 0){
      H[k,fnd]= 1}
}
return(H)
}


make_Dy_sparse <- function(nx,ny){
    #
    # function to take centered difference in the y direction
    # assumptions:
    #  (1) spatial grid on which differences are to be taken is
    #         has an extra column on each side and an extra row
    #         on each side
    #  (2) vector order of observations and grid starts in the upper
    #         left; e.g.,
    #       Prediction grid: nx=3, ny=2
    #             1 3 5
    #             2 4 6
    #         so, vector would be t(c(1,2,3,4,5,6))
    #       Observation grid: Nx = 5, Ny = 4
    #             1 5 9  13 17
    #             2 6 10 14 18
    #             3 7 11 15 19
    #             4 8 12 16 20
    #
    # Output: Dy - a (Nx*Ny) x (nx*ny) sparse matrix
    # Input:  nx - x-dim of prediction grid
    #         ny - y-dim of prediction grid
    #
    #*********************************************************
    #library('Matrix')

    Nx <- nx+2
    Ny <- ny+2
    N  <- Nx*Ny
    n  <- nx*ny

    Dy=Matrix(0,nrow=n,ncol=N,sparse=TRUE)
    #  Dy = matrix(0,n,N)

    strt <- ny
    cnt <- 0

    for (i in 1:nx){
        strt <- strt+2
        for (j in 1:ny){
            cnt <- cnt+1
            Dy[cnt,strt+(i-1)*ny+j] <- 1
            Dy[cnt,strt+(i-1)*ny+j+2] <- -1
        }
    }
    return(Dy)
}

make_Dx_sparse <- function(nx,ny){
    #
    # function to take centered difference in the x direction
    # assumptions:
    #  (1) spatial grid on which differences are to be taken is
    #         has an extra column on each side and an extra row
    #         on each side
    #  (2) vector order of observations and grid starts in the upper
    #         left; e.g.,
    #       Prediction grid: nx=3, ny=2
    #             1 3 5
    #             2 4 6
    #         so, vector would be t(c(1,2,3,4,5,6))
    #       Observation grid: Nx = 5, Ny = 4
    #             1 5 9  13 17
    #             2 6 10 14 18
    #             3 7 11 15 19
    #             4 8 12 16 20
    #
    # Output: Dx - a (Nx*Ny) x (nx*ny) sparse matrix
    # Input:  nx - x-dim of prediction grid
    #         ny - y-dim of prediction grid
    #
    #*********************************************************
    # library('Matrix')

    Nx <- nx+2
    Ny <- ny+2
    N  <- Nx*Ny
    n  <- nx*ny

    Dx=Matrix(0,nrow=n,ncol=N,sparse=TRUE)
    #Dx = matrix(0,n,N)

    strt <- -1
    cnt <- 0

    for (i in 1:nx){
        strt <- strt+2
        for (j in 1:ny){
            cnt <- cnt+1
            Dx[cnt,strt+(i-1)*ny+j] <- -1
            Dx[cnt,strt+(i-1)*ny+j+(2*Ny)] <- 1
        }
    }
    return(Dx)
}
