### InTriGauss.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 14 2017 (11:49) 
## Version: 
## last-updated: aug 28 2017 (09:50) 
##           By: Brice Ozenne
##     Update #: 93
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


##' @title Integration of a density over a triangle
##'
##' @description Integration of a density over a triangle
##' 
##' @param mu the mean vector
##' @param Sigma the variance-covariance matrix
##' @param df degree of freedom for a student distribution.
##' @param n number of points in one direction for the mesh
##' @param xmin the minimum value along the x axis
##' @param xmax the maximum value along the x axis
##' @param distribution Can be \code{"pmvnorm"} (normal distribution) or \code{"pvmt"} (student's t distribution)
##'
##' @examples
##' library(mvtnorm)
##' 
##' p <- 2
##' Sigma <- diag(p)
##' mu <- rep(0, p)
##' 
##' ## bivariate normal distribution
##' z2 <- qmvt(0.975, mean = mu, sigma = Sigma, df = 1e3)$quantile
##'
##' # compute integral
##' IntDensTri(mu = mu, Sigma = Sigma, n=5, xmin=0, xmax=3)-1/2
##' IntDensTri(mu = mu, Sigma = Sigma, n=5, xmin=0, xmax=5)-1/2
##' IntDensTri(mu = mu, Sigma = Sigma, n=5, xmin=0, xmax=100)-1/2
##'
##' IntDensTri(mu = mu, Sigma = Sigma, df = 5, n=5, xmin=0, xmax=3, distribution = "pmvt")-1/2
##' IntDensTri(mu = mu, Sigma = Sigma, df = 5, n=5, xmin=0, xmax=5, distribution = "pmvt")-1/2
##' IntDensTri(mu = mu, Sigma = Sigma, df = 5, n=5, xmin=0, xmax=100, distribution = "pmvt")-1/2
##'
##' ## trivariate normal distribution
##' p <- 3
##' Sigma <- diag(p)
##' mu <- rep(0, p)
##'
##' IntDensTri(mu = mu, Sigma = Sigma, n=5, xmin=c(0,0), xmax=c(10,10))-1/2
##' 
##' @author Brice Ozenne
##' @export 
IntDensTri <- function(mu, Sigma, df, n, xmin, xmax, 
                       distribution = "pmvnorm"){

    p <- length(mu)
    
    ## create the grid of points to integrate over
    grid <- createGrid(n, xmin[1], xmax[1], plot = FALSE)
    ## integrate the density over the grid
    ls.args <- list(mu,Sigma)    
    if(distribution=="pmvt"){
        names(ls.args) <- c("delta","sigma")
        ls.args$df <- df
    }else{
        names(ls.args) <- c("mean","sigma")
    }

    total.area <- 0
    if(p==2){
        grid$area <- apply(grid, 1, function(x){
            do.call(distribution, args = c(lower = list(c(x["x_min"],x["y_min"])),
                                           upper = list(c(x["x_max"],x["y_max"])),
                                           ls.args))
        })
        total.area <- grid[,sum(.SD$area*.SD$weight)]
    }else{

        gridZ <- data.frame(z_min = c(xmin[-1],-xmax[-1]),
                            z_max = c(xmax[-1],-xmin[-1]))
        nZ.grid <- NROW(gridZ)

        for(iZ in 1:nZ.grid){ # iZ <- 1
            grid$area <- apply(grid, 1, function(x){
                do.call(distribution, args = c(lower = list(c(x["x_min"],x["y_min"],gridZ[iZ,"z_min"])),
                                               upper = list(c(x["x_max"],x["y_max"],gridZ[iZ,"z_max"])),
                                               ls.args))
            })
            total.area <- total.area + grid[,sum(.SD$area*.SD$weight)]
        }

    }
        
    return(total.area)
}

#' @title Create a mesh for the integration
#' @description Create a mesh for the integration
#'
#' @param n the number of points for the mesh in one direction.
#' @param xmin the minimal x value.
#' @param xmax the maximal x value.
#' @param plot should the mesh be displayed
#' 
#' @examples
#' dt.res <- lavaSearch2:::createGrid(20, xmin = 0, xmax = 4)
createGrid <- function(n, xmin, xmax, plot = FALSE){

    if(xmin==0){
        by <- xmax/n
    }else{
        seqTry <- seq(0,xmax, length.out = n)
        by <- xmin/sum(seqTry<xmin)
    }
    
    seqPointsX <- seq(xmin,xmax, by = by)
    n.seqX <- length(seqPointsX)
    seqPointsY <- seq(0,xmax, by = by)
    delay <- max(0,which.min(abs(seqPointsX[1]-seqPointsY))-1)
    grid <- NULL
  
    for(iX in 1:n.seqX){ # iX <- 1
        grid <- rbind(grid,
                      cbind(x = seqPointsX[iX],
                            y = seqPointsY[1:(delay + iX)])
                      )
    }

    grid <- unique(as.data.table(grid))
    grid[, c("height") := max(.SD$y), by = "x"]
    grid[, c("x_min") := pmax(xmin,.SD$x-by/2)]
    grid[, c("y_min") := pmax(0,.SD$y-by/2)]
    grid[, c("x_max") := .SD$x+by/2]
    grid[, c("y_max") := .SD$y+by/2]
    grid[, c("weight") := 1]
    grid[.SD$y == .SD$height & .SD$y_max>.SD$x_min, c("weight") := 1/2]
    grid[, c("index") := 1:.N]
    #grid[x==0&y==0,c("y_min","ymax",weight) := c(0,0,0)]

    if(plot){
        plot(grid[,list(.SD$x,.SD$y)], pch = 20+2*grid$weight, xlim = c(0,max(grid$x_max)), ylim = c(0,max(grid$y_max)))
        graphics::points(grid[,list(.SD$x_min,.SD$y_min)], col = grDevices::rainbow(max(grid$index)), pch = 20)
        graphics::points(grid[,list(.SD$x_min,.SD$y_max)], col = grDevices::rainbow(max(grid$index)), pch = 20)
        graphics::points(grid[,list(.SD$x_max,.SD$y_min)], col = grDevices::rainbow(max(grid$index)), pch = 20)
        graphics::points(grid[,list(.SD$x_max,.SD$y_max)], col = grDevices::rainbow(max(grid$index)), pch = 20)
    }

    # duplicate the grid for negative y
    grid.YnegXpos <- copy(grid)
    grid.YnegXpos[, c("y","y_min","y_max"):= list(-.SD$y,-.SD$y_max,-.SD$y_min)]
    # duplicate the grid for negative x
    grid.pos <- rbind(grid,grid.YnegXpos)
    grid.neg <- copy(grid.pos)
    grid.neg[, c("x","x_min","x_max"):= list(-.SD$x,-.SD$x_max,-.SD$x_min)]

    fullGrid <- rbind(grid.pos,grid.neg)

    return(fullGrid)
}

#' 
#' MCintGaus(c(0,0),diag(1,2),n=1e4)
#' MCintGaus(c(0,0),diag(1,2),n=1e4, xmax = 10)
#'
#' 
## MCintGaus <- function(mu, Sigma, n, xmin, xmax = 5,
##                       distribution = "dmvnorm"){

##   p <- length(mu)
##   X <- GenerateSpherePoints(nrPoints = n, nrDim = p, r = xmax)
##   X <- X[apply(abs(X),1,which.max)==1,]
##   out <- sum(apply(X,1,dmvnorm,mean = mu, sigma = Sigma)*pi*xmax^2)/n
##   return(out)
  
## }

## # https://stackoverflow.com/questions/5016806/generating-multidimensional-data
## GenerateSpherePoints <- function(nrPoints,nrDim,center=rep(0,nrDim),r=1){
##   #generate the polar coordinates!
##   x <-  matrix(runif(nrPoints*nrDim,-pi,pi),ncol=nrDim)
##   x[,nrDim] <- x[,nrDim]/2
##   #recalculate them to cartesians
##   sin.x <- sin(x)
##   cos.x <- cos(x)
##   cos.x[,nrDim] <- 1  # see the formula for n.spheres
  
##   y <- sapply(1:nrDim, function(i){
##     if(i==1){
##       cos.x[,1]
##     } else {
##       cos.x[,i]*apply(sin.x[,1:(i-1),drop=F],1,prod)
##     }
##   })*sqrt(runif(nrPoints,0,r^2))
  
##   y <-  as.data.frame(
##     t(apply(y,1,'+',center))
##   )
  
##   names(y) <- make.names(seq_len(nrDim))
##   y
## }

#----------------------------------------------------------------------
### InTriGauss.R ends here
