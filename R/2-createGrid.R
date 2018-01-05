### createGrid.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 31 2017 (16:40) 
## Version: 
## last-updated: okt  3 2017 (19:01) 
##           By: Brice Ozenne
##     Update #: 41
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:


# {{{ createGrid
#' @title Create a mesh for the integration
#' @description Create a mesh for the integration
#'
#' @param n the number of points for the mesh in the x direction.
#' @param xmin the minimal x value.
#' @param xmax the maximal x value.
#' @param d.y the number of dimensions for the triangle.
#' @param d.z the number of dimensions.
#' @param zmax the maximal z value (in absolute value).
#' @param fine should the mesh be displayed
#' @param double should the grid be just outside the region of interest? Otherwise it will be just inside.
#'
#' @details This create a mesh for integrating over a triangular surface using rectangles.
#' The domain is define by constrains on three types of variables:
#' \itemize{
#' \item the x variable: [unidimensional] that varies freely between [-xmax,-xmin] U [xmin,xmax].
#' \item the y variables: [dimension d.y] constrained to be lower in absolute value than the x variable.
#' In 2D this corresponds to 2 triangles, and in higher dimension to cones/hypercones.
#' \item the z variables: [dimension d.z] constrained to vary between [-zmax,zmax].
#' }
#' The intersection of these three conditions define the domain.
#'
#' The mesh is obtained slicing the triangles using rectangles.
#' 
#' @examples
#'
#' createGrid <- lavaSearch2:::createGrid
#' 
#' ## no z 
#' gridInt_2d <- createGrid(5, d.y = 1, xmin = 0, xmax = 4, 
#'                          d.z = 0, fine = FALSE, double = FALSE)
#' gridExt_2d <- createGrid(5, d.y = 1, xmin = 0, xmax = 4, 
#'                          d.z = 0, fine = FALSE, double = TRUE)
#' 
#' gridInt_4d <- createGrid(5, d.y = 3, xmin = 0, xmax = 4, 
#'                          d.z = 0, fine = FALSE, double = FALSE)
#' gridExt_4d <- createGrid(5, d.y = 3, xmin = 0, xmax = 4, 
#'                          d.z = 0, fine = FALSE, double = TRUE)
#' 
#' gridInt_2d <- createGrid(5, d.y = 1, xmin = 0, xmax = 4, 
#'                          d.z = 0, fine = TRUE, double = FALSE)
#' 
#' ## no z
#' gridIntZ1_2d <- createGrid(5, d.y = 1, xmin = 0, xmax = 4, 
#'                            d.z = 1, zmax = 2, fine = FALSE, double = FALSE)
#' gridExtZ1_2d <- createGrid(5, d.y = 1, xmin = 0, xmax = 4, 
#'                            d.z = 1, zmax = 2, fine = FALSE, double = TRUE)
#'  
#' gridIntZ2_4d <- createGrid(5, d.y = 3, xmin = 0, xmax = 4, 
#'                            d.z = 2, zmax = 2, fine = FALSE, double = FALSE)
#' gridExtZ2_4d <- createGrid(5, d.y = 3, xmin = 0, xmax = 4, 
#'                            d.z = 2, zmax = 2, fine = FALSE, double = TRUE)
#' 
createGrid <- function(n,
                       xmin, xmax, d.y, 
                       d.z, zmax,
                       fine, double){


    ## find step along the x axis
    if(xmin[1]==0){
        by <- xmax/n
    }else{
        seqTry <- seq(0,xmax, length.out = n)
        by <- xmin/sum(seqTry<xmin)
    }
    seqPointsX <- seq(xmin,xmax, by = by)
    n.seqX <- length(seqPointsX)-1

    ## name variables
    y.minNames <- paste0("y",1:d.y,".min")
    y.maxNames <- paste0("y",1:d.y,".max")
    all.Names <- c("x.min","x.max",y.minNames,y.maxNames,"weight","index")

    seqNames.min <- c("x.min",y.minNames)
    seqNames.max <- c("x.max",y.maxNames)
    
    # {{{ main grid
    grid.main <- NULL
    for(iX in 1:n.seqX){ #  iX <- 1
        grid.main <- rbind(grid.main,
                           c(seqPointsX[iX], seqPointsX[iX+1],
                             rep(-seqPointsX[iX+double],d.y), rep(seqPointsX[iX+double],d.y),
                             1, iX)
                           )
    }    
    grid.main <- as.data.table(grid.main)
    setnames(grid.main, old = names(grid.main), new = all.Names )    
    ## gg.main <- ggplot(grid.main, aes(xmin = x.min, ymin = y1.min, xmax = x.max, ymax = y1.max, fill = index))
    ## gg.main <- gg.main + geom_rect()
    ## gg.main <- gg.main + geom_abline(slope = 1,color = "red") + geom_abline(slope = -1,color = "red")
    ## gg.main
    # }}}
    
    # {{{ fine grid
    grid.fine <- NULL    
    if(fine){

        ls.mp <- lapply(1:d.y, function(x){c(-1,1)})
        grid.mp <- expand.grid(ls.mp)
        n.mp <- NROW(grid.mp)
        ls.index <- list(1:d.y,
                         (d.y+1):(2*d.y))

        for(iX in 1:n.seqX){ # iX <- 1
            for(iMP in 1:n.mp){ # iMP <- 1

                index.neg <- which(grid.mp[iMP,]<0)
                ls.index2 <- ls.index
                for(i.neg in index.neg){
                    ls.index2[[1]][i.neg] <- ls.index[[2]][i.neg]
                    ls.index2[[2]][i.neg] <- ls.index[[1]][i.neg]    
                }
                M.seq <- cbind(rep(seqPointsX[iX],d.y)*grid.mp[iMP,],
                               rep(seqPointsX[iX+1],d.y)*grid.mp[iMP,])
                
                grid.fine <- rbind(grid.fine,
                                   c(seqPointsX[iX], seqPointsX[iX+1],
                                     as.numeric(M.seq[unlist(ls.index2)]),
                                     1/2, iX+1)
                                   )
            }
        }
        grid.fine <- as.data.table(grid.fine)
        setnames(grid.fine, old = names(grid.fine), new = all.Names)
        grid.fine[,c(all.Names) := lapply(.SD,as.numeric)]
        ## gg.fine <- ggplot(grid.fine, aes(xmin = x.min, ymin = y1.min, xmax = x.max, ymax = y1.max, fill = index))
        ## gg.fine <- gg.fine + geom_rect()
        ## gg.fine <- gg.fine + geom_abline(slope = 1,color = "red") + geom_abline(slope = -1,color = "red")
        ## gg.fine
    }
    # }}}

    # {{{ Merge grids and remove empty cells
    grid.all <- rbind(grid.main, grid.fine)    
    grid.all <- grid.all[grid.all[, .I[.SD$y1.min!=.SD$y1.max]]] # remove empty rectangles
    # }}}
    
    # {{{ duplicate grid  (negative x  and positive x)
    grid.all <- grid.all[ ,"index" := .SD$index - min(.SD$index) + 1]

    grid.all2 <- copy(grid.all)    
    grid.all2[, "x.min" := -grid.all$x.max] 
    grid.all2[, "x.max" := -grid.all$x.min] # DO NOT use just x.min, we need grid.all$x.min due to the previous line
    grid.all2[, "index" := .SD$index + max(.SD$index)]

    grid <- rbind(grid.all,grid.all2)
    # }}}

    # {{{ add z coordinate
    if(d.z>0){
        z.minNames <- paste0("z",1:d.z,".min")
        z.maxNames <- paste0("z",1:d.z,".max")
        
        grid[,c(z.minNames) := -abs(zmax)]
        grid[,c(z.maxNames) := abs(zmax)]
        seqNames.min <- c(seqNames.min, z.minNames)
        seqNames.max <- c(seqNames.max, z.maxNames)
    }
    # }}}

    return(list(grid = grid,
                seqNames.min = seqNames.min,
                seqNames.max = seqNames.max))
}
# }}}

# {{{ createGrid2
#' @title Create a mesh for the integration
#' @description Create a mesh for the integration
#'
#' @param n the number of points for the mesh in one direction.
#' @param xmin the minimal x value.
#' @param xmax the maximal x value.
#' @param plot should the mesh be displayed
#' 
#' @examples
#' dt.res <- lavaSearch2:::createGrid2(20, xmin = 0, xmax = 4)
createGrid2 <- function(n, xmin, xmax, plot = FALSE){

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
	indexYheight <- grid[,.I[.SD$y == .SD$height]]
    grid[indexYheight & .SD$y_max>.SD$x_min, c("weight") := 1/2]
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
# }}}

# {{{ MCintGaus
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
# }}}


#----------------------------------------------------------------------
### createGrid.R ends here
