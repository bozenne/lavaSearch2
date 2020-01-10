### biasCoxSnell.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: aug  2 2019 (10:20) 
## Version: 
## Last-Updated: jan 10 2020 (18:08) 
##           By: Brice Ozenne
##     Update #: 97
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .init_sscCoxSnell
.init_sscCoxSnell <- function(object,...){
    return(object,...)
}

## * .sscCoxSnell
.sscCoxSnell <- function(object){

    param <- object$sCorrect$ssc$param0
    name.param <- names(param)
    
    ## ** compute JJK    
    JJK <- .calcJJK(object)

    ## ** least squares
    Y <- (1/2) * sapply(name.param, function(iP){sum(JJK[,,iP] * object$sCorrect$vcov.param)})
    X <- object$sCorrect$information

    e.lm <- lm.fit(y = Y, x = X)
    newparam <- param - e.lm$coefficient
    
    ## ** export
    attr(newparam,"JJK") <- JJK
    attr(newparam,"lm") <- e.lm
    return(newparam)

}


## * .calcJJK
.calcJJK <- function(object){

    ## ** extract information
    dmu <- object$sCorrect$dmoment$dmu
    d2mu <- object$sCorrect$d2moment$d2mu
    dOmega <- object$sCorrect$dmoment$dOmega
    d2Omega <- object$sCorrect$d2moment$d2Omega

    missing.pattern <- object$sCorrect$missing$pattern
    name.pattern <- object$sCorrect$missing$name.pattern
    unique.pattern <- object$sCorrect$missing$unique.pattern
    n.pattern <- length(name.pattern)
    OmegaM1 <- object$sCorrect$moment$OmegaM1.missing.pattern
    
    name.param <- names(object$sCorrect$param)
    n.param <- length(name.param)
    n.cluster <- object$sCorrect$cluster$n.cluster
    
    grid.2meanD1.1varD1 <- object$sCorrect$skeleton$grid.2meanD1.1varD1
    grid.2meanD2.1meanD1 <- object$sCorrect$skeleton$grid.2meanD2.1meanD1
    grid.2varD2.1varD1 <- object$sCorrect$skeleton$grid.2varD2.1varD1
    n.grid.2meanD1.1varD1 <- NROW(grid.2meanD1.1varD1)
    n.grid.2meanD2.1meanD1 <- NROW(grid.2meanD2.1meanD1)
    n.grid.2varD2.1varD1 <- NROW(grid.2varD2.1varD1)

    ## ** prepare output    
    JJK <-  array(0, dim = c(n.param,n.param,n.param),
                  dimnames = list(name.param,name.param,name.param))

    ## ** loop over missing data pattern
    for(iP in 1:n.pattern){ ## iP <- 1
        iPattern <- name.pattern[iP]
        iIndex <- missing.pattern[[iPattern]]
        iY <- which(unique.pattern[iP,]==1)

        iOmegaM1 <- OmegaM1[[iPattern]]
        idmu <- .subsetList(dmu, indexRow = iIndex, indexCol = iY)
        idOmega <- .subsetList(dOmega, indexRow = iY, indexCol = iY)
        id2mu <- .subsetList2(d2mu, indexRow = iIndex, indexCol = iY)
        id2Omega <- .subsetList2(d2Omega, indexRow = iY, indexCol = iY)

        ## *** 1 second derivative and 1 first derivative regarding the variance
        if(n.grid.2varD2.1varD1>0){
            for(iGrid in 1:n.grid.2varD2.1varD1){ # iGrid <- 1
                iName1 <- grid.2varD2.1varD1[iGrid,"X"]
                iName2 <- grid.2varD2.1varD1[iGrid,"Y"]
                iName3 <- grid.2varD2.1varD1[iGrid,"Z"]

                ## term 1
                if(grid.2varD2.1varD1[iGrid,"d2XY"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2XY.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2XY.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName3]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                }

                ## term 2
                if(grid.2varD2.1varD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2XZ.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName2]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - 1/2 * sum(iDiag * n.cluster)
                }

                ## term 3
                if(grid.2varD2.1varD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2varD2.1varD1[iGrid,"d2YZ.Var2"]
                    iDiag <- diag(iOmegaM1 %*% id2Omega[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 %*% idOmega[[iName1]])
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + 1/2 * sum(iDiag * n.cluster)
                }
            }
        }
        
        ## *** 1 second derivative and 1 first derivative regarding the mean
        if(n.grid.2meanD2.1meanD1>0){
            for(iGrid in 1:n.grid.2meanD2.1meanD1){ # iGrid <- 1
                iName1 <- grid.2meanD2.1meanD1[iGrid,"X"]
                iName2 <- grid.2meanD2.1meanD1[iGrid,"Y"]
                iName3 <- grid.2meanD2.1meanD1[iGrid,"Z"]

                ## term 4
                if(grid.2meanD2.1meanD1[iGrid,"d2XY"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2XY.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2XY.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName3]])
                }

                ## term 5
                if(grid.2meanD2.1meanD1[iGrid,"d2XZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2XZ.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName2]])
                }

                ## term 6
                if(grid.2meanD2.1meanD1[iGrid,"d2YZ"]){
                    d2.Var1 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var1"]
                    d2.Var2 <- grid.2meanD2.1meanD1[iGrid,"d2YZ.Var2"]
                    JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + sum(id2mu[[d2.Var1]][[d2.Var2]] %*% iOmegaM1 * idmu[[iName1]])
                }
            }
        }

        ## *** 2 first derivative regarding the mean and one regarding the variance
        if(n.grid.2meanD1.1varD1>0){
            for(iGrid in 1:n.grid.2meanD1.1varD1){ # iGrid <- 1
                
                ## term 7
                iName1 <- grid.2meanD1.1varD1[iGrid,"Z"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Y"]
                value <- sum(idmu[[iName2]] %*% iOmegaM1 %*% idOmega[[iName1]] %*% iOmegaM1 * idmu[[iName3]])
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] + value


                ## term 8 
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Z"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"X"]
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value

                
                ## term 9
                iName1 <- grid.2meanD1.1varD1[iGrid,"X"]
                iName2 <- grid.2meanD1.1varD1[iGrid,"Y"]
                iName3 <- grid.2meanD1.1varD1[iGrid,"Z"]
                JJK[iName1,iName2,iName3] <- JJK[iName1,iName2,iName3] - value
            }
        }
        

    }

    return(JJK)
}

######################################################################
### biasCoxSnell.R ends here
