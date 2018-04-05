### calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (10:23) 
## Version: 
## Last-Updated: apr  5 2018 (15:54) 
##           By: Brice Ozenne
##     Update #: 191
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:


## * Documentation - calibrateType1
##' @title Simulation Study Assessing Bias and Type 1 Error
##' @description Perform a simulation study over one or several sample size
##' to assess the bias of the estimate
##' and the type 1 error of the Wald test and robust Wald test
##' @name calibrateType
##' 
##' @param object a \code{lvm} object defining the generative model.
##' @param null [character vector] names of the coefficient whose value will be tested against 0. 
##' @param n [integer vector, >0] sample size(s) considered in the simulation study.
##' @param n.rep [integer, >0] number of simulations per sample size.
##' @param coef.value [named numeric vector] values for each coefficient of the model.
##' Can also be \code{NULL}: in such a case the coefficients are set to default values decided by lava (usually 0 or 1).
##' @param dir.save [character] path to the directory were the results should be exported.
##' Can also be \code{NULL}: in such a case the results are not exported.
##' @param label.file [character] element to include in the file name.
##' @param n.true [integer, >0] sample size at which the estimated coefficients will be a reliable approximation of the true coefficients.
##' @param check.true [logical] should \code{coef.value} be compared to the estimate coefficients when the sample size equals \code{n.true}?
##' @param round.true [integer, >0] the number of decimal places to be used for the true value of the coefficients. No rounding is done if \code{NULL}.
##' @param bootstrap [logical] should bootstrap resampling be performed?
##' @param type.bootstrap [character vector]
##' @param n.bootstrap [integer, >0] the number of bootstrap sample to be used for each bootstrap.
##' @param seed [integer, >0] seed value that will be set at the beginning of the simulation to enable eproducibility of the results.
##' Can also be \code{NULL}: in such a case no seed is set.
##' @param trace [interger] should the execution of the function be trace. Can be 0, 1 or 2.
##' 
##' @return An object of class \code{calibrateType1}.
##' @seealso \code{link{autoplot.calibrateType1}} for a graphical display of the bias or of the type 1 error.
##' 
##' @author Brice Ozenne
##'
##' @examples
##' #### generative model ####
##' m.Sim <- lvm(c(Y1[mu1:sigma]~1*eta,
##'                Y2[mu2:sigma]~1*eta,
##'                Y3[mu3:sigma]~1*eta,
##'                eta~beta1*Group+beta2*Gender))
##' latent(m.Sim) <- ~eta
##' categorical(m.Sim, labels = c("M","F")) <- ~Gender
##'
##' vec.par <- c(mu2 = 0, ## mu1 is set to 0 by default
##'              mu3 = -0.37,
##'              "eta" = -2.74,
##'              beta1 = -1.20,
##'              beta2 = 0,
##'              sigma = 1.46,
##'              "eta~~eta" = 1.63)
##'
##' #### parameters to test ####
##' null <- c("mu2","beta2")
##'
##' #### launch simulation ####
##' \dontrun{
##' res <- calibrateType1(m.Sim, null = null, n = c(20,30,40), n.rep = 50, coef.value = vec.par)
##' autoplot(res, type = "bias")
##' autoplot(res, type = "type1error")
##' }
##' \dontshow{
##' res <- calibrateType1(m.Sim, null = null, n = c(20,30,40), n.rep = 2, coef.value = vec.par)
##' }

## * calibrateType1
##' @rdname calibrateType1
##' @export
calibrateType1 <- function(object, null, n, n.rep, coef.value,
                           dir.save = NULL, label.file = NULL,
                           n.true = 1e5, round.true = 2, check.true = TRUE,
                           bootstrap = FALSE, type.bootstrap = c("perc","stud","bca"), n.bootstrap = 1e3,
                           seed = NULL, trace = 2){

### ** prepare
    n.n <- length(n)

    ## coef names
    e0 <- estimate(object, data = lava::sim(object,n[n.n]))
    name.coef <- names(coef(e0))
    n.coef <- length(name.coef)
    
    ## coef type
    df.type <- coefType(e0, as.lava = FALSE)
    df.type <- df.type[df.type$name %in% name.coef,]
    type.coef <- setNames(df.type$detail, df.type$name)
    name.param <- df.type$param
    param2name <- setNames(df.type$name, df.type$param)

    ## coef true value
    coef.true <- setNames(rep(NA, n.coef), name.coef)
    if(!is.null(coef.value)){
        extraParam <- setdiff(sort(names(coef.value)), sort(name.param))
        missingParam <- setdiff(sort(name.param),sort(names(coef.value)))        
    }else{
        extraParam <- NULL
        missingParam <- name.param
    }
    missingCoef <- param2name[missingParam]

    if(!is.null(coef.value)){
        if(length(extraParam)>0){
            stop("Invalid argument \'coef.value\': some of the coefficient names do not match those of the estimated model \n",
                 "extra coefficients: \"",paste0(extraParam, collapse = "\" \""),"\"\n")
        }
        coef.true[param2name[names(coef.value)]] <- coef.value
    }
    if(length(missingParam)>0){
        message("Argument \'coef.value\' do not fully specify the value of the coefficients \n",
                "missing coefficients: \"",paste0(missingParam, collapse = "\" \""),"\"\n",
                "they will be set to their default value (zero or one) \n")
        if(trace>1){
            cat("* find missing coefficients ")
        }
        lavavalue <- coef(lava::estimate(object, data = lava::sim(object, n = n.true, p = coef.value)))
        if(!is.null(round.true)){
            lavavalue <- round(lavavalue, digits = round.true)
        }
        coef.true[missingCoef] <- lavavalue[missingCoef]
        if(trace>1){
            cat("- done \n")
        }
    }
    if(check.true){ 
        if(trace>1){
            cat("* check coefficient values match what can be estimated with n=",n.true," ",sep="")
        }
        
        if(length(missingParam)>0){
            coef.true.empirical <- lavavalue
        }else{
            coef.true.empirical <- coef(lava::estimate(object, data = lava::sim(object, n = n.true, p = coef.value)))
        }
        tol <- 1e-1
        ok.error <- setNames(rep(tol,n.coef), name.coef)
        ok.error[abs(coef.true)>tol*10] <- pmax(tol,coef.true[abs(coef.true)>tol*10]/10)
        if(any(abs(coef.true.empirical - coef.true) > ok.error)){
            warning("Some discrepancy between the estimated coefficient and the theoretical one for n=",n.true,"\n")
        }
        if(trace>1){
            cat("- done \n")
        }
    }


    ## null hypothesis
    n.null <- length(null)
    if(any(null %in% name.coef == FALSE)){
        incorrect.name <- which(null %in% name.coef == FALSE)
        if(all(null[incorrect.name] %in% names(param2name))){
            null[incorrect.name] <- as.character(param2name[null[incorrect.name]])
        }else{
            stop("Invalid argument \'null\': some of the coefficient names does not match those of the estimate model \n",
                 "incorrect names: \"",paste(null[incorrect.name], collapse = "\" \""),"\" \n")
        }
    }
    
    ## filename
    if(is.null(label.file)){label.file <- seed}
    filename_tempo.pvalue <- paste0("type1error-S",label.file,"(tempo).rds")
    filename_tempo.bias <- paste0("bias-S",label.file,"(tempo).rds")
    filename.pvalue <- gsub("\\(tempo\\)","",filename_tempo.pvalue)
    filename.bias <- gsub("\\(tempo\\)","",filename_tempo.bias)

    if(!is.null(dir.save)){
        lavaSearch2:::validPath(dir.save, type = "dir")
    }

### ** display
    if(trace>1){
        cat("* settings: \n")
        cat("  > simulation for n=",paste(n,collapse = " "),"\n",sep="")
        cat("  > model: \n")
        print(object)
        cat("  > coefficients: \n")
        print(coef.true)
        cat("  > bootstrap: ",bootstrap,"\n")
        if(!is.null(seed)){
            cat("  > seed: ",seed,"\n")
        }
        if(!is.null(dir.save)){
            cat("  > export results in ",dir.save,"\n")
        }
        
    }
    
### ** loop
    dt.pvalue <- NULL
    dt.bias <- NULL
    if(!is.null(seed)){
        set.seed(seed)
    }else{
        seed <- NA
    }

    if(trace>1){cat("* perform simulation: \n")}
    for(iN in 1:n.n){

        if(trace>0){cat("  > sample size=",n[iN],": ", sep = "")}
        n.tempo <- n[iN]

        for(iRep in 1:n.rep){
            if(trace>0){cat(iRep," ")}

            ## *** simulation
            dt.sim <- lava::sim(object, n = n.tempo, p = coef.value)

            ## *** model adjustement
            e.lvm <- lava::estimate(object, data = dt.sim)
            if(e.lvm$opt$convergence==1){next} ## exclude lvm that has not converged
            eS.lvm <- try(summary(e.lvm),silent = TRUE)
            if("try-error" %in% class(eS.lvm)){next}

            ## *** correction
            e.lvm.Satt <- e.lvm
            sCorrect(e.lvm.Satt) <- FALSE
            e.lvm.KR <- e.lvm
            suppressWarnings(sCorrect(e.lvm.KR, safeMode = TRUE) <- TRUE)
            ## check whether adjusted residuals could be computed (otherwise adjust.n=FALSE)
            test.warning <- inherits(attr(e.lvm.KR$sCorrect,"warning"),"try-error")
            coef.original <- coef(e.lvm)
            coef.corrected <- e.lvm.KR$sCorrect$param

            ## *** Wald test
            eS.Satt <- summary2(e.lvm.Satt)$coef
            eS.KR <- summary2(e.lvm.KR)$coef

            ## *** Robust Wald test        
            eS.robustSatt <- summary2(e.lvm.Satt, robust = TRUE)$coef
            eS.robustKR <- summary2(e.lvm.KR, robust = TRUE)$coef

            ## *** bootstrap
            if(bootstrap>0){
                e.boot <- eval(parse(text = "butils::bootReg(e.lvm, type = \"coef\", n.boot = n.bootstrap"))

                index.coef.boot <- match(null, name.coef)
                boot.perc <- summary(e.boot, p.value = TRUE, type = "perc", print = FALSE, index = index.coef.boot)
                boot.stud <- summary(e.boot, p.value = TRUE, type = "stud", print = FALSE, index = index.coef.boot)
                boot.bca <- summary(e.boot, p.value = TRUE, type = "bca", print = FALSE, index = index.coef.boot)
            }

            ## *** collect results (p.values)
            ls.iP <- list()

            ## z test
            ls.iP$p.Ztest <- eS.lvm$coef[null,"P-value"]
            ls.iP$p.robustZtest <-  2*(1-pnorm(abs(eS.robustSatt[null,"t-value"])))
        
            ## Satterwaite
            ls.iP$p.Satt <- eS.Satt[null,"P-value"]
            ls.iP$p.robustSatt <- eS.robustSatt[null,"P-value"]
        
            ## KR
            ls.iP$p.KR <- eS.KR[null,"P-value"]
            ls.iP$p.robustKR <- eS.robustKR[null,"P-value"]

            ## bootstrap
            if(bootstrap){
                ls.iP$p.bootPerc <- boot.perc[null,"p.value"]
                ls.iP$p.bootStud <- boot.stud[null,"p.value"]
                ls.iP$p.bootBca <- boot.bca[null,"p.value"]
            }
            ## metainformation
            iDT.pvalue <- cbind(data.frame(n = n.tempo,
                                           rep = iRep,
                                           seed = seed,
                                           nboot = n.bootstrap,
                                           niter = e.lvm.KR$sCorrect$opt$iterations,
                                           warning = test.warning,
                                           link = null,
                                           stringsAsFactors = FALSE),
                                do.call(cbind,ls.iP))
            rownames(iDT.pvalue) <- NULL
            dt.pvalue <- rbind(dt.pvalue, iDT.pvalue)

            ## *** collect result (bias)
            iDT.bias <- data.frame(n = n.tempo,
                                   rep = iRep,
                                   seed = seed,
                                   niter = e.lvm.KR$sCorrect$opt$iterations,
                                   warning = test.warning,
                                   estimate.truth = as.double(coef.true),
                                   estimate.ML = as.double(coef.original[name.coef]),
                                   estimate.MLcorrected = as.double(coef.corrected[name.coef]),
                                   name = names(coef.true),
                                   type = type.coef[name.coef],
                                   stringsAsFactors = FALSE)
            rownames(iDT.bias) <- NULL
            dt.bias <- rbind(dt.bias, iDT.bias)
        }

        ## *** export (tempo)
        if(!is.null(dir.save)){
            saveRDS(dt.pvalue, file = file.path(dir.save,filename_tempo.pvalue))
            saveRDS(dt.bias, file = file.path(dir.save,filename_tempo.bias))
        }
        if(trace>0){cat("\n")}
    }

    ## ** export
    if(!is.null(dir.save)){
        saveRDS(dt.pvalue, file = file.path(dir.save,filename.pvalue))
        saveRDS(dt.bias, file = file.path(dir.save,filename.bias))
    }

    out <- list(p.value = dt.pvalue,
                bias = dt.bias,
                coef.true = coef.true,
                null = null)
    class(out) <- append("calibrateType1",class(out))
    return(out)


}


######################################################################
### calibrateType1.R ends here
