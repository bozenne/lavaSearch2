### residplot.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 29 2017 (11:52) 
## Version: 
## last-updated: jan 15 2018 (11:46) 
##           By: Brice Ozenne
##     Update #: 119
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * Documentation - residplot
#' @title Plot the Residuals of a LVM Object
#' @description Plot the residuals for each outcome of a lvm object.
#' 
#' @name residplot
#'
#' @param object a lvm model.
#' @param res.variables the endogenous variable for which the residuals should be displayed.
#' @param obs.variables same as res.variables or a variable present in the model
#' @param sd.kernel the standard deviation of the kernel used to smooth the variance.
#' @param kernel the type of kernel used to smooth the variance.
#' @param plot.weights should the weights used to compute the variance be displayed?
#' @param ncol the number of columns in the graphical display.
#' @param smooth.mean should the mean be displayed across the obs.variables values?
#' @param smooth.sd should the variance be displayed across the obs.variables values?
#' @param plot should the fit be displayed in a graphical window.
#' @param ... additional arguments.
#' 
#' @return a list containing:
#' \itemize{
#' \item plot the display of the fit
#' \item data the dataset used to make the display
#' }
#' 
#' @examples 
#' m <- lvm(y1~eta+x1,y2~eta,y3~eta+x3)
#' latent(m) <- ~ eta
#' dd <- lava::sim(m,100) ## Simulate 100 observations from model
#' e <- estimate(m, dd) ## Estimate parameters
#'
#' res <- residplot(e)
#' residplot(e, obs.variables = "x1")
#'
#' m <- lvm(y~x)
#' distribution(m,~y) <- function(n,mean,x) rnorm(n,mean,exp(x)^.5)
#' d <- lava::sim(m,1e3)
#' 
#' e <- estimate(m,data = d)
#'
#' residplot(e)
#' # residplot(e, plot.weights = TRUE)
#' residplot(e, res.variables = "y", obs.variables = "x")
#' @export
residplot <- function (object, ...) {
  UseMethod("residplot", object)
}

## * residplot.lvmfit
#' @rdname residplot
#' @export
residplot.lvmfit <- function(object, res.variables = endogenous(object), obs.variables = res.variables,
                             sd.kernel = 0.5, kernel = "dnorm", plot.weights = FALSE, ncol = NULL,
                             smooth.mean = TRUE, smooth.sd = TRUE, plot = TRUE,...){

    fitted <- observed <- NULL ## [:for CRAN check] data.table
    
    if(any(obs.variables %in% lava::vars(object) == FALSE)){
        missing.var <- obs.variables[obs.variables %in% lava::vars(object) == FALSE]
        stop("variable \"",paste(missing.var, collapse = "\"\""),"\" not found \n",
             "all obs.variables must be in the model \n")
    }
    if(any(res.variables %in% lava::endogenous(object) == FALSE)){
        missing.var <- res.variables[res.variables %in% lava::endogenous(object) == FALSE]
        stop("outcome \"",paste(missing.var, collapse = "\"\""),"\" not found \n",
             "all res.variables must be endogenous \n")
    }
    if(length(obs.variables) !=  length(res.variables)){
        if(length(obs.variables)==1){
            obs.variables <- rep(obs.variables,length(res.variables))
        }else{
            stop("obs.variables and res.variables must have the same length \n")
        }
    }
    if(!identical(res.variables,obs.variables)){
        x.lab <- unique(obs.variables)
        if(length(x.lab)>1){
            stop("obs.variables must either be identical to res.variables \n",
                 "or contain only the name of one variable from the model \n")
        }
    }else{
        x.lab <- "observed outcome"
    }

    ## observed values
    res <- stats::model.frame(object)
    if(is.data.table(res)){
        data <- res[,.SD,.SDcols=obs.variables]
    }else{
        data <- as.data.table(res[,obs.variables])
    }
    names(data) <- res.variables
    data[, "XXXXIdXXXX" := 1:.N]

    dtL.data <- melt(data, id.vars = "XXXXIdXXXX",
                     variable.name = "endogenous")
    
    ## predicted values
    #predicted <- as.data.table(stats::predict(object)[,variables])
    #predicted[, "XXXXIdXXXX" := 1:.N]

    ## residuals values
    dt.residuals <- as.data.table(stats::residuals(object)[,res.variables])
    ## dt.residuals <- as.data.table(stats::stats::predict(object)[,res.variables])
    names(dt.residuals) <- res.variables
    dt.residuals[, "XXXXIdXXXX" := 1:.N]

    dtL.residuals <- melt(dt.residuals, id.vars = "XXXXIdXXXX",
                          variable.name = "endogenous")
    
    dtL.all <- merge(x = dtL.data, y = dtL.residuals, by = c("XXXXIdXXXX","endogenous"))
    setnames(dtL.all, old = c("value.x","value.y"), new = c("observed","fitted"))
    
    ## compute standard deviation
    setkeyv(dtL.all, "observed")
    dtL.all[, "sdY" := smoothSD(Y = fitted, time = observed,
                                plot.weights = plot.weights,
                                sd.kernel = sd.kernel, kernel = kernel),
            by = "endogenous"]
    
    ## display
    gg <- ggplot(dtL.all, aes_string(x = "observed", y = "fitted"))
    gg <- gg + geom_point()
    if(smooth.mean){
        gg <- gg + geom_smooth(method = "lm",aes(color = "mean"))
    }
    if(smooth.sd){
        sdY <- NULL ## [:for CRAN check] we don't want to move to aes_string because color is indeed a character
        gg <- gg + geom_line(aes(y = sdY, color = "standard devation"))
    }
    if(is.null(ncol)){
        ncol <- round(sqrt(length(res.variables)))
    }
    gg <- gg + facet_wrap(~endogenous, ncol = ncol, labeller = label_both)
    gg <- gg + xlab(x.lab) + ylab("residuals")
    if(plot){
        print(gg)
    }
    
    ## export
    return(invisible(list(data = dtL.all,
                          plot = gg)))

    
}

## * smoothSD
smoothSD <- function(Y, time, sd.kernel, kernel, plot.weights){

    n <- length(Y)    
    dt <- data.table(Y=Y,time=scale(time, center = stats::median(time), scale = stats::mad(time)), index = 1:n)
    dt[,time:=time-min(time)]
    setkey(dt,time)
    
    res <- sapply(1:n, function(i){ # i <- 25
        weight <- abs(do.call(kernel, args = list(x = dt$time[i]-dt$time, sd = sd.kernel)))
        #weight[weight<1e-6] <- 0
        weight <- weight/sum(weight)
        if(plot.weights){graphics::plot(weight)}
        w.mean <- stats::weighted.mean(Y, weight)
        w.sd <- sqrt(sum(weight * (Y - w.mean)^2))
        return(w.sd)
    })
    return(res)
}

#----------------------------------------------------------------------
### residplot.R ends here
