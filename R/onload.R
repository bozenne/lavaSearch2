.onLoad <- function(lib, pkg="lavaSearch2") {

    # available methods to compute the distribution of the max statistic
    lava::lava.options(calcDistMax.method = c("integration","bootstrap","boot-norm","boot-wild"))
}


.onAttach <- function(lib, pkg="lavaSearch2") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

lava_categorical2dummy <- get("categorical2dummy", envir = asNamespace("lava"), inherits = FALSE)
lava_estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)
