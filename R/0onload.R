## * .onLoad
.onLoad <- function(lib, pkg="lavaSearch2") {

    lava::lava.options(search.calc.quantile.int = FALSE, ## hidden argument for modelsearch2
                       search.type.information = "E", ## hidden argument for modelsearch2
                       ## search.perm.stat = "exact", ## hidden argument for modelsearch2 (otherwise "exact")
                       method.estimate2 = "ols", ## hidden argument for sCorrect
                       factor.dRvcov = 1/2,
                       ssc = "residuals",
                       df = "satterthwaite",
                       df.robust = 1
                       )

    lava::addhook("gaussian_weight.estimate.hook", hook = "estimate.hooks")

}

## * .onAttach
.onAttach <- function(lib, pkg="lavaSearch2") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

lava_categorical2dummy <- get("categorical2dummy", envir = asNamespace("lava"), inherits = FALSE)
lava_estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)
lava_procdata.lvm <- get("procdata.lvm", envir = asNamespace("lava"), inherits = FALSE)
