## * .onLoad
.onLoad <- function(lib, pkg="lavaSearch2") {

    # available methods to compute the distribution of the max statistic
    lava::lava.options(search.calcMaxDist = c("integration","boot-residual","boot-wild"),
                       search.p.adjust = c("fastmax", "max", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none","gof"),
                       search.calc.quantile.int = FALSE,
                       search.n.perm = 1e6,
                       method.estimate2 = "ols",
                       factor.dRvcov = 1/2
                       )
}

## * .onAttach
.onAttach <- function(lib, pkg="lavaSearch2") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

lava_categorical2dummy <- get("categorical2dummy", envir = asNamespace("lava"), inherits = FALSE)
lava_estimate.lvm <- get("estimate.lvm", envir = asNamespace("lava"), inherits = FALSE)
lava_procdata.lvm <- get("procdata.lvm", envir = asNamespace("lava"), inherits = FALSE)
