dNR <- function(rt,pars){
  stats::dnorm(rt,mean=pars[,"m"],sd=pars[,"s"])
}

pNR <- function(rt,pars){
  stats::pnorm(rt,mean=pars[,"m"],sd=pars[,"s"])
}


rNR <- function(lR,pars,p_types=c("m","s"),ok=rep(TRUE,dim(pars)[1])){
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR))
  dt <- matrix(Inf,ncol=nrow(pars)/nr,nrow=nr)
  pars <- pars[ok,,drop=FALSE]
  dt[ok] <- stats::rnorm(dim(pars)[1],mean=pars[,"m"],sd=pars[,"s"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#' The Normal Race Model
#'
#' Model file to estimate the Normal Race Model (NR) in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `LNR()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**            |
#'  |-----------|-----------|---------------|-----------|----------------------------|---------------------------|
#'  | *m*       | -         | \[-Inf, Inf\]   | 1         |                            | Location parameter           |
#'  | *s*       | log       | \[0, Inf\]      | log(1)    |                            | Scale parameter           |
#'
#' Because the NR is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' In `design()`, `matchfun` can be used to automatically create a latent match
#' (`lM`) factor with levels `FALSE` (i.e., the stimulus does not match the accumulator)
#' and `TRUE` (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in the model formula, typically for parameters related to
#' the rate of accumulation (see the example below).
#'
#' Typically this model is used with only one accumulator to instantiate a
#' standard normal-normal setup
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_NRmE <- design(data = forstmann,model=NR,matchfun=matchfun,
#' formula=list(m~lM + E,s~1),
#' contrasts=list(m=list(lM=ADmat)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export
#'


NR <- function() {
  list(
    type="RACE",
    p_types=c("m" = 0,"s" = log(1)),
    Ntransform=function(x,use=NULL) {
      if (is.null(use)) {
        x[,dimnames(x)[[2]] != "m"] <- exp(x[,dimnames(x)[[2]] != "m"])
      } else if (!all(use=="m")) {
        ok <- use[use != "m"]
        x[,ok] <- exp(x[,ok])
      }
      x
    },
    # p_vector transform scaling parameter by s=1 assumed in lnr.R
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (!is.null(attr(dadm,"adaptive"))) pars <- do_adaptive(pars,dadm)
      attr(pars,"ok") <- (pars[,"s"] > 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars,ok) {
       if (is.null(lR)) ok else rNR(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dNR(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pNR(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}






