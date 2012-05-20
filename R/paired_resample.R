#' Paired resampling
#'
#' Run a paired resampling scheme to obtain one- or two-tailed p-values.
#'
#' @aliases paired_resample
#' @param level1 numeric vector for the first level
#' @param level2 numeric vector for the second level
#' @param ntails option for one or two tailed, default is 2
#' @param nresamp number of samples to draw, default is 10000
#' @export paired_resample
#' @author Rob Cribbie and Phil Chalmers 
#' @examples
#' \dontrun{
#' x<-rnorm(10,mean=-.2)
#' y<-rnorm(10)
#' paired_resample(x,y,ntails=1,nresamp=10)
#' }
paired_resample<-function (level1, level2, ntails=2,nresamp=10000) {
 differ<-level1-level2
 sumdiff<-sum(differ)
 n<-length(differ)
 newdiff<-abs(differ)
 results<-matrix(0,nrow=nresamp, ncol=1)
 for (i in 1:nresamp) {
   newdat<-rbinom(n,1,.5)
   newdat<-replace(newdat,newdat==0,-1)
   newdat2<-newdiff*newdat
   results[i,1]<-sum(newdat2)
 }
 if (ntails==1) {
  pval<-sum(results<sumdiff)/nresamp
  pval2<-sum(results>sumdiff)/nresamp
  pvala<-min(pval,pval2)
 }
 if (ntails==2) {
  pvala<-sum(abs(results)>abs(sumdiff))/nresamp
 }
out=list(ntails=ntails,nresamp=nresamp,sum_diffs=sumdiff,pval=pvala)
out
}


