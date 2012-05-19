#' Linear comparisons
#'
#' Computes an F test for a linear comparison among group means. Uses procedures described in Maxwell &
#'  Delaney (2004; Designing Experiments and Analyzing Data), pages 1157-1162. This routine is suitable for
#'   one-way designs.
#'
#' @section If variances are assumed to differ across groups, the routine calculates MS-within and adjusts
#'  df-within using procedures developed by Welch (1938) & Shatterwaite (1946), as described in Maxwell &
#'   Delaney (2004; Designing Experiments and Analyzing Data), pages 165-168
#' @aliases linear.comparison
#' @param y Dependent variable
#' @param group Grouping variable
#' @param c.weights A list containing weights for multiple linear contrasts; or a numeric vector containing
#'   weights for 1 contrast
#' @param alpha compute 100*(1-alpha)\% confidence interval; default is 0.05
#' @param var.equal logical. Are the variances equal? Default is \code{TRUE}.
#' @export linear.comparison
#' @author Rob Cribbie and Phil Chalmers 
#' @examples
#' \dontrun{
#' data(contrasts)
#' attach(contrasts)
#'
#' linear.comparison(achiev,major,c.weights=list(c(1,0,-.5,-.5,0),
#'  c(1,-2,0,0,1),c(1,0,-1,-1,1),c(0,1,1,-1,-1)),var.equal=F)
#'  }
linear.comparison <-
function(y,group,c.weights,alpha=0.05,var.equal=TRUE){
    #
    # Computes an F test for a linear comparison among group means
    # Uses procedures described in Maxwell & Delaney (2004; Designing Experiments
    # and Analyzing Data), pages 157-1162.
    #
    # If variances are assumed to differ across groups, the routine calculates
    # MS_within and adjusts df_within using procedures developed by Welch (1938) & Shatterwaite (1946),
    # as described in Maxwell & Delaney (2004; Designing Experiments and Analyzing
    # Data), pages 165-168.
    #
    # This routine is suitable for one.way designs.
    #
    # y: dependent variable
    # group: grouping variable
    #c.weights: a list containing weights for multiple linear contrasts; or a numeric vector containg weights for 1 contrast
    #alpha: compute 100*(1-alpha)% confidence interval; default is 0.05
    #
    #
    # EXAMPLE:
    ## assuming there are 4 groups:
    #linear.contrast(y=myScores,group=myGroups,c.weights=list(c(-1,-1,-1,3)),var.equal=FALSE)
    #
    #my.contrasts <- list(c(-1,-1,-1,3),c(-1,-1,2,0),c(-1,1,0,0) );
    #linear.contrast(y=myScores,group=myGroups,c.weights=my.contrasts,var.equal=FALSE)
    #
    #

if(var.equal==TRUE){
#print("computing linear comparisons assuming equal variances among groups")
}
if(var.equal==FALSE){
#print("computing linear comparisons assuming unequal variances among groups")
}
if(class(c.weights)=="numeric"){
c.weights <- list(c.weights);
}
cw.list <- as.list(c.weights);
n.c <- length(cw.list)
tmp<-list();
for (kk in 1:n.c){
lc<-cw.list[[kk]];
# print(lc)
if (abs(sum(lc))>0.001){
stop("the sum of contrast weights in lc does not equal zero")
}
group.mean <-tapply(y,group,mean); # group means
if (length(group.mean) != length(lc)){
stop("the contrast vector is not the same as the number of groups")
}

if(var.equal==TRUE){
tmp[[kk]] <- lc.var.equal(y,group,lc,alpha,C=n.c)
}
if(var.equal==FALSE){
tmp[[kk]] <- lc.var.unequal(y,group,lc,alpha,C=n.c)
}
}
class(tmp) <- 'linear.comparison'
tmp
}

#' @S3method print linear.comparison
#' @rdname linear.comparison
#' @method print linear.comparison
print.linear.comparison <- function(x){
  cat('\n\t')
  X <- x[[1]]
  if(class(X) == "lc.var.unequal") cat("Linear comparison(s) with unequal variance \n\n")
    else cat("Linear comparison(s) with equal variance \n\n")
  for(i in 1:length(x)){
    X <- x[[i]]  
    cat('Contrast:',X$contrast, '\n')
    cat('Alpha:',X$alpha, '\n')
    cat('F(',X$df1,',', round(X$df2,3),') = ', round(X$F,3),', p = ',
        round(X$p.2tailed,3), '\n',sep='')
    cat('Cohen\'s d =',round(X$d.effect.size,3), '\n')  
    cat('Adjusted confidence interval:', round(X$adj.confint,3), '\n')
    cat('\n')
  }
}
