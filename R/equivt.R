#' Equivalence Test
#'
#' Test for the equivalence of two groups with a given tolerable interval. If a significant effect is returned 
#' for both t-values then there is evidence of equivalence. Allows for equal or unequal variances.
#'
#' @section \code{equivint} must be specified in the metric of the variables (i.e. not as a standard 
#' deviation, percent different, etc).
#' @aliases equivt
#' @param x Numeric vector for group 1.
#' @param y Numeric vector for group 2.
#' @param equivint Equivalence interval.
#' @param alpha Alpha value
#' @param varequiv Logical. Are there equal variances? Default is \code{FALSE}.
#' @param na.rm logical; remove missing data?
#' @param ... Additional arguments to be passed to the function.
#' @export equivt
#' @author Rob Cribbie and Phil Chalmers 
#' @examples
#' \dontrun{
#' data(nonnorm_hetvar)
#' attach(nonnorm_hetvar)
#' 
#' #equivalence within 2 points
#' eq <- 2
#' equivt(depres[Sex=="female"],depres[Sex=="male"], eq)
#' }
equivt <- function (x,y, equivint, alpha=.05,varequiv=FALSE, na.rm=TRUE, ...) {
   if (na.rm)    x <- x[!is.na(x)]
   if (na.rm)    y <- y[!is.na(y)]
   if (varequiv==FALSE) {
      t1<-(mean(x)-mean(y)-equivint)/sqrt((var(x)/length(x))+ (var(y)/length(y)))
      t2<-(mean(x)-mean(y)+equivint)/sqrt((var(x)/length(x))+ (var(y)/length(y)))
      dft<-(((var(x)/length(x))+(var(y)/length(y)))^2)/((var(x)^2/(length(x)^2*(length(x)-1)))+(var(y)^2/(length(y)^2*(length(y)-1))))
      probt1<-pt(t1,dft,lower.tail=TRUE)
      probt2<-pt(t2,dft,lower.tail=FALSE) 
      ifelse (probt1<alpha & probt2<alpha,
         decis<-"The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected",
         decis<-"The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
       }
    if (varequiv==TRUE) {
      t1<-(mean(x)-mean(y)-equivint)/sqrt(((((length(x)-1)*var(x)^2)+((length(y)-1)*var(y)^2))/(length(x)+length(y)-2))*(1/length(x)+1/length(y)))
      t2<-(mean(x)-mean(y)+equivint)/sqrt(((((length(x)-1)*var(x)^2)+((length(y)-1)*var(y)^2))/(length(x)+length(y)-2))*(1/length(x)+1/length(y)))
      dft<-length(x)+length(y)-2
      probt1<-pt(t1,dft,lower.tail=TRUE)
      probt2<-pt(t2,dft,lower.tail=FALSE)
      ifelse (probt1<alpha & probt2<alpha,
         decis<-"The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected",
         decis<-"The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
       }
       title<-"Schuirmann's Test of the Equivalence of Two Independent Groups"
       ei<-(c(equivint))
       names(ei)<-c("equivalence interval")
       tstats<-c(t1,t2)
       dfs<-c(dft,dft)
       pvals<-c(probt1,probt2)
       names(tstats)<-c("t1","t2")
       names(dfs)<-c("dft1","dft2")
       names(pvals)<-c("p_t1", "p_t2")
       out<-list (title,ei,tstats,dfs,pvals,decis)
	   class(out) <- 'equivt'
       out
}

#' @S3method print equivt
#' @rdname equivt
#' @method print equivt
print.equivt <- function(x){
  cat('\n\t',x[[1]],'\n\n',sep='')   
  cat('Equivalence interval:',x[[2]],'\n')
  tab <- data.frame(matrix(0,2,1))
  tab$t <- x[[3]]
  tab$df <- x[[4]]
  tab$p <- x[[5]]
  tab <- tab[,-1]
  row.names(tab) <- c('x', 'y')
  print(round(tab,3))
  cat('\n',x[[6]],'\n',sep='')
}

