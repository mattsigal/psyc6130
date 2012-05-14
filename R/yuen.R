yuen <-
function(x,y,tr=.2,alpha=.05){
#
#  Perform Yuen's test for trimmed means on the data in x and y.
#  The default amount of trimming is 20%
#  Missing values (values stored as NA) are automatically removed.
#
#  A confidence interval for the trimmed mean of x minus the 
#  the trimmed mean of y is computed and returned in yuen$ci.
#  The p-value is returned in yuen$p.value
#
#  For an omnibus test with more than two independent groups, 
#  use t1way.
#  This function uses winvar from chapter 2.
#
if(tr==.5)stop("Using tr=.5 is not allowed; use a method designed for medians")
if(tr>.25)print("Warning: with tr>.25 type I error control might be poor")
x<-x[!is.na(x)]  # Remove any missing values in x 
y<-y[!is.na(y)]  # Remove any missing values in y
h1<-length(x)-2*floor(tr*length(x))
h2<-length(y)-2*floor(tr*length(y))
q1<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
q2<-(length(y)-1)*winvar(y,tr)/(h2*(h2-1))
df<-(q1+q2)^2/((q1^2/(h1-1))+(q2^2/(h2-1)))
crit<-qt(1-alpha/2,df)
dif<-mean(x,tr)-mean(y,tr)
low<-dif-crit*sqrt(q1+q2)
up<-dif+crit*sqrt(q1+q2)
test<-abs(dif/sqrt(q1+q2))
yuen<-2*(1-pt(test,df))
ret <- list(ci=c(low,up),p.value=yuen,dif=dif,se=sqrt(q1+q2),teststat=test,crit=crit,df=df,trim=tr,
            M1 = mean(x,tr), M2 = mean(y,tr))
class(ret) <- 'yuen'
ret
}

print.yuen <- function(x){
  cat('\n\tYuen\'s t-test \n\n')
  cat('Proportion trimmed =', x$trim, '\n')
  cat('Trimmed means =', round(x$M1,3), round(x$M2,3), '\n')
  cat('t = ',round(x$teststat,3),', df = ', round(x$df,3), ', p = ', round(x$p.value,3),  '\n',sep='')
  cat('Confidence interval =', round(x$ci,5))
}