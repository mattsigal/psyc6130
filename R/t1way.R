t1way <-
function(x,tr=.2,grp=NA){
#
#  A heteroscedastic one-way ANOVA for trimmed means
#  using a generalization of Welch's method.
#
#  The data are assumed to be stored in $x$ in list mode.
#  Length(x) is assumed to correspond to the total number of groups.
#  By default, the null hypothesis is that all groups have a common mean.
#  To compare a subset of the groups, use grp to indicate which
#  groups are to be compared. For example, if you type the
#  command grp<-c(1,3,4), and then execute this function, groups
#  1, 3, and 4 will be compared with the remaining groups ignored.
#
#  Missing values are automatically removed.
#
if(is.matrix(x)){
y<-list()
for(j in 1:ncol(x))y[[j]]<-x[,j]
x<-y
}
if(is.na(grp[1]))grp<-c(1:length(x))
if(!is.list(x))stop("Data are not stored in a matrix or in list mode.")
J<-length(grp)  # The number of groups to be compared
h<-vector("numeric",J)
w<-vector("numeric",J)
xbar<-vector("numeric",J)
for(j in 1:J){
xx<-!is.na(x[[j]])
val<-x[[j]]
x[[j]]<-val[xx]  # Remove missing values
h[j]<-length(x[[grp[j]]])-2*floor(tr*length(x[[grp[j]]]))
   # h is the number of observations in the jth group after trimming.
w[j]<-h[j]*(h[j]-1)/((length(x[[grp[j]]])-1)*winvar(x[[grp[j]]],tr))
xbar[j]<-mean(x[[grp[j]]],tr)
}
u<-sum(w)
xtil<-sum(w*xbar)/u
A<-sum(w*(xbar-xtil)^2)/(J-1)
B<-2*(J-2)*sum((1-w/u)^2/(h-1))/(J^2-1)
TEST<-A/(B+1)
nu1<-J-1
nu2<-1./(3*sum((1-w/u)^2/(h-1))/(J^2-1))
sig<-1-pf(TEST,nu1,nu2)
ret <- list(TEST=TEST,nu1=nu1,nu2=nu2,p.value=sig,trim=tr)
class(ret) <- 't1way'
ret 
}

print.t1way <- function(x){
  cat('\n\tWelch\'s one way ANOVA with trimmed means \n\n')
  cat('Proportion trimmed =', x$trim, '\n')
  cat('F(',x$nu1,',',round(x$nu2,3),') = ',round(x$TEST,3),', p = ', 
      round(x$p.value, 3), '\n', sep='') 
}
