#' Trimmed means for 2-way ANOVA
#'
#' A heteroscedastic two-way ANOVA for trimmed means using a generalization of Welch's method.
#' 
#' @section By default, 
#'
#' @aliases t2way
#' @param J number of levels for factor 1
#' @param K number of levels for factor 2
#' @param data contains the raw data stored in list mode, or a matrix with columns corresponding to groups.
#'  If stored in list mode, data[[1]] contains the data for the first level of both factors: level 1,1.
#'  data[[2]] is assumed to contain the data for level 1 of the first factor and level 2 of the second 
#'  factor: level 1,2, etc. It is assumed that data has length JK, the total number of groups being tested.
#' @param tr proportion to be trimmed, default is .2
#' @param grp \code{1:p}
#' @param p is a numberic value J*K
#' @param MAT logical; if \code{TRUE} assumes data are stored in matrix with 3 columns with two of the
#'  columns indicated by the argument \code{lev.col} specifying the columns of x containing the values of the
#'  levels of the two factors.
#' @param lev.col something
#' @param var.col something
#' @param pr logical; print the results?
#' @export t2way
#' @author Rob Cribbie and Phil Chalmers 
#' @examples
#' \dontrun{
#' data(twowaythreeway)
#' attach(twowaythreeway)
#'                 
#' }
t2way<-function(J,K,data,tr=.2,grp=c(1:p),p=J*K,MAT=F,
lev.col=c(1:2),var.col=3,pr=T){
    #  Perform a J by K  (two-way) anova on trimmed means where
    #  all groups are independent.
    #
    #  The R variable data is assumed to contain the raw
    #  data stored in list mode, or a matrix with columns
    #  corresponding to groups. If stored in list mode, data[[1]] contains the data
    #  for the first level of both factors: level 1,1.
    #  data][2]] is assumed to contain the data for level 1 of the
    #  first factor and level 2 of the second factor: level 1,2
    #
    #  The default amount of trimming is tr=.2
    #
    #  It is assumed that data has length JK, the total number of
    #  groups being tested.
    #
    #  MAT=T, assumes data are stored in matrix with 3 columns
    #  with two of the columns indicated by the argument
    #  lev.col
    #  specifying the columns of x containing the values of the
    #  levels of the two factors.
    #  The outcome variable is in column
    #  var.col
    #  which defaults to column 3
    #  That is, this function calls selby2 for you.
    #   
if(is.data.frame(data))data=as.matrix(data)
if(tr==.5){
print("For medians, use med2way if there are no ties")
print("With ties, use linear contrasts in conjunction with medpb")
stop("")
}
if(MAT){
if(!is.matrix(data))stop("With MAT=T, data must be a matrix")
if(length(lev.col)!=2)stop("Argument lev.col should have 3 values")
temp=selby2(data,lev.col,var.col)
lev1=length(unique(temp$grpn[,1]))
lev2=length(unique(temp$grpn[,2]))
gv=apply(temp$grpn,2,rank)
gvad=10*gv[,1]+gv[,2]
grp=rank(gvad)
if(pr){
print(paste("Factor 1 has", lev1, "levels"))
print(paste("Factor 2 has", lev2, "levels"))
}
if(J!=lev1)warning("J is being reset to the number of levels found")
if(K!=lev2)warning("K is being reset to the number of levels found")
J=lev1
K=lev2
data=temp$x
}
if(is.matrix(data))data=listm(data)
if(!is.list(data))stop("Data are not stored in list mode")
if(p!=length(data)){
print("The total number of groups, based on the specified levels, is")
print(p)
print("The number of groups in data is")
print(length(data))
print("Warning: These two values are not equal")
}
tmeans<-0
h<-0
v<-0
for (i in 1:p){
data[[grp[i]]]=elimna(data[[grp[i]]])
tmeans[i]<-mean(data[[grp[i]]],tr)
h[i]<-length(data[[grp[i]]])-2*floor(tr*length(data[[grp[i]]]))
#    h is the effective sample size
v[i]<-(length(data[[grp[i]]])-1)*winvar(data[[grp[i]]],tr)/(h[i]*(h[i]-1))
#    v contains the squared standard errors
}
v<-diag(v,p,p)   # Put squared standard errors in a diag matrix.
ij<-matrix(c(rep(1,J)),1,J)
ik<-matrix(c(rep(1,K)),1,K)
jm1<-J-1
cj<-diag(1,jm1,J)
for (i in 1:jm1)cj[i,i+1]<-0-1
km1<-K-1
ck<-diag(1,km1,K)
for (i in 1:km1)ck[i,i+1]<-0-1
#  Do test for factor A
#cmat<-kron(cj,kron(ik,il))  # Contrast matrix for factor A
cmat<-kron(cj,ik)  # Contrast matrix for factor A
#alval<-c(999:1)/1000
alval<-c(1:999)/1000
for(i in 1:999){
irem<-i
Qa<-johan(cmat,tmeans,v,h,alval[i])
if(Qa$teststat>Qa$crit)break
}
A.p.value=irem/1000
# Do test for factor B
cmat<-kron(ij,ck)  # Contrast matrix for factor B
for(i in 1:999){
irem<-i
Qb<-johan(cmat,tmeans,v,h,alval[i])
if(Qb$teststat>Qb$crit)break
}
B.p.value=irem/1000
# Do test for factor A by B interaction
cmat<-kron(cj,ck)  # Contrast matrix for factor A by B
for(i in 1:999){
irem<-i
Qab<-johan(cmat,tmeans,v,h,alval[i])
if(Qab$teststat>Qab$crit)break
}
AB.p.value=irem/1000
tmeans=matrix(tmeans,J,K,byrow=T)
list(Qa=Qa$teststat,A.p.value=A.p.value,
Qb=Qb$teststat,B.p.value=B.p.value,
Qab=Qab$teststat,AB.p.value=AB.p.value,means=tmeans)
}

selby2<-function(m,grpc,coln=NA){
# Create categories according to the grpc[1] and grpc[2] columns
# of the matrix m. The function puts the values in column coln into
# a vector having list mode.
#
if(is.na(coln))stop("The argument coln is not specified")
if(length(grpc)>4)stop("The argument grpc must have length less than or equal to 4")
x<-vector("list")
ic<-0
if(length(grpc)==2){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2])]==c(cat1[i1],cat2[i2]))==2){
it<-it+1
temp[it]<-m[i,coln]
}
}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2]),1,2)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2]))
}
}}
}
if(length(grpc)==3){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
cat3<-selby(m,grpc[3],coln)$grpn
x<-vector("list")
ic<-0
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
for (i3 in 1:length(cat3)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2],grpc[3])]==c(cat1[i1],cat2[i2],cat3[i3]))==3){
it<-it+1
temp[it]<-m[i,coln]
}}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3]),1,3)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3]))
}}}}
}
if(length(grpc)==4){
cat1<-selby(m,grpc[1],coln)$grpn
cat2<-selby(m,grpc[2],coln)$grpn
cat3<-selby(m,grpc[3],coln)$grpn
cat4<-selby(m,grpc[4],coln)$grpn
x<-vector("list")
ic<-0
for (i1 in 1:length(cat1)){
for (i2 in 1:length(cat2)){
for (i3 in 1:length(cat3)){
for (i4 in 1:length(cat4)){
temp<-NA
it<-0
for (i in 1:nrow(m)){
if(sum(m[i,c(grpc[1],grpc[2],grpc[3],grpc[4])]==c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))==4){
it<-it+1
temp[it]<-m[i,coln]
}}
if(!is.na(temp[1])){
ic<-ic+1
x[[ic]]<-temp
if(ic==1)grpn<-matrix(c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]),1,4)
if(ic>1)grpn<-rbind(grpn,c(cat1[i1],cat2[i2],cat3[i3],cat4[i4]))
}}}}}
}
list(x=x,grpn=grpn)
}

johan<-function(cmat,vmean,vsqse,h,alpha=.05){
#
#  This function is used by other functions that come with this book,
#  and it can be used to test hypothesis not covered in the text.
#
#  The function performs Johansen's test of C mu = 0 for p independent groups,
#  where C is a k by p matrix of rank k and mu is a p by 1 matrix of
#  of unknown trimmed means.
#  The argument cmat contains the matrix C.
#  vmean is a vector of length p containing the p trimmed means
#  vsqe is a diagonal matrix containing the squared standard errors of the
#  the trimmed means in vmean.
#  h is a vector containing the effective sample sizes
#
yvec<-matrix(vmean,length(vmean),1)
if(!is.matrix(vsqse))vsqse<-diag(vsqse)
test<-cmat%*%vsqse%*%t(cmat)
invc<-solve(test)
test<-t(yvec)%*%t(cmat)%*%invc%*%cmat%*%yvec
R<-vsqse%*%t(cmat)%*%invc%*%cmat
A<-sum(diag((diag(R))^2/diag(h-1)))
df<-nrow(cmat)
crit<-qchisq(1-alpha,df)
crit<-crit+(crit/(2*df))*A*(1+3*crit/(df+2))
list(teststat=test[1],crit=crit[1])
}

winvar<-function(x,tr=.2,na.rm=FALSE){
#
#  Compute the gamma Winsorized variance for the data in the vector x.
#  tr is the amount of Winsorization which defaults to .2.
#
if(na.rm)x<-x[!is.na(x)]
y<-sort(x)
n<-length(x)
ibot<-floor(tr*n)+1
itop<-length(x)-ibot+1
xbot<-y[ibot]
xtop<-y[itop]
y<-ifelse(y<=xbot,xbot,y)
y<-ifelse(y>=xtop,xtop,y)
winvar<-var(y)
winvar
}

selby<-function(m,grpc,coln){
#
#
#  A commmon situation is to have data stored in an n by p matrix where
#  one or more of the columns are  group identification numbers.
#  This function groups  all values in column coln according to the
#  group numbers in column grpc and stores the  results in list mode.
#
#  More than one column of data can sorted
#
# grpc indicates the column of the matrix containing group id number
#
if(is.null(dim(m)))stop("Data must be stored in a matrix or data frame")
if(is.na(grpc[1]))stop("The argument grpc is not specified")
if(is.na(coln[1]))stop("The argument coln is not specified")
if(length(grpc)!=1)stop("The argument grpc must have length 1")
x<-vector("list")
grpn<-sort(unique(m[,grpc]))
it<-0
for (ig in 1:length(grpn)){
for (ic in 1:length(coln)){
it<-it+1
flag<-(m[,grpc]==grpn[ig])
x[[it]]<-m[flag,coln[ic]]
}}
list(x=x,grpn=grpn)
}

elimna<-function(m){
#
# remove any rows of data having missing values
#
if(is.null(dim(m)))m<-as.matrix(m)
ikeep<-c(1:nrow(m))
for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
elimna<-m[ikeep[ikeep>=1],]
elimna
}

kron<-function(m1,m2){
#  compute the Kronecker product of the two matrices m1 and m2.
#
m1<-as.matrix(m1) # Vectors of length p are converted to a p by 1 matrix
m2<-as.matrix(m2)
kron<-vector(mode="numeric",length=0)
for(i in 1:nrow(m1)){
m3<-m1[i,1]*m2
for(j in 2:ncol(m1))m3<-cbind(m3,m1[i,j]*m2)
if(i==1)kron<-m3
if(i>=2)kron<-rbind(kron,m3)
}
kron
}

t1way<-function(x,tr=.2,grp=NA){
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
list(TEST=TEST,nu1=nu1,nu2=nu2,p.value=sig)
}

