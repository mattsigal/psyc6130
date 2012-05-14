lc.var.equal <-
function(y,group,lc,alpha=0.05,C=1){

# Computes an F test for a linear comparison among group means
# assuming that variances do NOT differ among the groups.
# Uses procedures described in Maxwell & Delaney (2004; Designing Experiments
# and Analyzing Data), pages 157-1162.
#
# This routine is suitable for one.way designs.
#
# y: dependent variable
# group: grouping variable
#lc: weights for linear contrast
#alpha: compute 100*(1-alpha)% confidence interval; default is 0.05
#C: number of comparisons

if (abs(sum(lc))>0.001){
stop("the sum of contrast weights in lc does not equal zero")
}



c.2<-lc^2; # square the weights
group.n<-tapply(y,group,length); # group n's
group.mean <-tapply(y,group,mean); # group means
y.sd<-tapply(y,group,sd); # group standard deviations
y.sd.2<-y.sd^2; # sd squared

psi <- sum(lc*group.mean); # Equation 40; chapter 4
tmp.lm<-lm(y~group);
df2 <- tmp.lm$df.residual
SS.within <- sum(residuals(tmp.lm)^2);
MS.within <- SS.within/df2; # Eq 41; chapter 3

tmp.lm<-lm(y~1);
SS.total <- sum(residuals(tmp.lm)^2);
SS.between <- SS.total-SS.within;

df1 <- 1;
F <- (psi*psi) / ( MS.within * sum(c.2/group.n) ) # Eq. 41; chapter 4

F.crit <- qf((1-alpha),df1,df2);
psi.low <- psi - sqrt(F.crit)*sqrt(sum( (c.2/group.n)*y.sd.2 ));
psi.high <- psi + sqrt(F.crit)*sqrt(sum( (c.2/group.n)*y.sd.2 ));

F.w <-qf( (1-alpha/C),df1,df2);
psi.low.adj <- psi - sqrt(F.w)*sqrt(MS.within*sum(c.2/group.n) );
psi.high.adj <- psi + sqrt(F.w)*sqrt(MS.within*sum(c.2/group.n) );


SS.contrast <- (psi*psi)/ sum( (lc^2)/group.n);

d.effect.size <- (2*psi) / (sqrt(MS.within)*sum(abs(lc))); # Eq 4-52
R2.alerting <- SS.contrast / SS.between; # Eq 4-54
R2.effect.size <- SS.contrast / SS.total; # Eq 4-55
R2.contrast <- SS.contrast / (SS.contrast+SS.within); # Eq 4-56

t <- sqrt(F);
if(psi<0){t <- -1*sqrt(F)};

ret <- list(contrast=lc, F=F, t=t, df1=df1, df2=df2, p.2tailed=1-pf(F,df1,df2), psi=psi, 
	confinterval=c(psi.low,psi.high),adj.confint=c(psi.low.adj,psi.high.adj), alpha=alpha, 
	SS.contrast=SS.contrast, d.effect.size=d.effect.size, R2.alerting=R2.alerting, 
	R2.effect.size=R2.effect.size, R2.contrast=R2.contrast)
class(ret) <- 'lc.var.equal'	
ret
}

