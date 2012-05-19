#' Tukey and Games-Howell tests
#'
#' Tukey and Games-Howell post- hoc tests.
#' 
#' @aliases tukey
#' @param data the data
#' @param group logical; remove missing data?
#' @param method either \code{'Tukey'} or \code{'Games-Howell'}
#' @export tukey
#' @author Rob Cribbie and Phil Chalmers 
#' @examples
#' \dontrun{
#' data(nonnorm_hetvar)
#' attach(nonnorm_hetvar)
#' 
#' tukey(depres, Sex, method = 'Games-Howell')
#' }
tukey <- function(data, group, method=c("Tukey", "Games-Howell"))      
{
        OK <- complete.cases(data, group)                       
        data <- data[OK]
        group <- factor(group[OK])
        n <- tapply(data, group, length)                        
        a <- length(n)                                          
        phi.e <- sum(n)-a                                       
        Mean <- tapply(data, group, mean)                        
        Variance <- tapply(data, group, var)                     
        result1 <- cbind(n, Mean, Variance)                      
        rownames(result1) <- paste("Group", 1:a, sep="")
        method <- match.arg(method)
        if (method == "Tukey") {                                 
                v.e <- sum((n-1)*Variance)/phi.e                 
                t <- combn(a, 2, function(ij)                    
                                        abs(diff(Mean[ij]))/sqrt(v.e*sum(1/n[ij])) )
                p <- ptukey(t*sqrt(2), a, phi.e, lower.tail=FALSE)       
                Tukey <- cbind(t, p)                                     
                rownames(Tukey) <- combn(a, 2, paste, collapse=":")
                return(list(result1=result1, Tukey=Tukey, phi=phi.e, v=v.e))
        }
        else {                                                   
                t.df <- combn(a, 2, function(ij) {               
                                        t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
                                        df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
                                        return(c(t, df))} )
                t <- t.df[1,]
                df <- t.df[2,]
                p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)  
                Games.Howell <- cbind(t, df, p)                  
                rownames(Games.Howell) <- combn(a, 2, paste, collapse=":")
                return(list(result1=result1, Games.Howell=Games.Howell))
        }
}

