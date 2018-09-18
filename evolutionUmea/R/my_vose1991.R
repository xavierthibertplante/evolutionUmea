my_vose1991<-function(N) {
#' Return randomdeviates from Vose 1991 algorithm
#'
#' @param N Number of random deviate returned
#' @return Random deviate
#' @examples
#' p<-c(0.1,0.2,0.3,0.4)
#' my_initVose1991(p)
#' my_vose1991(4)
#' @references Michael D. Vose (1991) A linear algorithm for generating random numbers with a given distribution. _IEEE Transactions on Software Engineering_, *17*, 972-975.
#' @export
#   global n prob aalias
   myOut=rep(0,N)
   u=n*runif(N)+1
   j=floor(u)
   
   Z=u-j-prob[j]
   z2=which(Z>0)
   z1=which(Z<=0)
   myOut[z2]=aalias[j[z2]]
   myOut[z1]=j[z1]
  
   return(myOut)
}
