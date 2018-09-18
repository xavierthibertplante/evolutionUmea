my_initVose1991 <-function(p) {
#' Calculate coefficient for Vose 1991 algorithm
#'
#' @param p A vector of positive real values
#' @return Nothing: internal function
#' @examples
#' my_initVose1991(p)
#' @references Michael D. Vose (1991) A linear algorithm for generating random numbers with a given distribution. _IEEE Transactions on Software Engineering_, *17*, 972-975.
#' @export
# p must be normalized!!!
# added normalization
#   global n prob aalias
   p=p/sum(p)
   n<-length(p)
   prob<-10000*rep(1,n)
   aalias<-10000*rep(1,n)
   large<-which(p>1/n )
   small<-which(p<=1/n)
   s=length(small)+1
   l=length(large)+1
   while ((s != 1) && (l != 1)) {
      s=s-1
      j=small[s]
      l=l-1
      k=large[l]
      prob[j]=n*p[j]
      aalias[j]=k
      p[k]=p[k]+(p[j]-(1.0/n))
      if (p[k] > 1.0/n) {
         large[l] = k 
         l=l+1
      } else {
         small[s] = k 
         s=s+1
      }
   }
   while  (s>1) {
      s=s-1
      prob[small[s]]=1
   }
   while  (l>1) {
      l=l-1
      prob[large[l]]=1
   }

   n<<-n
   prob<<-prob
   aalias<<-aalias


}
