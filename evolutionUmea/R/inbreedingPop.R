inbreedingPop<-function(pop){
#' @export
   NN<-length(as.vector(pop))/2
   lst<-unique(as.vector(pop))
   tt<-0
   for (i in 1:length(lst)) {
      x<-sum(pop==lst[i])
      tt<-tt+((x/(2*NN))*((x-1)/((2*NN)-1)))
   }
   return(tt)
}
