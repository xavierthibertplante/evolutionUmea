inbreedingCoefficient<-function(N,it,initF){
#' Inbreeding Coefficient
#' 
#' @param N Population size (intteger greater than zero). Default 100.
#' @param initF Initial inbreeding coefficient (real between zero and one). Default 0.
#' @param it Number of generations (integer greater than zero). Default 50. 
#' @return A graph of inbreeding coefficient. Ft is the theoretical expectation, Fa is the individual autozygosity, Fp is the population level autozygosity
#' @examples
#' inbreedingCoefficient()
#' inbreedingCoefficient(N=25)
#' inbreedingCoefficient(N=25,it=500)
#' inbreedingCoefficient(N=30,initF=0.6)
#' @export
if(missing(N)) N<-100
if(missing(it)) it<-50
if(missing(initF)) initF<-0
# p73  of populus book
L<-1
selfing=F

OUT<-NULL

par(mfrow=c(1,1))

pop<-matrix(1:(N*2*L),ncol=N,nrow=2*L)
if(initF !=  0 ) {
   x<-initF*2*N
   pop[1:2,1:floor(x/2)] <- (-1)
   tt<-inbreedingPop(pop)
   k<-floor(x/2)+1
   while ( (tt < initF) & (k<=(N))) {
      pop[1,k] <- (-1)
      tt<-inbreedingPop(pop)
      k<-k+1
      tt
   } 
   tMin=log(1-initF)/log(1-(1/(2*N)))
   tMax=tMin+it-1
   F_t<-1-(1-(1/(2*N)))^(tMin:tMax)
} else {
   F_t<-1-(1-(1/(2*N)))^(1:it)
}
pop2<-pop

F_a<-1:it # place holder
F_p<-1:it # place holder


for (i in 1:it) {  
   # stats utput frequency
   tmpOUT<-rep(NA,2)
   # individual autozygosity
   tmpOUT[1]<-sum(pop[1,]==pop[2,])/N
   # population autozygosity
   tmpOUT[2]<-inbreedingPop(pop)
   OUT<-cbind(OUT,tmpOUT)
   # time iteration
   for (j in 1:N) { # replace each individuals
      # find parents
      x<-sample(1:N,2,replace=selfing)
      mother<-x[1]
      father<-x[2]
      # generate gametes
      chM<-pop[2*(1:L)-sample(0:1,L,replace=T),mother]
      chF<-pop[2*(1:L)-sample(0:1,L,replace=T),father]
      # combine into offspring
      pop2[2*(1:L)-1,j]<-chM
      pop2[2*(1:L),j]<-chF
   }
   pop<-pop2 # replace the populations
}

plot(F_t,xlab="Generations (t)",ylab="Inbreeding Coefficient",ylim=c(0,1),type='l',main="Inbreeding Coefficient vs Time",col=2)
lines(OUT[1,],col=3)
lines(OUT[2,],col=4)
legend("bottomright",legend=c(expression("F"["t"]),expression("F"["a"]),expression("F"["p"])),pch=19,col=c(2,3,4))


}
