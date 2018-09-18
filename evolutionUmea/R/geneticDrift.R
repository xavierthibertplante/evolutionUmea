geneticDrift<-function(N,L,p,it,selfing){
#' Genetic Drift
#' 
#' @param N Population size (integer greater than zero). Default 10
#' @param L Number of loci (integer greater than zero).  Default 1.
#' @param p Initial frequency of A (real between zero and one). Default p=0.5.
#' @param it Number of generations (integer greater than zero). Default 150. 
#' @param selfing If selfing allowed (TRUE or FALSE). Default FALSE.
#' @return A graph of allele frequency.
#' @examples
#' geneticDrift()
#' geneticDrift(p=0.5)
#' geneticDrift(N=100,L=3,p=0.3,it=500,selfing=TRUE)
#' @export
if(missing(N)) N<-10
if(missing(L)) L<-1
if(missing(p)) p<-0.5
if(missing(it)) it<-3*N
if(missing(selfing)) selfing<-F

OUT<-NULL

pop<-matrix(ifelse(runif(N*2*L)>p,0,1),ncol=N,nrow=2*L)
pop2<-pop

par(mfrow=c(1,1))

for (i in 1:it) {  
   # stats utput frequency
   tmpOUT<-rep(NA,L)
   for (j in 1:L) {
      tmpOUT[j]<-mean(pop[(2*j):(2*j-1),])
   }
   OUT<-cbind(OUT,tmpOUT)
   # time iteration
   # find parents
   parentsA<-sample(1:N,N,replace=T)
   parentsB<-sample(1:N,N,replace=T)
   if (!selfing) {
      # prevent selfing
      tmp<-parentsA==parentsB
      while(sum(tmp) !=0) {
         parentsB[tmp]<-sample(1:N,sum(tmp),replace=T)
         tmp<-parentsA==parentsB
      }
   }
   # combine into offspring
   pop2[2*(1:L)-1,1:N]<-pop[2*(1:L)-sample(0:1,L,replace=T),parentsA]
   pop2[2*(1:L),1:N]  <-pop[2*(1:L)-sample(0:1,L,replace=T),parentsB]
#   }
   pop<-pop2 # replace the populations
}

plot(OUT[1,],xlab="Generations (t)",ylab="Allelic Frequency of A (p)",ylim=c(0,1),type='l',main="Genetic Drift")
if(L > 1) {
   for (i in 2:L) {
      lines(OUT[i,],col=i)
   }
}
}
