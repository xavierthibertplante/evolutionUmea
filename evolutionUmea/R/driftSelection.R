driftSelection<-function(N,p,it,AA,Aa,aa){
#' Drift Selection
#' 
#' @param N Population size (integer greater than zero). Default 100.
#' @param p Initial frequency of A (real between zero and one). Default 0.5.
#' @param AA Relative fitness of the genotype AA (real between zero and one). Default 0.5.
#' @param Aa Relative fitness of the genotype Aa (real between zero and one).  Default 1.0.
#' @param aa Relative fitness of the genotype aa (real between zero and one). Default 0.5
#' @param it Number of generations (integer greater than zero). Default 500. 
#' @return A graph of allele frequency.
#' @examples
#' driftSelection()
#' driftSelection(p=0.5)
#' driftSelection(AA=0.5,Aa=1.0,aa=0.5,p=0.5,it=150)
#' @export
if(missing(N)) N<-100
if(missing(p)) p<-0.5
if(missing(it)) it<-500
if(missing(AA)) AA<-0.8
if(missing(Aa)) Aa<-1.0
if(missing(aa)) aa<-0.9
L<-1
selfing=T

OUT<-NULL

pop<-matrix(ifelse(runif(N*2*L)<p,0,1),ncol=N,nrow=2*L)
pop2<-pop

par(mfrow=c(1,1))

for (i in 1:it) {  
   # stats output frequency
   tmpOUT<-rep(NA,L)
   for (j in 1:L) {
      tmpOUT[j]<-mean(pop[(2*j):(2*j-1),])
   }
   OUT<-cbind(OUT,tmpOUT)
   # time iteration
   # calculate fitness 
   fitT<-colSums(pop)
   fit<-rep(NA,N)
   fit[fitT==0]<-AA
   fit[fitT==1]<-Aa
   fit[fitT==2]<-aa
   my_initVose1991(fit)
   # find parents
   parentsA<-my_vose1991(N)
   parentsB<-my_vose1991(N)
   if (!selfing) {
      # prevent selfing
      tmp<-parentsA==parentsB
      while(sum(tmp) !=0) {
         parentsB[tmp]<-my_vose1991(sum(tmp))
         tmp<-parentsA==parentsB
      }
   }
   # combine into offspring
   pop2[2*(1:L)-1,1:N]<-pop[2*(1:L)-sample(0:1,L,replace=T),parentsA]
   pop2[2*(1:L),1:N]  <-pop[2*(1:L)-sample(0:1,L,replace=T),parentsB]
#   }
   pop<-pop2 # replace the populations
}

plot((1-OUT[1,]),xlab="Generations (t)",ylab="Allelic Frequency of A (p)",ylim=c(0,1),type='l',main=paste("Drift and Selection \n", "N=",N,"p=",p,"it=",it,"AA=",AA,"Aa=",Aa,"aa=",aa))
if(L > 1) {
   for (i in 2:L) {
      lines((1-OUT[i,]),col=i)
   }
}
}
