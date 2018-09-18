populationStructure<-function(N,p,m,it,Nd,selfing){
#' Autosomal Selection
#' 
#' @param N Population Size in each deme (integer greater than zero). Default 100.
#' @param p Initial frequency of A (real between zero and one). Default 0.5.
#' @param m Migration rate (real between zero and one). Default 0.01.
#' @param it Number of generations (integer greater than zero). Default 300. 
#' @param Nd Number of demes (integergreater than zero). Default 6.
#' @param selfing Is selfing possible (TRUE or FALSE).  Default FALSE.
#' @return A graph of allele frequency and a graph of F-statistics.
#' @examples
#' populationStructure()
#' populationStructure(p=0.5)
#' populationStructure(N=100,selfing=TRUE)
#' populationStructure(Nd=10,selfing=FALSE)
#' @export
if(missing(N)) N<-100
if(missing(p)) p<-0.5
if(missing(m)) m<-0.01
if(missing(it)) it<-300
if(missing(Nd)) Nd<-6
if(missing(selfing)) selfing<-F

initP<-p

par(mfrow=c(1,1))

#if(missing(L)) L<-1
L<-1

if (length(p) > 1) {
   if (length(p) != Nd) {
      p<-rep(p[1],Nd)
   }
}  else  {
   p<-rep(p,Nd)
}
OUT<-NULL
Fis<-NULL
Fst<-NULL
Fit<-NULL

pop<-array(dim=c(2*L,N,Nd))
for (i in 1:Nd) { 
pop[,,i]<-matrix(ifelse(runif(N*2*L)>p[i],0,1),ncol=N,nrow=2*L)
}
pop2<-pop

for (i in 1:it) {  
   tmpOUT<-rep(NA,Nd)
   tmpFis<-rep(NA,Nd)
   tmpFst<-rep(NA,Nd)
   tmpFit<-rep(NA,Nd)
   tmpHs<-rep(NA,Nd)
   tmpHi<-rep(NA,Nd)
#  Total heterozygosity observed
#   Ht_obs<-sum(pop[2*(1:L)-1,,]!=pop[2*(1:L),,])/(L*N*Nd)
#  Total heterozygosity expected
   tmpP<-sum(pop[,,])/(2*L*N*Nd)
   Ht<-2*tmpP*(1-tmpP)
   for (k in 1:Nd) { # For each populations
      # stats output frequency
      j<-1
      tmpOUT[k]<-mean(pop[(2*j):(2*j-1),,k])
      # F-stats 
      # Deme heterozygosity observed
      Hi<-sum(pop[2*(1:L)-1,,k]!=pop[2*(1:L),,k])/(L*N)
      #  Total heterozygosity expected
      Hs<-sum(pop[,,k])/(2*L*N)
      tmpHs[k]<-2*Hs*(1-Hs)
      tmpHi[k]<-Hi
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
      pop2[2*(1:L)-1,1:N,k]<-pop[2*(1:L)-sample(0:1,L,replace=T),parentsA,k]
      pop2[2*(1:L),1:N,k]  <-pop[2*(1:L)-sample(0:1,L,replace=T),parentsB,k]
   #   }
      pop[,,k]<-pop2[,,k] # replace the populations
   }
   OUT<-cbind(OUT,tmpOUT)
   Fis<-cbind(Fis,(mean(tmpHs,na.rm=T)-mean(tmpHi,na.rm=T))/mean(tmpHs,na.rm=T))
#   Fst<-cbind(Fst,tmpFst)
   Fst<-cbind(Fst,(Ht-mean(tmpHs,na.rm=T))/Ht)
   Fit<-cbind(Fit,(Ht-mean(tmpHi,na.rm=T))/Ht)
   # migration
   nMig<-rbinom(1, N*Nd, m/2) # m/2 we will move two individuals at the time
   if (nMig > 0) {
      for  (j in 1:nMig) {
         deme<-sample(1:Nd,2,replace=F)
         ind<-sample(1:N,2,replace=T)
         tmpIND<-pop[(2*(1:L)-1):(2*(1:L)),ind[1],deme[1]]
         pop[(2*(1:L)-1):(2*(1:L)),ind[1],deme[1]]<-
            pop[(2*(1:L)-1):(2*(1:L)),ind[2],deme[2]]
         pop[(2*(1:L)-1):(2*(1:L)),ind[2],deme[2]]<-tmpIND
      }
   }
}

par(mfrow=c(2,1))
plot(OUT[1,],xlab="Generations (t)",ylab="Allelic Frequency (p)",ylim=c(0,1),type='l',main=paste("Population Structure: Demic Allele Frequencies \n", "N=",N,"p=",initP,"m=",m,"it=",it,"Nd=",Nd,"selfing=",selfing))
if(k > 1) {
   for (i in 2:k) {
      lines(OUT[i,],col=i)
   }
}
Fis<-colMeans(Fis,na.rm=T)
Fst<-colMeans(Fst,na.rm=T)
Fit<-colMeans(Fit,na.rm=T)
Fis[is.infinite(Fis)]<-NA
Fst[is.infinite(Fst)]<-NA
Fit[is.infinite(Fit)]<-NA
Fis[is.nan(Fis)]<-NA
Fst[is.nan(Fst)]<-NA
Fit[is.nan(Fit)]<-NA
lowLim<-min(Fis,Fst,Fit,na.rm=T)
highLim<-max(Fis,Fst,Fit,na.rm=T)
plot(Fis,xlab="Generations (t)",ylab="F-statistics",
ylim=c(lowLim,highLim),
type='l',main="Population Structure: Wright's F_Statistics",col=2)
lines(Fst,col=3)
lines(Fit,col=4)
lines(rep(0,it),col=1,lty=2)
legend("topleft",legend=c(expression("F"["is"]),expression("F"["st"]),expression("F"["it"])),pch=19,col=c(2,3,4))

}
