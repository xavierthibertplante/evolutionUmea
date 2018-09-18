autosomalSelection<-function(AA,Aa,aa,p,it) {
#' Autosomal Selection
#' 
#' @param AA Relative fitness of the genotype AA (real between zero and one). Default 1.
#' @param Aa Relative fitness of the genotype Aa (real between zero and one).  Default 0.9.
#' @param aa Relative fitness of the genotype aa (real between zero and one). Default 0.8
#' @param p Initial frequency of A (real between zero and one). Default c(0.10 0.27 0.44 0.61 0.78 0.95).
#' @param it Number of generations (integer greater than zero). Default 150. 
#' @return A graph of allele frequency, a graph of change of allele frequency as function of allele frequency and, if only one initial frequency, a graph of genotype frequency.
#' @examples
#' autosomalSelection()
#' autosomalSelection(p=0.5)
#' autosomalSelection(p=c(0.1,0.9))
#' autosomalSelection(AA=0.5,Aa=1.0,aa=0.5,p=0.5,it=150)
#' @export
if(missing(AA)) AA<-1
if(missing(Aa)) Aa<-0.9
if(missing(aa)) aa<-0.8
if(missing(p)) p<-seq(from=0.1,to=0.95,length.out=6)
if(missing(it)) it<-150

par(mfrow=c(1,1))

OUT<-NULL
OUT<-rbind(OUT,p)

if (length(p) == 1 ) {
   OUT2<-NULL
   OUT2<-rbind(OUT2,c(p^2,2*p*(1-p),(1-p)^2))
}


for (i in 1:(it-1)) {
   q<-1-p
   p<-p* (p*AA+q*Aa)/((p^2)*AA+2*p*q*Aa+q^2*aa)
   OUT<-rbind(OUT,p)
   if (length(p) == 1 ) {
      OUT2<-rbind(OUT2,c(p^2,2*p*(1-p),(1-p)^2))
   }
}


if (length(p) == 1 ) {
   par(mfrow=c(3,1))
} else  {
   par(mfrow=c(2,1))
}
# p vs t
plot(x=1,y=1,col=0,xlim=c(1,it),ylim=c(0,1),main="Autosomal Selection",
xlab="Generations (t)", ylab="Allele A Frequency (p)")
for (i in 1:ncol(OUT)){
   points(x=1:it,y=OUT[,i],col=i)
}

# derivative
tmpp<-seq(from=0,to=1,by=0.01)
tmpq<-1-tmpp
dp<-tmpp* (tmpp*AA+tmpq*Aa)/((tmpp^2)*AA+2*tmpp*tmpq*Aa+tmpq^2*aa)  -tmpp
plot(tmpp,dp,xlab="Allelic Frequency of A (p)",ylab=expression(paste(Delta,"p")), main="Allelic Frequency Change")


if (length(p) == 1 ) {
   plot(x=1,y=1,col=0,xlim=c(1,it),ylim=c(0,1),main="Genotypic Frequency Trajectory",
   xlab="Generations (t)", ylab="Genotype Frequencies")
   for (i in 1:3) {
      points(x=1:it,y=OUT2[,i],col=i)
   }
   legend("right",legend=c(expression("p"^2),"2pq",expression("q"^2)),pch=19,col=c(1,2,3))
}



#return(OUT)

}
