evMagicComplete<-function(magic,sigma_s,sigma_a,K,b,mig,h) {
#' evMagicComplete
#' 
#' @param magic Magic trait (TRUE) or not (FALSE).  Default FALSE
#' @param sigma_s Strength of disruptive selection. Default 0.8
#' @param sigma_a Strength of non-random mating. Default 0.1
#' @param K Carying capacity. Default 512
#' @param b Average nuber of offspring per female. Default 4
#' @param mig Migration rate. Default 0.025
#' @param h Mutation effect size. Default 0.125
#' @return A multi panel graph of the evolution of naturally selected trait, neutral trait, preference trait, assortative mating metric, intial mating network, and final mating network.
#' @examples
#' evMagicComplete()
#' evMagicComplete(magic=TRUE)
#' evMagicComplete(magic=FALSE)
#' @export
require(igraph)
#set.seed(12345678)

if(missing(magic)) magic<-F
if(missing(sigma_s)) sigma_s<-0.8
if(missing(sigma_a)) sigma_a<-0.1
if(missing(K)) K<-512
if(missing(b)) b<-4
if(missing(mig)) mig<-0.025
if(missing(h)) h<-0.125


popS<-100 # initial population size
#sigma_s<-0.8
#sigma_a<-0.1
#K<-512
#b<-4
#mig<-0.025
#h<-0.125
#magic=F

x<-c(rep(0,popS/2),rep(0,popS/2))
c<-rep(0,popS)
y<-rep(0,popS)
env<-c(rep(0,popS/2),rep(1,popS/2))
bm1<-b-1

trackC0<-NULL
trackC1<-NULL
trackX0<-NULL
trackX1<-NULL
trackGamma0<-NULL
trackGamma1<-NULL
trackY0<-NULL
trackY1<-NULL

for (it in 1:300){
################### fitness
w<-exp((-(x-env)^2)/(2*sigma_s^2))
#################### survival
nu<-rep(NA,popS)
cnd<-env==0
nu[cnd]<-K*w[cnd]/((K*w[cnd])+bm1*sum(cnd))
cnd<-env==1
nu[cnd]<-K*w[cnd]/((K*w[cnd])+bm1*sum(cnd))
################### mortality
death<-runif(popS)
cnd<-(death<=nu)
x<-x[cnd]
c<-c[cnd]
env<-env[cnd]
y<-y[cnd]
popS<-length(env)
################### assortative mating
m<-matrix(NA,ncol=popS,nrow=popS)
colnames(m)=1:popS
rownames(m)=1:popS
n<-matrix(NA,ncol=popS,nrow=popS)
colnames(n)=1:popS
rownames(n)=1:popS
x1<-NULL
c1<-NULL
env1<-NULL
y1<-NULL
for (i in 1:popS){
   if (magic) {
      if (c[i] >= 0) {
         m[i,]<-exp((-c[i]^2*(x[i]-x)^2)/(2*sigma_a^2))
      } else {
         m[i,]<-2-exp((-c[i]^2*(x[i]-x)^2)/(2*sigma_a^2))
      }
   } else {
      if (c[i] >= 0) {
         m[i,]<-exp((-c[i]^2*(y[i]-y)^2)/(2*sigma_a^2))
      } else {
         m[i,]<-2-exp((-c[i]^2*(y[i]-y)^2)/(2*sigma_a^2))
      }
   }
   n[i,]<-m[i,]
   m[i,]<-m[i,]/(sum(m[i,]))
   m[i,]<-ifelse(env[i]==env,m[i,],0)
   tmp<-sample(m[i,],1,prob=m[i,])
   tmp2<-which(m[i,]==tmp)
   if (length(tmp2)>1) {
      f<-sample(tmp2,1) # father
   } else {
      f<-tmp2
   }
   for (j in 1:rpois(1,b)){
      x1<-c(x1,min(max(((x[i]+x[f])/2)+rnorm(1)*h,-100),100))
      c1<-c(c1,min(max(((c[i]+c[f])/2)+rnorm(1)*h,0),1))
      y1<-c(y1,min(max(((y[i]+y[f])/2)+rnorm(1)*h,-100),100))
      env1<-c(env1,env[i]) # kid close to the mother
   }
}
netX2<-x
netY2<-y
netC2<-c
netEnv2<-env
################### population switch
popS<-length(x1)
x<-x1
c<-c1
env<-env1
y<-y1
################### migration
lst<-sample(1:popS,popS*mig)
env[lst]<-ifelse(env[lst]==1,0,1)

# Calculate population measure
barx0<-mean(x[env==0])
barx1<-mean(x[env==1])
barc0<-mean(c[env==0])
barc1<-mean(c[env==1])
bary0<-mean(y[env==0])
bary1<-mean(y[env==1])
trackC0<-c(trackC0,barc0)
trackC1<-c(trackC1,barc1)
trackX0<-c(trackX0,barx0)
trackX1<-c(trackX1,barx1)
trackY0<-c(trackY0,bary0)
trackY1<-c(trackY1,bary1)
if (magic) {
   if (barc0 >= 0) { 
      psid<-exp((-barc0^2*(barx0-barx1)^2)/(2*sigma_a^2))
   } else {
      psid<-2-exp((-barc0^2*(barx0-barx1)^2)/(2*sigma_a^2))
   }
} else {
   if (barc0 >= 0) { 
      psid<-exp((-barc0^2*(bary0-bary1)^2)/(2*sigma_a^2))
   } else {
      psid<-2-exp((-barc0^2*(bary0-bary1)^2)/(2*sigma_a^2))
   }
}
psis<-1 # exp((-barc0^2*(barx0-barx0)^2)/(2*sigma_a^2))
gamma0<-psid/(psid+psis)
if (magic) {
   if (barc1 >= 0) { 
      psid<-exp((-barc1^2*(barx1-barx0)^2)/(2*sigma_a^2))
   } else {
      psid<-2-exp((-barc1^2*(barx1-barx0)^2)/(2*sigma_a^2))
   }
} else {
   if (barc1 >= 0) { 
      psid<-exp((-barc1^2*(bary1-bary0)^2)/(2*sigma_a^2))
   } else {
      psid<-2-exp((-barc1^2*(bary1-bary0)^2)/(2*sigma_a^2))
   }
}
psis<-1 # exp((-barc1^2*(barx1-barx1)^2)/(2*sigma_a^2))
gamma1<-psid/(psid+psis)
trackGamma0<-c(trackGamma0,gamma0)
trackGamma1<-c(trackGamma1,gamma1)
# plot the graph   
if(it==2) {
   sze<-length(n[,1])
   tmp<-sort(sample(1:sze,floor(sze/10),replace=F))
   p<-n[tmp,tmp]
   netX1<-netX2[tmp]
   netY1<-netY2[tmp]
   netC1<-netC2[tmp]
   netEnv1<-netEnv2[tmp]
   sze<-length(netX1)
   net=graph.adjacency(p,mode="undirected",weighted=TRUE,diag=FALSE)
   if (magic) {
      layout1<-cbind(netX1+(1-netC1)*rnorm(length(sze)),netX1+(1-netC1)*rnorm(sze))
   } else {
      layout1<-cbind(netY1+(1-netC1)*rnorm(length(sze)),netY1+(1-netC1)*rnorm(sze))
   }
   V(net)$color<-netEnv1+2
   V(net)$size<-10.0
}
}
# subsample
sze<-length(n[,1])
tmp<-sort(sample(1:sze,floor(sze/10),replace=F))
p<-n[tmp,tmp]
netX2<-netX2[tmp]
netY2<-netY2[tmp]
netC2<-netC2[tmp]
netEnv2<-netEnv2[tmp]
p<-as.matrix(p)
net2=graph.adjacency(p,mode="undirected",weighted=TRUE,diag=FALSE)
V(net2)$color<-netEnv2+2
V(net2)$size<-10.0
sze<-length(netX2)
if (magic) {
   layout2<-cbind(netX2+(1-netC2)*rnorm(length(sze)),netX2+(1-netC2)*rnorm(sze))
} else {
   layout2<-cbind(netY2+(1-netC2)*rnorm(length(sze)),netY2+(1-netC2)*rnorm(sze))
}

# plot the graph   
par(mfrow=c(2,3))
plot(trackX0,ylim=c(min(c(trackX0,trackX1)),max(trackX0,trackX1)),
main="Ecological trait",ylab="x",xlab="Time",col=2)
points(trackX1,col=3)

plot(trackY0,ylim=c(min(c(trackY0,trackY1)),max(trackY0,trackY1)),
main="Neutral trait",ylab="y",xlab="Time",col=2)
points(trackY1,col=3)

plot.igraph(net,vertex.label=NA,layout=layout1, edge.color="black",edge.width=E(net)$weight, main="Initial network")
#test<-dist(10-9*nInit)
#plot(hclust(test))

plot(trackC0,ylim=c(min(c(trackC0,trackC1)),max(trackC0,trackC1)),
main="Preference",ylab="c",xlab="Time",col=2)
points(trackC1,col=3)

plot(trackGamma0,ylim=c(0,1),
main="Assortative mating",ylab="Gamma",xlab="Time",col=2)
points(trackGamma1,col=3)

plot.igraph(net2,vertex.label=NA,layout=layout2, edge.color="black",edge.width=E(net2)$weight,main="Final network")

#test<-dist(10-9*p)
#plot(hclust(test))


} # function
