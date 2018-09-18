saveFigure<-function(fileN){
#' Save Figure
#' 
#' @param fileN File name witout extension. Default "default".
#' @return Save a figure uder the current directory. If the file name exist an _k will be added and k will increment by one on each call of saveFigure.
#' @examples
#' saveFigure()
#' saveFigure("populationStructure")
#' saveFigure("inbreedingCoefficient")
#' saveFigure("test")
#' @export

if(missing(fileN)) fileN<-"default"

if (file.exists(paste(fileN,".jpg",sep=""))) {
   k<-1
   while (file.exists(paste(fileN,"_",k,".jpg",sep=""))) {
      k<-k+1
   }
   fileN<-paste(fileN,"_",k,".jpg",sep="")
} else {
   fileN<-paste(fileN,".jpg",sep="")
}

dev.copy(jpeg,fileN)
dev.off()

}
