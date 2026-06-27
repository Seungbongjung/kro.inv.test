#' Kronecker-Invariant test
#'
#'
#'
#'
#'
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @import RSpectra, covKCD
#' @export
large.eig.null=function(n,p1,p2,center=FALSE,trans=TRUE,sigma.known=FALSE,samp.num=1000,samp.iter=10,seed=100){

  gamma1=p1/sqrt(n); gamma2=p2/sqrt(n);

  set.seed(seed)

  for(i in 1:samp.iter){

    if(sigma.known==FALSE){

      if(trans==TRUE){

      }else{

      }

    }else{

    }

    if(i%%samp.iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  return(test.stat)
}
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @import RSpectra, covKCD
#' @export
large.eig.alt=function(n,p1,p2,core,center=FALSE,trans=TRUE,sigma.known=FALSE,samp.num=1000,seed=100){

}
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @import RSpectra
#' @export
large.eig.power=function(dat,center=TRUE,trans=TRUE,para.test=TRUE,samp.num=1000,seed=100){

}

