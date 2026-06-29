#' Inverse symmetric square root
#'
#' Compute the inverse of the symmetric square root of a positive definite matrix.
#'
#' @param cov a positive definite matrix.
#'
#' @return the inverse of the symmetric square root of given a positive definite matrix \code{cov}.
#'
#' @author Bongjung Sung
#' @examples
#'
#' # generate a positive definite matrix
#' set.seed(100)
#' X=matrix(rnorm(16*32),ncol=16)
#' S=crossprod(X,X)/32
#' sym.inv.root(S)
#'
#' @export
sym.inv.root=function(cov){

  cov.eig=eigen(cov)
  Q=cov.eig$vectors; D=cov.eig$values
  return(Q%*%diag(1/sqrt(D))%*%t(Q))

}

#' Symmetric square root
#'
#' Compute the symmetric square root of a positive definite matrix.
#'
#' @param cov a positive definite matrix.
#'
#' @return the symmetric square root of given a positive definite matrix \code{cov}.
#'
#' @author Bongjung Sung
#' @examples
#'
#' # generate a positive definite matrix
#' set.seed(100)
#' X=matrix(rnorm(16*32),ncol=16)
#' S=crossprod(X,X)/32
#' sym.root(S)
#'
#' @export
sym.root=function(cov){

  cov.eig=eigen(cov)
  Q=cov.eig$vectors; D=cov.eig$values
  return(Q%*%diag(sqrt(D))%*%t(Q))

}

#' Sample covariance matrix of the data tensor.
#'
#' Compute the sample covariance matrix of a given data tensor.
#'
#' @param dat the \eqn{n \times p1 \times p2} data tensor.
#' @param center logical, whether to center the data or not.
#'
#' @return the sample covariance matrix of a given data tensor \code{dat}.
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' X=array(rnorm(96),dim=c(8,3,4))
#' dat2cov(X,center=TRUE)
#'
#' @export
dat2cov=function(dat,center){

  # dat: a n x p1 x p2 tensor
  # n: the sample size
  n=dim(dat)[1]

  # center: whether to center the data or not
  if(center==TRUE){
    dat=apply(dat,c(2,3),scale,scale=FALSE)
  }

  dat=aperm(dat,perm=c(2,3,1))
  dat=matrix(dat,ncol=n)
  return(tcrossprod(dat,dat)/n)

}

#' Parameters of the transformation to obtain the Tracy-Widom law
#'
#' Compute parameters associated with transforming the largest eigenvalue of the sample core to obtain \code{TW}-law.
#'
#' @param K1 the \eqn{p1 \times p1} covariance matrix.
#' @param K2 the \eqn{p2 \times p2} covariance matrix.
#' @param n the sample size.
#'
#' @return a vector of parameters that are associated with transforming the largest eigenvalue of the sample core to obtain \code{TW}-law.
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' p1=10; p2=10; n=400; p=p1*p2
#' X=matrix(rnorm(p*n),ncol=p)
#' S=crossprod(X,X)/n
#' S.kcd=covKCD::covKCD(S,p1,p2)
#' K1=S.kcd$K1
#' K2=S.kcd$K2
#' tw.para(K1,K2,n)
#'
#'@export
tw.para=function(K1,K2,n){

  p1=ncol(K1); p2=ncol(K2)
  p=p1*p2
  K1.eig=1/eigen(K1)$values
  K2.eig=1/eigen(K2)$values

  K.eig=as.numeric(kronecker(K1.eig,K2.eig))
  K.eig.max=max(K1.eig)*max(K2.eig)

  f=function(x){
    num=(sum((x*K.eig/(1-x*K.eig))^2)/p-n/p)^2
    return(num)
  }

  xi=optimize(f,lower=1e-03,upper=1/K.eig.max)$minimum
  E.plus=1/xi*(1+p/n*sum((xi*K.eig/(1-xi*K.eig))/p))
  gamma0=(p/n)*sum((K.eig/(1-xi*K.eig))^3)/p+1/(xi^3)
  gamma0=gamma0^(-1/3)

  return(c(gamma0,E.plus))

}

#' Parameters of the transformation for the extended LRT statistic
#'
#' Compute the parameters associated with transforming the extended LRT statistic based on the sample core.
#'
#' @param p1 the row dimension
#' @param p2 the column dimension
#' @param n the sample size
#'
#' @return a parameter that is associated with transforming the extended LRT statistic based on the sample core.
#'
#' @author Bongjung Sung
#' @examples
#'
#' p1=20; p2=16; n=400
#' elrt.para(p1,p2,n)
#'
#'@export
elrt.para=function(p1,p2,n){

  p=p1*p2; y=p/n
  a1=(y+2-sqrt(y^2+4))/(2*sqrt(y))
  if(y<1){
    int.log=-1/2*(2*a1/sqrt(y)+(1-y)/y*log(1-a1*sqrt(y))+(1-y)/y*log(sqrt(y)/a1-y)-(y+1)/y*log(sqrt(y)/a1))
    return(p*int.log)
  }else{
    int.log=-1/2*(2*a1/sqrt(y)+(y-1)/y*log(1-a1/sqrt(y))+(y-1)/y*log(sqrt(y)/a1-1)-(y+1)/y*log(sqrt(y)/a1))
    return(p*int.log)
  }

}

