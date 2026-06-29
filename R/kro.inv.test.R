#' Test statistic based on the largest eigenvalue of the sample core under the null
#'
#' Implement the Monte-Carlo simulations of the test statistic based on the largest eigenvalue of the sample core \eqn{\hat{C}} for Gaussian populations under the null hypothesis of separability.
#' Here the sample core \eqn{\hat{C}} refers to the core component of the sample covariance matrix \eqn{S}.
#' The test statistic may be transformed to approximate the Tracy-Widom (\code{TW}) law. If the population covariance matrix under the null is assumed to be known, the transformation may be based on the kronecker MLE (the separable component of \eqn{S}).
#' Otherwise, the transformation is based on the quantities depending only on \eqn{(n,p_1,p_2)}. For the details, see (20) of Sung and Hoff (2025).
#' Unless transformed, the function will return the values of \eqn{\lambda_1 (\hat{C})}.
#'
#' @param n the sample size.
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param trans logical,whether to transform the test statistic to approximate \code{TW}-law or not; TRUE by default.
#' @param sigma.known logical, whether to assume the known population covariance matrix or not; FALSE by default. If TRUE, the transformed test statistic is returned, regardless of \code{trans}.
#' @param samp.num the number of iterations for the Monte-Carlo simulation; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return a vector of \code{samp.num} Monte-Carlo simulated test statistics based on the largest eigenvalue of the sample core.
#'
#' @author Bongjung Sung
#'
#' @examples
#'
#' p1=40; p2=40; n=6400
#' y=p1*p2/n
#' set.seed(100)
#' test.stat=large.eig.null(n,p1,p2,center=FALSE)
#' test.stat
#' hist(test.stat,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Null Distribution")
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @export
large.eig.null=function(n,p1,p2,center=TRUE,trans=TRUE,sigma.known=FALSE,samp.num=1000,iter=10){

  p=p1*p2

  test.stat=NULL

  for(i in 1:samp.num){

    X=matrix(rnorm(n*p),ncol=p)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.kcd=covKCD::covKCD(S,p1,p2)

    if(n>=p){
      null.stat=RSpectra::eigs_sym(S.kcd$C,1)$values
    }else{
      K1.inv.root=sym.inv.root(S.kcd$K1)
      K2.inv.root=sym.inv.root(S.kcd$K2)
      X.whit=tcrossprod(X,kronecker(K2.inv.root,K1.inv.root))
      null.stat=(RSpectra::svds(X.whit,1)$d)^2/n
    }

    if(sigma.known==FALSE){
      if(trans==TRUE){
        y=p/n
        E.plus=(1+sqrt(y))^2
        gamma0=(sqrt(y)/(1+sqrt(y))^4)^(1/3)
        null.stat=gamma0*n^(2/3)*(null.stat-E.plus)
      }
      test.stat=c(test.stat,null.stat)
    }else{
      para=tw.para(S.kcd$K1,S.kcd$K2,n)
      null.stat=para[1]*n^(2/3)*(null.stat-para[2])
      test.stat=c(test.stat,null.stat)
    }

    if(i%%iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  return(test.stat)

}


#' Test statistic based on the largest eigenvalue of the sample core under the alternative regime
#'
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @import RSpectra
#' @import covKCD
#' @export
large.eig.alt=function(n,p1,p2,sigma,center=TRUE,trans=TRUE,sigma.known=FALSE,samp.num=1000,iter=10){



}
#' Empirical power of the test based on the largest eigenvalue of the sample core.
#'
#'
#'
#' @param dat the \eqn{n \times p1 \times p2} data tensor.
#' @param alpha the level of test; 0.05 by default.
#' @param center logical,
#' @param trans logical,
#' @param samp.num
#' @param iter
#'
#' @return \code{large.eig.power} returns a list with the following elements:
#' \describe{
#' \item{test.stat}{the row covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{null.stat}{the column covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{para.pval}{the factor matrix of a rank-1 core of dimension \eqn{p1 \times p1} ;}
#' \item{nonpara.pval}{If \code{lambda.gen}=TRUE, the non-spiked eigenvalue \eqn{\lambda \in (0,1)}. Otherwise, \code{NULL}.}
#' }
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=20; n=1600
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @import RMTstat
#' @export
large.eig.power=function(dat,alpha=0.05,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  # dimension of data
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]; p=p1*p2

  # computing the test statistic

  if(n>=p){
    S=dat2cov(dat,center=center)
    S.core=covKCD::covKCD(S,p1,p2)$C
    test.stat=RSpectra::eigs_sym(S.core,1)$values
  }else{
    dat=aperm(dat,perm=c(2,3,1))
    X=matrix(dat,ncol=n)
    X=scale(X,center=center,scale=FALSE)
    K1.inv.root=sym.inv.root(S.kcd$K1)
    K2.inv.root=sym.inv.root(S.kcd$K2)
    X.whit=crossprod(kronecker(K2.inv.root,K1.inv.root),X)
    test.stat=(RSpectra::svds(X.whit,1)$d)^2/n
  }

  # Monte-Carlo simulations of null test statistics
  print("Simulating the null test statistics.")
  null.stat=large.eig.null(n,p1,p2,center=FALSE,trans=trans,samp.num=samp.num,iter=iter)

  # Monte-Carlo based empirical power
  para.pval=mean(test.stat<null.stat)
  # Nonparametric power by comparing the test statistic against the Tracy-Widom law.
  npara.pval=1-RMTstat::ptw(test.stat,beta=1)

  return(list(test.stat=test.stat,null.stat=null.stat,para.pval=para.pval,npara.pval=npara.pval))

}
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @export
elrt.null=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

}
#' Simulating the extended LRT statistics under the alternative regime.
#'
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @export
elrt.alt=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

}
#' Empirical power based on the exended LRT.
#
#
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @import RSpectra
#' @import covKCD
#' @export
elrt.power=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

}
#' Simulating the separable expansion test statistics under the null hypothesis
#'
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @export
sep.exp.null=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  p=p1*p2

  test.stat=NULL
  set.seed(seed)

  for(i in 1:samp.iter){

    X=matrix(rnorm(n*p),ncol=p)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.core=covKCD::covKCD(S,p1,p2)$C
    test.stat=c(test.stat,sum(S.core^2)/p-1)

    if(i%%samp.iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  if(trans==TRUE){
    test.stat=n*test.stat-p-1
  }

  return(test.stat)

}
#' Simulating the separable expansion test statistics under the alternative regime
#'
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#'
#' @import RSpectra
#' @import covKCD
#' @export
sep.exp.alt=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  p=p1*p2

  test.stat=NULL
  set.seed(seed)

  for(i in 1:samp.iter){

    X=matrix(rnorm(n*p),ncol=p)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.core=covKCD::covKCD(S,p1,p2)$C
    test.stat=c(test.stat,sum(S.core^2)/p-1)

    if(i%%samp.iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  if(trans==TRUE){
    test.stat=n*test.stat-p-1
  }

  return(test.stat)

}


#' Empirical power of the test based on the separable expansion
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=50
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#'
#' @import RSpectra
#' @import covKCD
#' @export
sep.exp.power=function(dat,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  p=p1*p2

  test.stat=NULL
  set.seed(seed)

  for(i in 1:samp.iter){

    X=matrix(rnorm(n*p),ncol=p)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.core=covKCD::covKCD(S,p1,p2)$C
    test.stat=c(test.stat,sum(S.core^2)/p-1)

    if(i%%samp.iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }
  return(test.stat)

}

