#' Partial-isotropy rank-1 core
#'
#' Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-1 core for square matrix-variate data.
#'
#' If a core covariance matrix is of rank-1, the data should be a square matrix. Thus, the row and column dimensions must coincide.
#' The covariance matrix \eqn{\Sigma} with a partial-isotropy rank-1 core takes the form of
#' \deqn{\Sigma = (K_2 \otimes K_1)^{1/2}((1- \lambda) \mathrm{vec}(A) \mathrm{vec}(A)^\top + \lambda I_{p_1^2})(K_2 \otimes K_1)^{1/2, \top}}
#' for \eqn{\lambda \in (0,1)}, positive definite matrices \eqn{K_1, K_2} of the same dimensions \eqn{p_1 \times p_1}, and \eqn{A \in \mathbb{R}^{p_1 \times p_1}} such that \eqn{\mathrm{vec}(A) \mathrm{vec}(A)^\top} is a rank-1 core.
#' Here \eqn{M^{1/2}} denotes either symmetric square root or the Cholesky factor of a positive definite matrix \eqn{M}.
#' The separable and core components of \eqn{\Sigma}, denoted \eqn{K} and \eqn{C}, are
#' \deqn{K=K_2\otimes K_1,\quad C=(1- \lambda) \mathrm{vec}(A) \mathrm{vec}(A)^\top + \lambda I_{p_1^2}.}
#' For the exact formula, see Theorem 1 of Sung and Hoff (2025).
#'
#' @param p1 the row and column dimensions.
#' @param lambda.gen logical, whether to generate a non-spiked eigenvalue or not.; TRUE by default.
#'
#' @return \code{pi.rank1.core} returns a list with the following elements:
#' \describe{
#' \item{K1}{the row covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{K2}{the column covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{A}{the factor matrix of a rank-1 core of dimension \eqn{p1 \times p1} ;}
#' \item{lambda}{If \code{lambda.gen}=TRUE, the non-spiked eigenvalue \eqn{\lambda \in (0,1)}. Otherwise, \code{NULL}.}
#' }
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' p1=10
#' pi.rank1.core(p1)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import pracma
#' @import stats
#' @export
pi.rank1.core=function(p1,lambda.gen=TRUE){

  p=p1^2; n=2*p1

  # randomly generate the non-spiked eigenvalue in (0,1) if lambda.gen==TRUE
  if(lambda.gen==TRUE){
    lambda=stats::runif(1,0,1)
  }else{
    lambda=NULL
  }

  # randomly generate the row and column covariance matrices K1 and K2.
  K1=matrix(stats::rnorm(p1*n),ncol=p1)
  K1=crossprod(K1,K1)/n

  K2=matrix(stats::rnorm(p1*n),ncol=p1)
  K2=crossprod(K2,K2)/n

  # randomly generate the factor matrix for rank-1 core.
  A=pracma::randortho(p1)*sqrt(p1)

  return(list(K1=K1,K2=K2,A=A,lambda=lambda))

}


#' Partial-isotropy rank-2 core
#'
#' Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-2 core for matrix-variate data.
#'
#' If a core covariance matrix is of rank-2, the dimension \eqn{(p_1, p_2)} should satisfy one of the followings: \eqn{p_1 = p_2} or \eqn{p_1 - p_2 | p_1} \eqn{(p_1 > p_2)}, or \eqn{p_2 - p_1 | p_2} \eqn{(p_2 > p_1)}.
#' The covariance matrix \eqn{\Sigma} with a partial-isotropy rank-2 core takes the form of
#' \deqn{\Sigma = (K_2 \otimes K_1)^{1/2}((1- \lambda) A A^\top + \lambda I_{p_1 p_2})(K_2 \otimes K_1)^{1/2, \top}}
#' for \eqn{\lambda \in (0,1)}, positive definite matrices \eqn{K_1} and \eqn{K_2} of the dimensions \eqn{p_1 \times p_1} and \eqn{p_2 \times p_2}, respectively, and \eqn{A \in \mathbb{R}^{p_1p_2\times 2}} whose \eqn{i}th column is a vectorization of \eqn{p_1 \times p_2} matrix \eqn{A_i} such that \eqn{AA^\top} is a rank-2 core.
#' Here \eqn{M^{1/2}} denotes either symmetric square root or the Cholesky factor of a positive definite matrix \eqn{M}.
#' The separable and core components of \eqn{\Sigma}, denoted \eqn{K} and \eqn{C}, are
#' \deqn{K=K_2\otimes K_1,\quad C=(1- \lambda) A A^\top + \lambda I_{p_1 p_2}.}
#' For the exact formula, see Theorem 1 of Sung and Hoff (2025).
#'
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param lambda.gen logical, whether to generate a non-spiked eigenvalue or not; TRUE by default.
#'
#' @return \code{pi.rank2.core} returns a list with the following elements:
#' \describe{
#' \item{K1}{the row covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{K2}{the column covariance matrix of dimension \eqn{p2 \times p2};}
#' \item{A}{an array of factor matrices of a rank-2 core of dimension \eqn{p1 \times p2 \times 2};}
#' \item{lambda}{If \code{lambda.gen}=TRUE, the non-spiked eigenvalue \eqn{\lambda \in (0,1)}. Otherwise, \code{NULL}.}
#' }
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' # It must be that p1=p2 or p1-p2|p1 (p1>p2) or p2-p1|p2 (p2>p1).
#' p1=12; p2=10
#' pi.rank2.core(p1,p2)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import pracma
#' @import stats
#' @export
pi.rank2.core=function(p1,p2,lambda.gen=TRUE){

  p.max=max(p1,p2); p.min=min(p1,p2)
  if(p.max!=p.min && (p.max%%(p.max-p.min))!=0){
    stop("Incorrect choice of (p1,p2): if p1>=p2, p1=p2 or p1-p2|p1 with the converse when p1<p2.")
  }

  # randomly generate the non-spiked eigenvalue in (0,1) if lambda.gen==TRUE
  if(lambda.gen==TRUE){
    lambda=stats::runif(1,0,1)
  }else{
    lambda=NULL
  }

  p=p1*p2

  n1=2*p1; n2=2*p2

  # randomly generate the row and column covariance matrices K1 and K2
  K1=matrix(stats::rnorm(n1*p1),ncol=p1)
  K1=crossprod(K1,K1)/n1
  K2=matrix(stats::rnorm(n2*p2),ncol=p2)
  K2=crossprod(K2,K2)/n2

  # randomly generate the factor matrices for rank-2 core
  if(p.max==p.min){

    U=pracma::randortho(p1); V=pracma::randortho(p1)
    D1=diag(sort(stats::runif(p1,0,1),decreasing=TRUE))
    D2=diag(sqrt(1-diag(D1)^2)*sample(c(-1,1),size=p1,replace=TRUE))
    A=array(0,dim=c(p1,p1,2))
    A[,,1]=sqrt(p1)*tcrossprod(tcrossprod(U,D1),V)
    A[,,2]=sqrt(p1)*tcrossprod(tcrossprod(U,D2),V)

  }else{

    m=pracma::gcd(p.max,p.min)
    k=p.min/m

    A=array(0,dim=c(p.max,p.min,2))

    U=pracma::randortho(p.max)
    D1=diag(rep(sqrt(k:1),each=m))
    D2=diag(rep(sqrt(1:k),each=m))
    V=pracma::randortho(p.min)
    A1=tcrossprod(tcrossprod(U[,1:p.min],D1),V)
    A2=tcrossprod(tcrossprod(U[,(m+1):p.max],D2),V)
    A1=sqrt(m)*A1; A2=sqrt(m)*A2
    A[,,1]=A1; A[,,2]=A2

    if(p1<p2){
      A=aperm(A,perm=c(2,1,3))
    }

  }

  return(list(K1=K1,K2=K2,A=A,lambda=lambda))

}


#' Partial-isotropy rank-r core.
#'
#' Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-\eqn{r} core for matrix-variate data.
#'
#' If a core covariance matrix is of rank-\eqn{r} for general \eqn{r}, the dimension \eqn{(p_1, p_2)} should satisfy one of the followings: \eqn{p_1/p_2 + p_2/p_1 <r} or \eqn{p_1 = p_2 r} \eqn{(p_1 \geq p_2)}, or \eqn{p_2 = p_1 r} \eqn{(p_2 \geq p_1)}.
#' The covariance matrix \eqn{\Sigma} with a partial-isotropy rank-\eqn{r} core takes the form of
#' \deqn{\Sigma = (K_2 \otimes K_1)^{1/2}((1- \lambda) A A^\top + \lambda I_{p_1 p_2})(K_2 \otimes K_1)^{1/2, \top}}
#' for \eqn{\lambda \in (0,1)}, positive definite matrices \eqn{K_1} and \eqn{K_2} of the dimensions \eqn{p_1 \times p_1} and \eqn{p_2 \times p_2}, respectively, and \eqn{A \in \mathbb{R}^{p_1p_2\times r}} whose \eqn{i}th column is a vectorization of \eqn{p_1 \times p_2} matrix \eqn{A_i} such that \eqn{AA^\top} is a rank-\eqn{r} core.
#' Here \eqn{M^{1/2}} denotes either symmetric square root or the Cholesky factor of a positive definite matrix \eqn{M}.
#' The separable and core components of \eqn{\Sigma}, denoted \eqn{K} and \eqn{C}, are
#' \deqn{K=K_2\otimes K_1,\quad C=(1- \lambda) A A^\top + \lambda I_{p_1p_2}.}
#' For the exact formula when \eqn{p_1=p_2 r} or \eqn{p_2=p_1 r}, see Theorem 1 of Sung and Hoff (2025).
#'
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param r the partial-isotropy rank.
#' @param lambda.gen logical, whether to generate a non-spiked eigenvalue or not; TRUE by default.
#'
#' @return \code{pi.rank_r.core} returns a list with the following elements:
#' \describe{
#' \item{K1}{the row covariance matrix of dimension \eqn{p1 \times p1};}
#' \item{K2}{the column covariance matrix of dimension \eqn{p2 \times p2};}
#' \item{A}{an array of factor matrices of a rank-r core of dimension \eqn{p1 \times p2 \times r};}
#' \item{lambda}{If \code{lambda.gen}=TRUE, the non-spiked eigenvalue \eqn{\lambda \in (0,1)}. Otherwise, \code{NULL}.}
#' }
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' # It must be that p1/p2+p2/p1<r or p1=p2r or p2=p1r
#' p1=12; p2=8; r=3
#' pi.rank_r.core(p1,p2,r)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import pracma
#' @import covKCD
#' @import stats
#' @export
pi.rank_r.core=function(p1,p2,r,lambda.gen=TRUE){

  p.max=max(p1,p2); p.min=min(p1,p2)
  thres=p.max/p.min+p.min/p.max
  if(r<=thres && (p.max/p.min!=r)){
    stop("Incorrect choice of (p1,p2,r): p1/p2+p2/p1<r or p1=p2r if p1>=p2 with the converse when p1<p2.")
  }

  # randomly generate the non-spiked eigenvalue in (0,1) if lambda.gen==TRUE
  if(lambda.gen==TRUE){
    lambda=stats::runif(1,0,1)
  }else{
    lambda=NULL
  }

  p=p1*p2

  n1=2*p1; n2=2*p2

  # randomly generate the row and column covariance matrices K1 and K2
  K1=matrix(stats::rnorm(n1*p1),ncol=p1)
  K1=crossprod(K1,K1)/n1
  K2=matrix(stats::rnorm(n2*p2),ncol=p2)
  K2=crossprod(K2,K2)/n2

  # randomly generate the factor matrices for rank-r core
  if(thres<r){

    A=matrix(stats::rnorm(p*r),ncol=r)
    A.kcd=covKCD::covKCD(tcrossprod(A,A),p1,p2)
    K1.inv.root=sym.inv.root(A.kcd$K1)
    K2.inv.root=sym.inv.root(A.kcd$K2)
    A=crossprod(kronecker(K2.inv.root,K1.inv.root),A)
    A=array(A,dim=c(p1,p2,r))

  }else{

    A=pracma::randortho(p.max)*sqrt(p.min)
    if(p1>=p2){
      A=array(A,dim=c(p1,p2,r))
    }else{
      A=array(A,dim=c(p2,p1,r))
      A=aperm(A,perm=c(2,1,3))
    }
    A=matrix(A,ncol=r)

  }

  return(list(K1=K1,K2=K2,A=A,lambda=lambda))

}


#' A covariance matrix with a partial-isotropy rank-r core.
#'
#' Given a list of parameters constituting the covariance matrix with a partial-isotropy rank-\eqn{r} core, the function assembles the covariance matrix using these parameters.
#'
#' @param para.list the list of parameters constituting the covariance matrix with a partial-isotropy rank-\eqn{r} core.
#' @param lambda0 the pre-specified value of the non-spiked eigenvalue \code{lambda} in a partial-isotropy core. If it is not specified in \code{para.list} or does not belong to (0,1), it will be randomly generated; 0.5 by deafult.
#' @param root the choice of the square root of positive definite matrices; must be either \code{"sym"} (symmetric square root) or \code{"chol"} (Cholesky factor), \code{"sym"} by default.
#'
#' @return a \eqn{p1p2 \times p1p2} covariance matrix \eqn{\Sigma} with a partial-isotropy rank-\eqn{r} core.
#' The attribute \eqn{\lambda} of \eqn{\Sigma} denotes the value of the non-spiked eigenvalue.
#'
#' @author Bongjung Sung
#' @examples
#'
#' set.seed(100)
#' # generate a list of parameters for a covariance matrix of a partial-isotropy core.
#' p1=14; p2=10; r=3
#' para.list=pi.rank_r.core(p1,p2,r)
#' # assembles the covariance matrix using the above list of parameters.
#' pi.core(para.list)
#'
#' @import stats
#' @export
pi.core=function(para.list,lambda0=0.5,root="sym"){

  if(root!="sym" && root!="chol"){
    stop("root should be sym (symmetric square root) or chol (Cholesky factor).")
  }

  K1=para.list$K1
  K2=para.list$K2
  A=para.list$A
  lambda=para.list$lambda

  # If a non-spiked eigenvalue lambda is not specified in para.list, then it is either
  # randomly generated or set to be the value as lambda0 provided that lambda0 is in (0,1).
  # This value can be readily verified with the attribute of the resulting object (Sigma).
  if(is.null(lambda)){
    if(0<lambda0 && lambda0<1){
      lambda=lambda0
    }else{
      print("Invalid choice of lambda0: it should be in (0,1). A non-spiked eigenvalue is thus randomly generated from Unif(0,1).")
      lambda=stats::runif(1,0,1)
    }
  }

  p1=ncol(K1); p2=ncol(K2); p=p1*p2

  if(root=="sym"){
    K1.root=sym.root(K1); K2.root=sym.root(K2)
  }else{
    K1.root=t(chol(K1)); K2.root=t(chol(K2))
  }

  K.root=kronecker(K2.root,K1.root)
  r=length(A)/p
  A=matrix(A,nrow=p,ncol=r)
  core=K.root%*%A
  Sigma=(1-lambda)*tcrossprod(core,core)+lambda*kronecker(K2,K1)
  attr(Sigma,"lambda")=lambda
  return(Sigma)

}
