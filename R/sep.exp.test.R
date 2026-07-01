#' Test statistic based on the separable expansion test under the null
#'
#' Implement the Monte-Carlo simulations of the test statistic based on the separable expansion of the sample core \eqn{\hat{C}} for Gaussian populations as benchmark cases under the null hypothesis of separability.
#' Here the sample core \eqn{\hat{C}} refers to the core component of the sample covariance matrix \eqn{S}. The test statistic is given by
#' \deqn{ T(Y)=|| \hat{C} ||_F^2/p-1,}
#'  where \eqn{p = p_1 p_2} for row and column dimensions, \eqn{p_1} and \eqn{p_2}, respectively.
#' The test statistic may be transformed as \eqn{nT(Y)-p-1} if \code{trans}=TRUE.
#' For the details, see Section 3.3 and 5.1 of Sung and Hoff (2025). Also, if the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, the data center should be centered.
#'
#' @param n the sample size.
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param trans logical,whether to transform the test statistic; TRUE by default.
#' @param samp.num the number of iterations for the Monte-Carlo simulation; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return the vector of \code{samp.num} Monte-Carlo simulated null test statistics based on the separable expansion of the sample core.
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=5; p2=3; n=60
#'
#' set.seed(100)
#' sep.exp.null(n,p1,p2,center=FALSE,samp.num=20)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import covKCD
#' @import stats
#' @export
sep.exp.null=function(n,p1,p2,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  p=p1*p2

  test.stat=NULL

  for(i in 1:samp.num){

    # generate the data under the null
    # assuming Sigma=I_p by virtue of the Kronecker-invariance.
    X=matrix(stats::rnorm(n*p),ncol=p)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.core=covKCD::covKCD(S,p1,p2)$C

    null.stat=sum(S.core^2)/p-1

    if(trans==TRUE){
      null.stat=null.stat*n-p-1
    }
    test.stat=c(test.stat,null.stat)

    if(i%%iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  return(test.stat)

}


#' Test statistic based on the separable expansion test under the alternative regime.
#'
#' Implement the Monte-Carlo simulations of the test statistic based on the separable expansion of the sample core \eqn{\hat{C}} for Gaussian populations as benchmark cases under the alternative regime.
#' Here the sample core \eqn{\hat{C}} refers to the core component of the sample covariance matrix \eqn{S}. The test statistic is given by
#' \deqn{ T(Y)=|| \hat{C} ||_F^2/p-1,}
#'  where \eqn{p = p_1 p_2} for row and column dimensions, \eqn{p_1} and \eqn{p_2}, respectively.
#' The test statistic may be transformed as \eqn{nT(Y)-p-1} if \code{trans}=TRUE. For the details, see Section 3.3 and 5.1 of Sung and Hoff (2025).
#' The population covariance matrix should be specified in \code{sigma}. Also, if the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, the data center should be centered.
#'
#' @param n the sample size.
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param sigma the population covariance matrix.
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param trans logical,whether to transform the test statistic; TRUE by default.
#' @param samp.num the number of iterations for the Monte-Carlo simulation; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return the vector of \code{samp.num} Monte-Carlo simulated test statistics based on the separable expansion of the sample core under the alternative.
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=5; p2=3; r=4; n=100
#'
#' set.seed(100)
#' para.list=pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
#' Sigma=pi.core(para.list,lambda0=0.95)
#'
#' sep.exp.alt(n,p1,p2,Sigma,center=FALSE,samp.num=20)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import covKCD
#' @import stats
#' @export
sep.exp.alt=function(n,p1,p2,sigma,center=TRUE,trans=TRUE,samp.num=1000,iter=10){

  p=p1*p2

  if(ncol(sigma)!=p){
    stop("Incorrectly specified covariance matrix: the dimension should match with the data.")
  }
  if((matrixcalc::is.positive.definite(sigma))!=TRUE){
    stop("Incorrectly specified covariance matrix: it should be strictly positive definite.")
  }

  # symmetric root of covariance matrix
  sigma.root=sym.root(sigma)

  test.stat=NULL

  for(i in 1:samp.num){

    X=tcrossprod(matrix(stats::rnorm(n*p),ncol=p),sigma.root)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.core=covKCD::covKCD(S,p1,p2)$C

    alt.stat=sum(S.core^2)/p-1

    if(trans==TRUE){
      alt.stat=alt.stat*n-p-1
    }
    test.stat=c(test.stat,alt.stat)

    if(i%%iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  return(test.stat)

}


#' Empirical power of the test based on the separable expansion the sample core under Gaussian populations
#'
#' Evaluate the empirical power of the test based on the separable expansion of the sample core \eqn{\hat{ C }} with \eqn{n} i.i.d. random matrices generated according to \eqn{N_{p_1 \times p_2} (0, \Sigma)} with given \eqn{\Sigma} (specified in \code{sigma}).
#' The test statistic is given by \eqn{T(Y) = || \hat{C} ||_F^2/p-1}. The power is evaluated with respect to \eqn{nT(Y)-p-1}.
#' Given a level \eqn{\alpha \in (0,1)} (specified in \code{alpha}), the parametric power is evaluated. To this end, the p-value (\code{para.pval}) is first evaluated for each test statistic by comparing it to the Monte-Carlo approximated null distribution (\code{null.stat}).
#' Then the power (\code{para.power}) is evaluated as the proportion of these p-values smaller than \eqn{\alpha}. For details, see Section 5.2 of Sung and Hoff (2025).
#' If the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, you should center the data (\code{center=TRUE}).
#'
#' @param n the sample size.
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param sigma the population covariance matrix.
#' @param alpha the level of the test.
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param null.samp.num the number of iterations for the Monte-Carlo simulation of the test statistics under the null; 1000 by default.
#' @param alt.samp.num the number of iterations for the Monte-Carlo simulation of the test statistics under the alternative; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return \code{sep.exp.power} returns a list of the following elements:
#' \describe{
#' \item{alt.stat}{the vector of \code{alt.samp.num} Monte-Carlo simulated test statistics under \eqn{N_{p_1 \times p_2}(0, \code{sigma})} after some transformation;}
#' \item{null.stat}{the vector of \code{null.samp.num} Monte-Carlo simulated test statistics under the null after some transformation;}
#' \item{para.pval}{the vector of p-values for each test statistic in \code{alt.stat} evaluated based on Monte-Carlo approximated empirical null distribution (\code{null.stat});}
#' \item{para.power}{the proportion of the p-values in \code{para.pval} smaller than \code{alpha};}
#' }
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=12; p2=10; r=4; n=100
#'
#' set.seed(100)
#' para.list=pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
#' Sigma=pi.core(para.list,lambda0=0.98)
#'
#' sep.exp.power(n,p1,p2,Sigma,center=FALSE,null.samp.num=20,alt.samp.num=20)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @export
sep.exp.power=function(n,p1,p2,sigma,alpha=0.05,center=TRUE,null.samp.num=1000,alt.samp.num=1000,iter=10){

  print(paste("Simulating the test statistics under the alternative."))
  alt.stat=sep.exp.alt(n,p1,p2,sigma,center=center,samp.num=alt.samp.num,iter=iter)

  print(paste("Simulating the test statistics under the null."))
  null.stat=sep.exp.null(n,p1,p2,center=FALSE,samp.num=null.samp.num,iter=iter)

  # computing p-values
  para.pval=apply(matrix(alt.stat,ncol=alt.samp.num),2,FUN=function(x){mean(x<null.stat)})

  # computing power
  para.power=mean(para.pval<alpha)

  return(list(alt.stat=alt.stat,null.stat=null.stat,para.pval=para.pval,para.power=para.power))

}


#' Empirical power of the test based on the separable expansion of the sample core with the given data
#'
#' Evaluate the empirical power of the test based on the separable expansion of the sample core \eqn{\hat{ C }} with the tensor data, \eqn{Y \in \mathbb{R}^{n \times p_1 \times p_2}}.
#' The test statistic is given by \eqn{T(Y) = || \hat{C} ||_F^2/p-1}. The power is evaluated with respect to \eqn{nT(Y)-p-1}.
#' Given a level \eqn{\alpha \in (0,1)} (specified in \code{alpha}), the parametric power is evaluated. To this end, the p-value (\code{para.pval}) is first evaluated for each test statistic by comparing it to the Monte-Carlo approximated null distribution (\code{null.stat}).
#' Then the power (\code{para.power}) is evaluated as the proportion of these p-values smaller than \eqn{\alpha}. For details, see Section 5.2 of Sung and Hoff (2025).
#' If the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, you should center the data (\code{center=TRUE}).
#'
#' @param dat the \eqn{n \times p1 \times p2} data tensor (sample x row x column).
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param samp.num  the number of iterations for simulating the transformed null test statistic; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return \code{sep.exp.power.dat} returns a list with the following elements:
#' \describe{
#' \item{test.stat}{the computed test statistic based on \code{dat} after some transformation;}
#' \item{null.stat}{the vector of \code{samp.num} Monte-Carlo simulated null test statistics after some transformation;}
#' \item{para.pval}{the p-value evaluated based on Monte-Carlo approximated empirical null distribution (\code{null.stat});}
#' }
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=10; p2=12; r=4; n=100
#' p=p1*p2
#'
#' set.seed(100)
#' para.list=pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
#' Sigma=pi.core(para.list,lambda0=0.99)
#' Sigma.root=sym.root(Sigma)
#'
#' dat=crossprod(Sigma.root,matrix(rnorm(n*p),ncol=n))
#' dat=array(dat,dim=c(p1,p2,n))
#' dat=aperm(dat,perm=c(3,1,2))
#' sep.exp.power.dat(dat,center=FALSE,samp.num=20)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import covKCD
#' @export
sep.exp.power.dat=function(dat,center=TRUE,samp.num=1000,iter=10){

  # dimension of data
  dat.dim=dim(dat)
  n=dat.dim[1]; p1=dat.dim[2]; p2=dat.dim[3]; p=p1*p2

  # computing the test statistic

  S=dat2cov(dat,center=center)
  S.core=covKCD::covKCD(S,p1,p2)$C

  test.stat=n/p*sum(S.core^2)-p-1

  # Monte-Carlo simulations of null test statistics
  print("Simulating the null test statistics.")
  null.stat=sep.exp.null(n,p1,p2,center=FALSE,samp.num=samp.num,iter=iter)

  # Monte-Carlo based empirical power
  para.pval=mean(test.stat<null.stat)

  return(list(test.stat=test.stat,null.stat=null.stat,para.pval=para.pval))

}

