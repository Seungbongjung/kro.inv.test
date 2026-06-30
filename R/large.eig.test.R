#' Test statistic based on the largest eigenvalue of the sample core under the null
#'
#' Implement the Monte-Carlo simulations of the test statistic based on the largest eigenvalue of the sample core \eqn{\hat{C}} for Gaussian populations as benchmark cases under the null hypothesis of separability.
#' Here the sample core \eqn{\hat{C}} refers to the core component of the sample covariance matrix \eqn{S}.
#' The test statistic may be transformed to approximate the Tracy-Widom (\code{TW}) law if \code{trans}=TRUE. This transformation is based on the quantities depending only on \eqn{(n,p_1,p_2)}.
#' If the population covariance matrix under the null is assumed to be known (\code{sigma.known}=TRUE), the transformation is done regardless of \code{trans} based on the Kronecker MLE (the separable component of \eqn{S}).
#' For the details, see (20) of Sung and Hoff (2025). Unless transformed, the function will return the values of \eqn{\lambda_1 (\hat{C})}.
#' Also, if the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, the data center should be centered.
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
#' @return the vector of \code{samp.num} Monte-Carlo simulated test statistics based on the largest eigenvalue of the sample core.
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=40; p2=40; n=6400
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

    # generate the data under the null
    # assuming Sigma=I_p by virtue of the Kronecker-invariance.
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

      # transform the test statistic using the Kronecker MLE

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
#' Implement the Monte-Carlo simulations of the test statistic based on the largest eigenvalue of the sample core \eqn{\hat{C}} for Gaussian populations as benchmark cases under the alternative regime.
#' Here the sample core \eqn{\hat{C}} refers to the core component of the sample covariance matrix \eqn{S}.
#' The test statistic may be transformed to approximate the Tracy-Widom (\code{TW}) law if \code{trans}=TRUE, particularly under the local alternative regime (see Theorem 2 of Sung and Hoff (2025)). This transformation is based on the quantities depending only on \eqn{(n,p_1,p_2)}.
#' If the population covariance matrix of the mean is assumed to be known (\code{sigma.known}=TRUE), the transformation is done regardless of \code{trans} based on the Kronecker MLE (the separable component of \eqn{S}).
#' For the details, see (20) of Sung and Hoff (2025). Unless transformed, the function will return the values of \eqn{\lambda_1 (\hat{C})}.
#' Also, if the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, the data center should be centered.
#'
#' @param n the sample size.
#' @param p1 the row dimension.
#' @param p2 the column dimension.
#' @param sigma the population covariance matrix.
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param trans logical,whether to transform the test statistic to approximate \code{TW}-law or not; TRUE by default.
#' @param sigma.known logical, whether to assume the known population covariance matrix or not; FALSE by default. If TRUE, the transformed test statistic is returned, regardless of \code{trans}.
#' @param samp.num the number of iterations for the Monte-Carlo simulation; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return the vector of \code{samp.num} Monte-Carlo simulated test statistics based on the largest eigenvalue of the sample core under \eqn{}.
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=640
#'
#' set.seed(100)
#' para.list=pi.rank2.core(p1,p2,lambda.gen=FALSE)
#'
#' # local alternative
#' lambda1=1-(0.5/(320/2+0.5))
#' sigma1=pi.core(para.list,lambda0=lambda1)
#' test.stat1=large.eig.alt(n,p1,p2,sigma=sigma1,center=FALSE)
#' hist(test.stat1,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Distribution")
#'
#' # dense alternative
#' lambda2=1-(1.2/(320/2+1.2))
#' sigma2=pi.core(para.list,lambda0=lambda2)
#' test.stat2=large.eig.alt(n,p1,p2,sigma=sigma2,center=FALSE)
#' hist(test.stat2,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Distribution")
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @import matrixcalc
#' @export
large.eig.alt=function(n,p1,p2,sigma,center=TRUE,trans=TRUE,sigma.known=FALSE,samp.num=1000,iter=10){

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

    X=tcrossprod(matrix(rnorm(n*p),ncol=p),sigma.root)
    if(center==TRUE){
      X=scale(X,center=center,scale=FALSE)
    }

    S=crossprod(X,X)/n
    S.kcd=covKCD::covKCD(S,p1,p2)

    if(n>=p){
      alt.stat=RSpectra::eigs_sym(S.kcd$C,1)$values
    }else{
      K1.inv.root=sym.inv.root(S.kcd$K1)
      K2.inv.root=sym.inv.root(S.kcd$K2)
      X.whit=tcrossprod(X,kronecker(K2.inv.root,K1.inv.root))
      alt.stat=(RSpectra::svds(X.whit,1)$d)^2/n
    }

    if(sigma.known==FALSE){
      if(trans==TRUE){
        y=p/n
        E.plus=(1+sqrt(y))^2
        gamma0=(sqrt(y)/(1+sqrt(y))^4)^(1/3)
        alt.stat=gamma0*n^(2/3)*(alt.stat-E.plus)
      }
      test.stat=c(test.stat,alt.stat)
    }else{

      # transform the test statistic using the Kronecker MLE

      para=tw.para(S.kcd$K1,S.kcd$K2,n)
      null.stat=para[1]*n^(2/3)*(alt.stat-para[2])
      test.stat=c(test.stat,alt.stat)
    }

    if(i%%iter==0){
      print(paste(i,"th iteration done.",sep=""))
    }

  }

  return(test.stat)

}


#' Empirical power of the test based on the largest eigenvalue of the sample core under Gaussian populations
#'
#' Evaluate the empirical power of the test based on \eqn{\lambda_1 (\hat{C})} with \eqn{n} i.i.d. random matrices generated according to \eqn{N_{p_1 \times p_2} (0, \Sigma)} for the sample core \eqn{\hat{C}} with given \eqn{\Sigma} (specified in \code{sigma}).
#' The power is evaluated with the transformed \eqn{\lambda_1 (\hat{C})} so that it may follow the Tracy-Widom law under the null and some local alternative regimes (see Theorem 2 of Sung and Hoff (2025)).
#' Given a level \eqn{\alpha \in (0,1)} (specified in \code{alpha}), both parametric and nonparametric power are evaluated. For the parametric power, the p-value (\code{para.pval}) is first evaluated for each test statistic by comparing it to the Monte-Carlo approximated null distribution (\code{null.stat}).
#' Then the power (\code{para.power}) is evaluated as the proportion of these p-values smaller than \eqn{\alpha}. Similarly, to compute the nonparametric power, the p-value (\code{nonpara.pval}) is first evaluated for each test statistic by comparing it to the Tracy-Widom law.
#' Finally, the power (\code{nonpara.power}) is evaluated as an analogy to \code{para.power}. For details, see Section 5.2 of Sung and Hoff (2025).
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
#' @return \code{large.eig.power} returns a list of the following elements:
#' \describe{
#' \item{alt.stat}{the vector of \code{alt.samp.num} Monte-Carlo simulated test statistics under \eqn{N_{p_1 \times p_2}(0, \code{sigma}) after some transformation};}
#' \item{null.stat}{the vector of \code{null.samp.num} Monte-Carlo simulated test statistics under the null after some transformation;}
#' \item{para.pval}{the vector of p-values for each test statistic in \code{alt.stat} evaluated based on Monte-Carlo approximated empirical null distribution (\code{null.stat});}
#' \item{para.power}{the proportion of the p-values in \code{para.pval} smaller than \code{alpha};}
#' \item{nonpara.pval}{the vector of p-values for each test statistic in evaluated by comparing the test statistic to the Tracy-Widom law;}
#' \item{nonpara.power}{the proportion of the p-values in \code{nonpara.pval} smaller than \code{alpha}.}
#' }
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=20; p2=16; n=640
#'
#' set.seed(100)
#' para.list=pi.rank2.core(p1,p2,lambda.gen=FALSE)
#'
#' # local alternative
#' lambda1=1-(0.5/(320/2+0.5))
#' sigma1=pi.core(para.list,lambda0=lambda1)
#' large.eig.power(n,p1,p2,sigma1,center=FALSE)
#'
#' # dense alternative
#' lambda2=1-(0.8/(320/2+0.8))
#' sigma2=pi.core(para.list,lambda0=lambda2)
#' large.eig.power(n,p1,p2,sigma2,center=FALSE)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RMTstat
#' @export
large.eig.power=function(n,p1,p2,sigma,alpha=0.05,center=TRUE,null.samp.num=1000,alt.samp.num=1000,iter=10){

  print(paste("Simulating the test statistics under the alternative."))
  alt.stat=large.eig.alt(n,p1,p2,sigma,center=center,samp.num=alt.samp.num,iter=iter)

  print(paste("Simulating the test statistics under the null."))
  null.stat=large.eig.null(n,p1,p2,center=FALSE,samp.num=null.samp.num,iter=iter)

  # computing p-values
  para.pval=apply(matrix(alt.stat,ncol=alt.samp.num),2,FUN=function(x){mean(x<null.stat)})
  nonpara.pval=1-RMTstat::ptw(alt.stat,beta=1)

  # computing power
  para.power=mean(para.pval<alpha)
  nonpara.power=mean(nonpara.pval<alpha)

  return(list(alt.stat=alt.stat,null.stat=null.stat,para.pval=para.pval,nonpara.pval=nonpara.pval,
              para.power=para.power,nonpara.power=nonpara.power))

}

#' Empirical power of the test based on the largest eigenvalue of the sample core with the given data
#'
#' Evaluate the empirical power of the test based on \eqn{\lambda_1 (\hat{C})} with the tensor data, \eqn{Y \in \mathbb{R}^{n \times p_1 \times p_2}}.
#' The power is evaluated with the transformed \eqn{\lambda_1 (\hat{C})} so that it may follow the Tracy-Widom law under the null and some local alternative regimes (see Theorem 2 of Sung and Hoff (2025)).
#' The parametric p-value (\code{para.pval}) is evaluated for each test statistic by comparing it to the Monte-Carlo approximated null distribution (\code{null.stat}).
#' On the other hand, the nonparametric p-value (\code{nonpara.pval}) is evaluated for each test statistic by comparing it to the Tracy-Widom law.
#' For details, see Section 5.2 of Sung and Hoff (2025). If the mean is assumed to be known, you may not center the data (\code{center=FALSE}). Otherwise, you should center the data (\code{center=TRUE}).
#'
#' @param dat the \eqn{n \times p1 \times p2} data tensor (sample x row x column).
#' @param center logical, whether to center the data or not; TRUE by default.
#' @param samp.num  the number of iterations for simulating the transformed null test statistic; 1000 by default.
#' @param iter the unit number at which to print the number of current iterations; 10 by default.
#'
#' @return \code{large.eig.power.dat} returns a list with the following elements:
#' \describe{
#' \item{test.stat}{the computed test statistic based on \code{dat} after some transformation;}
#' \item{null.stat}{the vector of \code{samp.num} Monte-Carlo simulated null test statistics after some transformation;}
#' \item{para.pval}{the p-value evaluated based on Monte-Carlo approximated empirical null distribution (\code{null.stat});}
#' \item{nonpara.pval}{the p-value evaluated by comparing the test statistic to the Tracy-Widom law.}
#' }
#'
#' @author Bongjung Sung
#'
#' @examples
#' p1=30; p2=18; r=4; n=200
#' p=p1*p2
#'
#' set.seed(100)
#' para.list=pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
#'
#' # local alternative
#' Sigma1=pi.core(para.list,lambda0=0.999)
#' Sigma.root1=sym.root(Sigma1)
#' dat1=crossprod(Sigma.root1,matrix(rnorm(n*p),ncol=n))
#' dat1=array(dat1,dim=c(p1,p2,n))
#' dat1=aperm(dat1,perm=c(3,1,2))
#' large.eig.power.dat(dat1,center=FALSE)
#'
#' # dense alternative
#' Sigma2=pi.core(para.list,lambda0=0.8)
#' Sigma.root2=sym.root(Sigma2)
#' dat2=crossprod(Sigma.root2,matrix(rnorm(n*p),ncol=n))
#' dat2=array(dat2,dim=c(p1,p2,n))
#' dat2=aperm(dat2,perm=c(3,1,2))
#' large.eig.power.dat(dat2,center=FALSE)
#'
#' @references
#' Sung, B. and Hoff, P. (2025). Testing Separability of High-Dimensional Covariance Matrices.
#' arXiv preprint arXiv:2506.17463.
#'
#' @import RSpectra
#' @import covKCD
#' @import RMTstat
#' @export
large.eig.power.dat=function(dat,center=TRUE,samp.num=1000,iter=10){

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
    X=t(X)
    X=scale(X,center=center,scale=FALSE)
    S.kcd=covKCD::covKCD(crossprod(X,X)/n,p1,p2)
    K1.inv.root=sym.inv.root(S.kcd$K1)
    K2.inv.root=sym.inv.root(S.kcd$K2)
    X.whit=tcrossprod(kronecker(K2.inv.root,K1.inv.root),X)
    test.stat=(RSpectra::svds(X.whit,1)$d)^2/n
  }

  y=p/n
  E.plus=(1+sqrt(y))^2
  gamma0=(sqrt(y)/(1+sqrt(y))^4)^(1/3)
  test.stat=gamma0*n^(2/3)*(test.stat-E.plus)

  # Monte-Carlo simulations of null test statistics
  print("Simulating the null test statistics.")
  null.stat=large.eig.null(n,p1,p2,center=FALSE,samp.num=samp.num,iter=iter)

  # Monte-Carlo based empirical power
  para.pval=mean(test.stat<null.stat)
  # Nonparametric power by comparing the test statistic against the Tracy-Widom law.
  nonpara.pval=1-RMTstat::ptw(test.stat,beta=1)

  return(list(test.stat=test.stat,null.stat=null.stat,para.pval=para.pval,nonpara.pval=nonpara.pval))

}

