kro.inv.test
================
Bongjung Sung
2026-06-30

### Kronecker-Invariant Tests for High-Dimensional Separability Testing

Suppose $Y_1,\ldots,Y_n$ are i.i.d. $p_1\times p_2$ random matrices with
$V[Y_1]=\Sigma$. Using this data, we test the separability of $\Sigma$,
i.e.,

$$H_0: \Sigma=\Sigma_2\otimes\Sigma_1 \quad \text{versus} \quad H_1:\Sigma\neq \Sigma_2\otimes \Sigma_1,$$

where $\Sigma_1$ and $\Sigma_2$ denote the row and column covariance
matrices, respectively.

Suppose $K^{1/2}CK^{1/2,\top}$ be the Kronecker-core decomposition of
$\Sigma$. Here $K$ and $C$ are referred to as the separable and core
components of $\Sigma$. According to Hoff et al. (2023), $K$ accounts
for the most separable part of $\Sigma$, while the remaining part $C$ is
obtained by whitening $\Sigma$ through $K$. A key insight is that
$\Sigma$ is separable if and only if $C=I_{p_1p_2}$. Thus, testing the
separability of $\Sigma$ is equivalent to testing

$$H_0: C=I_p \quad \text{versus} \quad H_1: C\neq I_p$$

, motivating the use of the sample core $\hat{C}$ for separability
testing, namely, the core component of the sample covariance matrix $S$.

The empirical spectral distribution of $\hat{C}$ is invariant to the
unknown value of $K$. Thus, any test statistic based on the eigenvalues
of $\hat{C}$ has a distribution satisfying this invariance, referred to
as Kronecker-invariance.

We implement three Kronecker-invariant tests for Gaussian populations as
benchmark cases, based on the largest eigenvalue, the extended LRT for
the sphericity, and the separable expansion of the sample core
$\hat{C}$. For each test, we provide functions that implement
Monte-Carlo simulations of the corresponding test statistics under both
the null and alternative hypotheses. The empirical power is evaluated
based on a Monte-Carlo simulation scheme. The nonparametric power is
also provided for the largest eigenvalue test, whose the asymptotic null
distribution was derived, by comparing the test statistic to its
asymptotic null distribution, which is the Tracy-Widom law.

The functions that randomly generate a list of the parameters
constituting the covariance matrix with a partial-isotropy core are also
provided. For details, see Sung and Hoff (2025).

### Reference

Sung, Bongjung and Hoff, Peter (2025). [Testing Separability of
High-Dimensional Covariance Matrices](https://arxiv.org/abs/2506.17463)

Hoff, Peter, McCormack, Andrew and Zhang, Anru R. (2023). [Core
shrinkage covariance estimation for matrix-variate
data](https://arxiv.org/abs/2207.12484)

### Installation

``` r
# Development version 
devtools::install_github("Seungbongjung/kro.inv.test")
```

### Example

``` r
set.seed(100)

# Monte-Carlo simulations of the test statistic based on the largest eigenvalue of the sample core under the null. 

p1=40; p2=40; n=6400
test.stat=kro.inv.test::large.eig.null(n,p1,p2,center=FALSE)
test.stat
hist(test.stat,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Null Distribution")

# Monte-Carlo simulations of the test statistic based on the largest eigenvalue of the sample core under the alternative.

p1=20; p2=16; n=640
para.list=kro.inv.test::pi.rank2.core(p1,p2,lambda.gen=FALSE)

# local alternative
lambda1=1-(0.5/(320/2+0.5))
sigma1=kro.inv.test::pi.core(para.list,lambda0=lambda1)
test.stat1=kro.inv.test::large.eig.alt(n,p1,p2,sigma=sigma1,center=FALSE)
hist(test.stat1,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Distribution")

# dense alternative
lambda2=1-(1.2/(320/2+1.2))
sigma2=kro.inv.test::pi.core(para.list,lambda0=lambda2)
test.stat2=kro.inv.test::large.eig.alt(n,p1,p2,sigma=sigma2,center=FALSE)
hist(test.stat2,freq=FALSE,breaks=25,xlab="Test Statistic",ylab="Density",main="Monte-Carlo Approximated Distribution")

# Empirical power of the test based on the largest eigenvalue of the sample core under Gaussian populations

p1=20; p2=16; n=640

set.seed(100)
para.list=kro.inv.test::pi.rank2.core(p1,p2,lambda.gen=FALSE)

# local alternative
lambda1=1-(0.5/(320/2+0.5))
sigma1=kro.inv.test::pi.core(para.list,lambda0=lambda1)
kro.inv.test::large.eig.power(n,p1,p2,sigma1,center=FALSE)

# dense alternative
lambda2=1-(0.8/(320/2+0.8))
sigma2=kro.inv.test::pi.core(para.list,lambda0=lambda2)
kro.inv.test::large.eig.power(n,p1,p2,sigma2,center=FALSE)

# Empirical power of the test based on the largest eigenvalue of the sample core with the given data.

p1=30; p2=18; r=4; n=200
p=p1*p2

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)

# local alternative
Sigma1=kro.inv.test::pi.core(para.list,lambda0=0.999)
Sigma.root1=kro.inv.test::sym.root(Sigma1)
dat1=crossprod(Sigma.root1,matrix(rnorm(n*p),ncol=n))
dat1=array(dat1,dim=c(p1,p2,n))
dat1=aperm(dat1,perm=c(3,1,2))
kro.inv.test::large.eig.power.dat(dat1,center=FALSE)

# dense alternative
Sigma2=kro.inv.test::pi.core(para.list,lambda0=0.8)
Sigma.root2=kro.inv.test::sym.root(Sigma2)
dat2=crossprod(Sigma.root2,matrix(rnorm(n*p),ncol=n))
dat2=array(dat2,dim=c(p1,p2,n))
dat2=aperm(dat2,perm=c(3,1,2))
kro.inv.test::large.eig.power.dat(dat2,center=FALSE)

# Monte-Carlo simulations of the test statistic based on the extended LRT, applied to the sample core, under the null. 

p1=16; p2=10; n=100

test.stat=kro.inv.test::elrt.null(n,p1,p2,center=FALSE)
hist(test.stat,xlab="Test Statistic",ylab="Density",breaks=25,freq=FALSE,main="Monte-Carlo Approximated Null Distribution")

# Monte-Carlo simulations of the test statistic based on the extended LRT, applied to the sample core, under the alternative. 

p1=16; p2=10; r=4; n=200

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.95)

test.stat=kro.inv.test::elrt.alt(n,p1,p2,Sigma,center=FALSE)
hist(test.stat,xlab="Test Statistic",ylab="Density",breaks=25,freq=FALSE,main="Monte-Carlo Approximated Distribution")
 
# Empirical power of the test based on the extended LRT, applied to the sample core, under Gaussian populations

p1=20; p2=16; r=4; n=200

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.98)

kro.inv.test::elrt.power(n,p1,p2,Sigma,center=FALSE)

# Empirical power of the test based on the extended LRT, applied to the sample core, with the given data.

p1=20; p2=16; r=4; n=200
p=p1*p2

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.98)
Sigma.root=kro.inv.test::sym.root(Sigma)

dat=crossprod(Sigma.root,matrix(rnorm(n*p),ncol=n))
dat=array(dat,dim=c(p1,p2,n))
dat=aperm(dat,perm=c(3,1,2))
kro.inv.test::elrt.power.dat(dat,center=FALSE)

# Monte-Carlo simulations of the test statistic based on the separable expansion of the sample core under the null. 

p1=20; p2=10; n=300

test.stat=kro.inv.test::sep.exp.null(n,p1,p2,center=FALSE)
hist(test.stat,xlab="Test Statistic",ylab="Density",breaks=25,freq=FALSE,main="Monte-Carlo Approximated Null Distribution")

# Monte-Carlo simulations of the test statistic based on the separable expansion of the sample core under the alternative. 

p1=20; p2=15; r=4; n=100

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.95)

test.stat=kro.inv.test::sep.exp.alt(n,p1,p2,Sigma,center=FALSE)
hist(test.stat,xlab="Test Statistic",ylab="Density",breaks=25,freq=FALSE,main="Monte-Carlo Approximated Distribution")

# Empirical power of the test based on the separable expansion of the sample core under Gaussian populations

p1=20; p2=15; r=4; n=100

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.98)

kro.inv.test::sep.exp.power(n,p1,p2,Sigma,center=FALSE)
 
# Empirical power of the test based on the separable expansion of the sample core with the given data.

p1=20; p2=16; r=5; n=100
p=p1*p2

para.list=kro.inv.test::pi.rank_r.core(p1,p2,r,lambda.gen=FALSE)
Sigma=kro.inv.test::pi.core(para.list,lambda0=0.99)
Sigma.root=kro.inv.test::sym.root(Sigma)

dat=crossprod(Sigma.root,matrix(rnorm(n*p),ncol=n))
dat=array(dat,dim=c(p1,p2,n))
dat=aperm(dat,perm=c(3,1,2))
kro.inv.test::sep.exp.power.dat(dat,center=FALSE)

# Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-1 core.

p1=10
kro.inv.test::pi.rank1.core(p1)

# Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-2 core.

# It must be that p1=p2 or p1-p2|p1 (p1>p2) or p2-p1|p2 (p2>p1).
p1=12; p2=10
kro.inv.test::pi.rank2.core(p1,p2)

# Randomly generate a list of parameters for the covariance matrix  with a partial-isotropy rank-r core.

# It must be that p1/p2+p2/p1<r or p1=p2r or p2=p1r
p1=12; p2=8; r=3
kro.inv.test::pi.rank_r.core(p1,p2,r)

# Assemble the covariance matrix with a partial-isotropy rank-r core given a list of parameters.

# generate a list of parameters for a covariance matrix of a partial-isotropy core.
p1=14; p2=10; r=3
para.list=kro.inv.test::pi.rank_r.core(p1,p2,r)
# assembles the covariance matrix using the above list of parameters.
kro.inv.test::pi.core(para.list)
```
