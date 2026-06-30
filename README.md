kro.inv.test
================
Bongjung Sung
2026-06-29

### kro.inv.test

Suppose $Y_1,\ldots,Y_n$ are i.i.d. $p_1\times p_2$ random matrices with
$V[Y_1]=\Sigma$. Using this data, we test the separability of $\Sigma$,
i.e.,

$$ 
H_0: \Sigma=\Sigma_2\otimes\Sigma_1 \quad \text{versus} \quad H_1:\Sigma\neq \Sigma_2\otimes \Sigma_1,$$
where $\Sigma_1$ and $\Sigma_2$ denote the row and column covariance
matrices, respectively.

Suppose $K^{1/2}CK^{1/2,\top}$ be the Kronecker-core decomposition of
$\Sigma$. Here $K$ and $C$ are referred to as the separable and core
components of $\Sigma$. According to Hoff et al. (2023), $K$ accounts
for the most separable part of $\Sigma$, while the remaining part $C$ is
obtained by whitening $\Sigma$ through $K$. A key insight is that
$\Sigma$ is separable if and only if $C=I_{p_1p_2}$. Thus, testing the
separability of $\Sigma$ is equivalent to testing $$ 
H_0: C=I_p \quad \text{versus} \quad H_1: C\neq I_p$$
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
constituting the covariance matrix with a partial-isotropy rank-$r$ core
are also provided. For details, see Sung and Hoff (2025).

### Reference

Sung, Bongjung and Hoff, Peter (2025). [Testing Separability of
High-Dimensional Covariance Matrices](https://arxiv.org/abs/2506.17463)

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

#

#

# 

#

# 

# 

#

# 

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
