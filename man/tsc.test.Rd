\name{tsc.test}
\alias{tsc.test}
\title{Likelihood-ratio Tests for Two-Sample Comparisons}
\description{The functions performs the two-sample comparisons using the following exact test procedures: the exact likelihood-ratio test (LRT) for equality of two normal populations proposed in [4]; the combined test based on the LRT and Shapiro-Wilk (S-W) test for normality via the Bonferroni correction technique; the newly proposed density-based empirical likelihood (DBEL) ratio test. To calculate p-values of the DBEL procedures, three procedures are used: (a) the traditional Monte Carlo (MC) method implemented in C++, (b) a new interpolation method based on regression techniques to operate with tabulated critical values of the test statistic; (c) a Bayesian type method that uses the tabulated critical values as the prior information and  MC generated DBEL-test-statistic's values as data.}
\usage{tsc.test(x,y,method="DBEL",t_m=2,mc=3000)
}
\arguments{
  \item{x}{numeric vector of data for x (missing values are not allowed).}
  \item{y}{numeric vector of data for y (missing values are not allowed).}
  \item{method}{a character string specifying the method for obtaining the test statistic and corresponding p-value. It must be one of "TAS", "TAS&SW" or "DBEL" (default). "TAS" indicates using the exact LRT; "TAS&SW" indicates using the combined test based on the LRT and the S-W test via the Bonferroni correction technique; "DBEL" indicates using the DBEL ratio test. Similar to the DBEL ratio test procedures can be found in [1], [2], [3].}  \item{t_m}{indicates a method for obtaining the p-value when the "DBEL" method is used. It must have values 1, 2 (default), 3, where t_m=1 corresponds to a traditional MC method; t_m=2 corresponds to the interpolation method based on regression techniques and tabulated critical values, this method is similar to that described in [2]; t_m=3 corresponds to a Bayesian type method that combines the method t_m=1 and t_m=2 in a manner similar to that proposed in [5].
}
  \item{mc}{number of monte carlo simulations used to obtain p-value when method="DBEL" and t_m=1 (mc=3000 is default).}

}
\details{The function performs the two-sample comparison using exact procedures: for the LRT to test H_0: X~N, Y~N, E(X)=E(Y), Var(X)=Var(Y) vs. H_1: X~N, Y~N, E(X) is not = E(Y), Var(X) is not = Var(Y); for the LRT combined with the S-W test to test H_0: X~N, Y~N, E(X)=E(Y), Var(X)=Var(Y) vs. H_1: X, Y do not follow normal distributions, E(X) is not = E(Y), Var(X) is not = Var(Y); for the DBEL ratio test to test H_0: X~N, Y~N, E(X)=E(Y), Var(X)=Var(Y) vs. H_1: X, Y do not follow normal distributions, E(X) is not = E(Y), Var(X) is not = Var(Y) (Here X~N means X distributed following a normal distribution).
}
\value{
  Returns a vector of length 2 with a value of the test statistic and the corresponding p-value.
  \item{test_stat}{the value of the test statistic.}
  \item{p_value}{the p-value for the test.}
}
\references{
  [1] Jeffrey C. Miecznikowski, Albert Vexler, Lori A. Shepherd (2013), dbEmpLikeGOF: An R Package for Nonparametric Likelihood-ratio Tests for Goodness-of-Fit and Two-Sample Comparisons Based on Sample Entropy.  Journal of Statistical Software 54(3) 1-19.

  [2] Albert Vexler, Hovig Tanajian, Alan D.Hutson (2014), Density-Based Empirical Likelihood Procedures for Testing Symmetry of Data Distributions and K-Sample Comparisons. The Stata Journal 14(2) 304-328.

  [3] Albert Vexler, Gregory Gurevich (2010), Empirical Likelihood Ratios Applied to Goodness-of-fit Tests Based on Sample Entropy. Computational Statistics & Data Analysis 54(2) 531-545.

  [4] Lingyun Zhang, Xinzhong Xu, Gemai Chen (2012), The Exact Likelihood Ratio Test for Equality of Two Normal Populations. The American Statistician 66(3) 180-184.

  [5] Albert Vexler, Young Min Kim, Jihnhee Yu, Nicole A. Lazar, Alan D. Hutson (2014), Computing Critical Values of Exact Tests by Incorporating Monte Carlo Simulations Combined with Statistical Tables. Scandinavian Journal of Statistics 41(4) 1013-1030.}
\author{Yang Zhao, Albert Vexler, Alan Hutson}\examples{##Ex.1 
x <- rnorm(57,0,1)
y <- rnorm(67,0,1)
##two-sample comparisons test for whether x and y are from normal distributions, 
##  and whether the mean of x is the equal to the mean of y, 
##  and whether the variance of x is the equal to the variance of y.
##  method in [4] is used to obtain the test statistic and corresponding p-value.
test_lrt<-tsc.test(x,y,method="TAS") 
##  combined method based on LRT and S-W via the Bonferroni techinque
##  is used to obtain the p-value. 
test_comb<-tsc.test(x,y,method="TAS&SW") 
##DBEL method is used to obtain the test statistics. 
##Monte carlo method is used to obtain the p-value with 1000 monte carlo simulations.
test_dbel1<-tsc.test(x,y,method="DBEL",t_m=1,mc=1000)
##DBEL method is used to obtain the test statistic. 
##The interpolation method based on the regression technique and tabulated critical values 
##  is used to obtain the p-value.
test_dbel2<-tsc.test(x,y,method="DBEL",t_m=2)
##DBEL method is used to obtain the test statistic. 
##The Bayesian method is used to obtain the p-value.
test_dbel3<-tsc.test(x,y,method="DBEL",t_m=3)

##Ex.2
A<- rnorm(15,0,1)
B<- runif(31,-1,1)
test_lrt1<-tsc.test(A,B,method="TAS") #p-value is 0.3656844.
test_comb1<-tsc.test(A,B,method="TAS&SW") #p-value is 0.02588757.
test_dbel4<-tsc.test(A,B,method="DBEL",t_m=1,mc=1000) #p-value is 0.001.
test_dbel5<-tsc.test(A,B,method="DBEL",t_m=2) #p-value is 0.001774751.
test_dbel6<-tsc.test(A,B,method="DBEL",t_m=3) #p-value is 0.008112455.
##B is not from the normal distribution, so the null hypothesis should be rejected.
##The LRT method does not reject H_0, since this methods work just for X~N and Y~N. 
}
\keyword{methods}