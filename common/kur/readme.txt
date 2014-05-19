These files implement the outlier detection procedure described in the
paper "Robust covariance matrix estimation and multivariate
outlier detection", by D. Peña and F.J. Prieto, Technometrics, 43, 286-300
(2001).

The following (small) modifications have been introduced with respect to the
preceding procedure:

- Different cutoffs have been introduced for the different projection directions
- The number of observations retained after the initial pass of the algorithm is
  allowed to be less than 0.5*(n+p+1), although it should be larger than 2*p
  In no case the number of observations that are finally considered to be outliers
  is allowed to be less than 0.5*(n-p-1)
- A limit has been introduced on the largest size of the data sets to be analyzed,
  to ensure reasonable performances under Matlab
- When suspect outliers are examined to be regrouped, the cutoff is (approximately)
  chosen from a Hotelling-t2 distribution. This improves the type I errors for small
  values of n (the number of observations)

The files required by the procedure are:

kur_rce                      the main file
kur_nw                       a file that computes the orthogonal projection directions
max_kur                      a file to compute one maximization direction
min_kur                      a file to compute one minimization direction
val_kur                      a file that evaluates the projected kurtosis coefficient
kur_rcex                     an alternative main file (see note * below)

The last four files should not be used independently.

Instructions to use the different files are included in the help for each file.
In particular, kur_rce expects as input a matrix having each observation as one of
its rows. The size of the matrix will be the number of observations times the
number of variables. The maximum size has been set to 50 variables and 750 observations,
due to Matlab limitations in handling larger data sets.

All files have been written using Matlab 4.2.

* Note that kur_rce makes use of the following Matlab routines:

    betacdf.m  betainv.m  betapdf.m  distchck.m  finv.m

  These routines are from the Statistics toolbox. If the toolbox is not available, the
  routine kur_rcex could be used (but type I errors will increase for small n).
