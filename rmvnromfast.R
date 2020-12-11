library(Rcpp)
library(RcppArmadillo)
Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")


cores <- parallel::detectCores(logical = FALSE)
sourceCpp("mvrnormc.cpp")


v1 	<- 0.5968
v2	<- 2.3847
cov	<- 0.147

covMat <- matrix(
  c(v1, cov, cov, v2),
  nrow = 2,
  ncol = 2
)

x <- matrix(c(1,2), nrow=1)

library(fastmvtnorm)
dmvnrm_arma_mc(x, mean=c(0,0), sigma = covMat)
fastdmvnorm(x, mean=c(0,0), sigma = covMat)
dmvnrm_arma_mc(x, mean=c(0,0), sigma = covMat)

RE_seq_1 = seq(from = -2 * covMat[1, 1], to = 2 * covMat[1, 1], length.out = 500)
RE_seq_2 = seq(from = -2 * covMat[2, 2], to = 2 * covMat[2, 2], length.out = 500)
RE_W_mat <- outer(X = RE_seq_1, Y = RE_seq_2, FUN = Vectorize(function(x, y) dmvnrm_arma_mc(matrix(c(x, y), nrow=1), mean=c(0,0), sigma = covMat)))




dmvnrm_arma_mc(c(1, 2), mean=c(0,0), sigma = covMat)


x <- matrix(rnorm(10000 * 10), ncol = 10)

x <- matrix (c(RE_seq_1 , RE_seq_2))
Sx <- cov(x)




system.time(RE_W_mat1 <- outer(X = RE_seq_1, Y = RE_seq_2, FUN = Vectorize(function(x, y) fastdmvnorm(matrix(c(x, y), nrow=1), mean=c(0,0), sigma = covMat))))

system.time(RE_W_mat2 <- outer(X = RE_seq_1, Y = RE_seq_2, FUN = Vectorize(function(x, y) binorm_pdf(matrix(c(x, y), nrow=1), sigma = covMat))))

system.time(RE_W_mat3 <- outer(X = RE_seq_1, Y = RE_seq_2, FUN = Vectorize(function(x, y) mvtnorm::dmvnorm(matrix(c(x, y), nrow=1), sigma = covMat))))

