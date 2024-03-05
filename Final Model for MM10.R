#inital values
p=1
lambda=2
mu1=1
mu2=0.5
#level10+
a <- c(10*p*mu1,p*mu2,0,0,0,0,0,0,0,0,0)
b <- c(10*(1-p)*mu1,9*p*mu1+1*(1-p)*mu2,2*p*mu2,0,0,0,0,0,0,0,0)
c <- c(0,9*(1-p)*mu1,8*p*mu1+2*(1-p)*mu2,3*p*mu2,0,0,0,0,0,0,0)
d <- c(0,0,8*(1-p)*mu1,7*p*mu1+3*(1-p)*mu2,4*p*mu2,0,0,0,0,0,0)
e <- c(0,0,0,7*(1-p)*mu1,6*p*mu1+4*(1-p)*mu2,5*p*mu2,0,0,0,0,0)
f <- c(0,0,0,0,6*(1-p)*mu1,5*p*mu1+5*(1-p)*mu2,6*p*mu2,0,0,0,0)
g <- c(0,0,0,0,0,5*(1-p)*mu1,4*p*mu1+6*(1-p)*mu2,7*p*mu2,0,0,0)
h <- c(0,0,0,0,0,0,4*(1-p)*mu1,3*p*mu1+7*(1-p)*mu2,8*p*mu2,0,0)
i <- c(0,0,0,0,0,0,0,3*(1-p)*mu1,2*p*mu1+8*(1-p)*mu2,9*p*mu2,0)
j <- c(0,0,0,0,0,0,0,0,2*(1-p)*mu1,1*p*mu1+9*(1-p)*mu2,10*p*mu2)
k <- c(0,0,0,0,0,0,0,0,0,1*(1-p)*mu1,10*(1-p)*mu2)
Qbackward <- matrix(c(a, b, c, d, e, f, g, h, i, j, k), nrow = length(a), byrow = FALSE)
Qbackward
Qstay<-diag(c(-(lambda+10*mu1+0*mu2), -(lambda+9*mu1+1*mu2), -(lambda+8*mu1+2*mu2), -(lambda+7*mu1+3*mu2), -(lambda+6*mu1+4*mu2), -(lambda+5*mu1+5*mu2), -(lambda+4*mu1+6*mu2), -(lambda+3*mu1+7*mu2), -(lambda+2*mu1+8*mu2), -(lambda+1*mu1+9*mu2), -(lambda+0*mu1+10*mu2)))
Qstay
Qforward<-diag(c(lambda, lambda, lambda, lambda, lambda, lambda, lambda, lambda, lambda, lambda, lambda))
Qforward
library(matlib)
QstayI  <- solve(Qstay)
QstayI
V<- -(Qforward%*%QstayI)
V
W<- -(Qbackward%*%QstayI)
W
solve_recurrence_matrix <- function(V, W, num_iterations) {
  R <- array(0, dim = c(num_iterations, 11, 11))
  R[1,,] <- matrix(0, nrow = 11, ncol = 11)
  for (k in 1:(num_iterations - 1)) {
    R[k + 1,,] <- V + (R[k,,])%*%(R[k,,]) %*% W
  }
  
  result_list <- lapply(1:num_iterations, function(i) R[i,,])
  return(result_list)
}
num_iterations <- 50
result_list <- solve_recurrence_matrix(V, W, num_iterations)
for (i in 1:num_iterations) {
  print(paste("Iteration", i))
  print(result_list[[i]])
}

i <- 50
R_at_i <- result_list[[i]]
print(paste("R matrix at iteration", i))
print(R_at_i)
#level10
Qbackward10<- Qbackward
Qstay10<-Qstay
a10 <- c(p*lambda,0,0,0,0,0,0,0,0,0)
b10 <- c((1-p)*lambda,p*lambda,0,0,0,0,0,0,0,0)
c10 <- c(0,(1-p)*lambda,p*lambda,0,0,0,0,0,0,0)
d10 <- c(0,0,(1-p)*lambda,p*lambda,0,0,0,0,0,0)
e10 <- c(0,0,0,(1-p)*lambda,p*lambda,0,0,0,0,0)
f10 <- c(0,0,0,0,(1-p)*lambda,p*lambda,0,0,0,0)
g10 <- c(0,0,0,0,0,(1-p)*lambda,p*lambda,0,0,0)
h10 <- c(0,0,0,0,0,0,(1-p)*lambda,p*lambda,0,0)
i10 <- c(0,0,0,0,0,0,0,(1-p)*lambda,p*lambda,0)
j10 <- c(0,0,0,0,0,0,0,0,(1-p)*lambda,p*lambda)
k10 <- c(0,0,0,0,0,0,0,0,0,(1-p)*lambda)
Qforward10 <- matrix(c(a10, b10, c10, d10, e10, f10, g10, h10, i10, j10, k10), nrow = length(a10), byrow = FALSE)
Qforward10
inv10<-solve(Qstay10+R_at_i%*%Qbackward10)
inv10
R10<- -(Qforward10%*%(inv10))
R10
#level9
a99 <- c(10*mu1,1*mu2,0,0,0,0,0,0,0,0,0)
b99 <- c(0,9*mu1,2*mu2,0,0,0,0,0,0,0,0)
c99 <- c(0,0,8*mu1,3*mu2,0,0,0,0,0,0,0)
d99 <- c(0,0,0,7*mu1,4*mu2,0,0,0,0,0,0)
e99 <- c(0,0,0,0,6*mu1,5*mu2,0,0,0,0,0)
f99 <- c(0,0,0,0,0,5*mu1,6*mu2,0,0,0,0)
g99 <- c(0,0,0,0,0,0,4*mu1,7*mu2,0,0,0)
h99 <- c(0,0,0,0,0,0,0,3*mu1,8*mu2,0,0)
i99 <- c(0,0,0,0,0,0,0,0,2*mu1,9*mu2,0)
j99 <- c(0,0,0,0,0,0,0,0,0,1*mu1,10*mu2)
Qbackward9 <- matrix(c(a99, b99, c99, d99, e99, f99, g99, h99, i99, j99), nrow = length(a99), byrow = FALSE)
Qbackward9
a9 <- c(p*lambda,0,0,0,0,0,0,0,0)
b9 <- c((1-p)*lambda,p*lambda,0,0,0,0,0,0,0)
c9 <- c(0,(1-p)*lambda,p*lambda,0,0,0,0,0,0)
d9 <- c(0,0,(1-p)*lambda,p*lambda,0,0,0,0,0)
e9 <- c(0,0,0,(1-p)*lambda,p*lambda,0,0,0,0)
f9 <- c(0,0,0,0,(1-p)*lambda,p*lambda,0,0,0)
g9 <- c(0,0,0,0,0,(1-p)*lambda,p*lambda,0,0)
h9 <- c(0,0,0,0,0,0,(1-p)*lambda,p*lambda,0)
i9 <- c(0,0,0,0,0,0,0,(1-p)*lambda,p*lambda)
j9 <- c(0,0,0,0,0,0,0,0,(1-p)*lambda)
Qforward9 <- matrix(c(a9, b9, c9, d9, e9, f9, g9, h9, i9, j9), nrow = length(a9), byrow = FALSE)
Qforward9
Qstay9<-diag(c(-(lambda+9*mu1+0*mu2), -(lambda+8*mu1+1*mu2), -(lambda+7*mu1+2*mu2), -(lambda+6*mu1+3*mu2), -(lambda+5*mu1+4*mu2), -(lambda+4*mu1+5*mu2), -(lambda+3*mu1+6*mu2), -(lambda+2*mu1+7*mu2), -(lambda+1*mu1+8*mu2), -(lambda+0*mu1+9*mu2)))
Qstay9
inv9<-solve(Qstay9+R10%*%Qbackward9)
inv9
R9<- -(Qforward9%*%(inv9))
R9
#level8
a88 <- c(9*mu1,1*mu2,0,0,0,0,0,0,0,0)
b88 <- c(0,8*mu1,2*mu2,0,0,0,0,0,0,0)
c88 <- c(0,0,7*mu1,3*mu2,0,0,0,0,0,0)
d88 <- c(0,0,0,6*mu1,4*mu2,0,0,0,0,0)
e88 <- c(0,0,0,0,5*mu1,5*mu2,0,0,0,0)
f88 <- c(0,0,0,0,0,4*mu1,6*mu2,0,0,0)
g88 <- c(0,0,0,0,0,0,3*mu1,7*mu2,0,0)
h88 <- c(0,0,0,0,0,0,0,2*mu1,8*mu2,0)
i88 <- c(0,0,0,0,0,0,0,0,1*mu1,9*mu2)
Qbackward8 <- matrix(c(a88, b88, c88, d88, e88, f88, g88, h88, i88), nrow = length(a88), byrow = FALSE)
Qbackward8
a8 <- c(p*lambda,0,0,0,0,0,0,0)
b8 <- c((1-p)*lambda,p*lambda,0,0,0,0,0,0)
c8 <- c(0,(1-p)*lambda,p*lambda,0,0,0,0,0)
d8 <- c(0,0,(1-p)*lambda,p*lambda,0,0,0,0)
e8 <- c(0,0,0,(1-p)*lambda,p*lambda,0,0,0)
f8 <- c(0,0,0,0,(1-p)*lambda,p*lambda,0,0)
g8 <- c(0,0,0,0,0,(1-p)*lambda,p*lambda,0)
h8 <- c(0,0,0,0,0,0,(1-p)*lambda,p*lambda)
i8 <- c(0,0,0,0,0,0,0,(1-p)*lambda)
Qforward8 <- matrix(c(a8, b8, c8, d8, e8, f8, g8, h8, i8), nrow = length(a8), byrow = FALSE)
Qforward8
Qstay8<-diag(c(-(lambda+8*mu1+0*mu2), -(lambda+7*mu1+1*mu2), -(lambda+6*mu1+2*mu2), -(lambda+5*mu1+3*mu2), -(lambda+4*mu1+4*mu2), -(lambda+3*mu1+5*mu2), -(lambda+2*mu1+6*mu2), -(lambda+1*mu1+7*mu2), -(lambda+0*mu1+8*mu2)))
Qstay8
inv8<-solve(Qstay8+R9%*%Qbackward8)
inv8
R8<- -(Qforward8%*%(inv8))
R8
#level7
a77 <- c(8*mu1,1*mu2,0,0,0,0,0,0,0)
b77 <- c(0,7*mu1,2*mu2,0,0,0,0,0,0)
c77 <- c(0,0,6*mu1,3*mu2,0,0,0,0,0)
d77 <- c(0,0,0,5*mu1,4*mu2,0,0,0,0)
e77 <- c(0,0,0,0,4*mu1,5*mu2,0,0,0)
f77 <- c(0,0,0,0,0,3*mu1,6*mu2,0,0)
g77 <- c(0,0,0,0,0,0,2*mu1,7*mu2,0)
h77 <- c(0,0,0,0,0,0,0,1*mu1,8*mu2)
Qbackward7 <- matrix(c(a77, b77, c77, d77, e77, f77, g77, h77), nrow = length(a77), byrow = FALSE)
Qbackward7
a7 <- c(p*lambda,0,0,0,0,0,0)
b7 <- c((1-p)*lambda,p*lambda,0,0,0,0,0)
c7 <- c(0,(1-p)*lambda,p*lambda,0,0,0,0)
d7 <- c(0,0,(1-p)*lambda,p*lambda,0,0,0)
e7 <- c(0,0,0,(1-p)*lambda,p*lambda,0,0)
f7 <- c(0,0,0,0,(1-p)*lambda,p*lambda,0)
g7 <- c(0,0,0,0,0,(1-p)*lambda,p*lambda)
h7 <- c(0,0,0,0,0,0,(1-p)*lambda)
Qforward7 <- matrix(c(a7, b7, c7, d7, e7, f7, g7, h7), nrow = length(a7), byrow = FALSE)
Qforward7
Qstay7<-diag(c(-(lambda+7*mu1+0*mu2), -(lambda+6*mu1+1*mu2), -(lambda+5*mu1+2*mu2), -(lambda+4*mu1+3*mu2), -(lambda+3*mu1+4*mu2), -(lambda+2*mu1+5*mu2), -(lambda+1*mu1+6*mu2), -(lambda+0*mu1+7*mu2)))
Qstay7
inv7<-solve(Qstay7+R8%*%Qbackward7)
inv7
R7<- -(Qforward7%*%(inv7))
R7
#level6
a66 <- c(7*mu1,1*mu2,0,0,0,0,0,0)
b66 <- c(0,6*mu1,2*mu2,0,0,0,0,0)
c66 <- c(0,0,5*mu1,3*mu2,0,0,0,0)
d66 <- c(0,0,0,4*mu1,4*mu2,0,0,0)
e66 <- c(0,0,0,0,3*mu1,5*mu2,0,0)
f66 <- c(0,0,0,0,0,2*mu1,6*mu2,0)
g66 <- c(0,0,0,0,0,0,1*mu1,7*mu2)
Qbackward6 <- matrix(c(a66, b66, c66, d66, e66, f66, g66), nrow = length(a66), byrow = FALSE)
Qbackward6
a6 <- c(p*lambda,0,0,0,0,0)
b6 <- c((1-p)*lambda,p*lambda,0,0,0,0)
c6 <- c(0,(1-p)*lambda,p*lambda,0,0,0)
d6 <- c(0,0,(1-p)*lambda,p*lambda,0,0)
e6 <- c(0,0,0,(1-p)*lambda,p*lambda,0)
f6 <- c(0,0,0,0,(1-p)*lambda,p*lambda)
g6 <- c(0,0,0,0,0,(1-p)*lambda)
Qforward6 <- matrix(c(a6, b6, c6, d6, e6, f6, g6), nrow = length(a6), byrow = FALSE)
Qforward6
Qstay6<-diag(c(-(lambda+6*mu1+0*mu2), -(lambda+5*mu1+1*mu2), -(lambda+4*mu1+2*mu2), -(lambda+3*mu1+3*mu2), -(lambda+2*mu1+4*mu2), -(lambda+1*mu1+5*mu2), -(lambda+1*mu1+6*mu2)))
Qstay6
inv6<-solve(Qstay6+R7%*%Qbackward6)
inv6
R6<- -(Qforward6%*%(inv6))
R6
#level5
a55 <- c(6*mu1,1*mu2,0,0,0,0,0)
b55 <- c(0,5*mu1,2*mu2,0,0,0,0)
c55 <- c(0,0,4*mu1,3*mu2,0,0,0)
d55 <- c(0,0,0,3*mu1,4*mu2,0,0)
e55 <- c(0,0,0,0,2*mu1,5*mu2,0)
f55 <- c(0,0,0,0,0,1*mu1,6*mu2)
Qbackward5<- matrix(c(a55, b55, c55, d55, e55, f55), nrow = length(a55), byrow = FALSE)
Qbackward5
Qstay5<-diag(c(-(lambda+5*mu1+0*mu2), -(lambda+4*mu1+1*mu2), -(lambda+3*mu1+2*mu2), -(lambda+2*mu1+3*mu2), -(lambda+1*mu1+4*mu2), -(lambda+0*mu1+5*mu2)))
Qstay5
a5 <- c(p*lambda,0,0,0,0)
b5 <- c((1-p)*lambda,p*lambda,0,0,0)
c5 <- c(0,(1-p)*lambda,p*lambda,0,0)
d5 <- c(0,0,(1-p)*lambda,p*lambda,0)
e5 <- c(0,0,0,(1-p)*lambda,p*lambda)
f5 <- c(0,0,0,0,(1-p)*lambda)
Qforward5 <- matrix(c(a5, b5, c5, d5, e5,f5), nrow = length(a5), byrow = FALSE)
Qforward5
inv5<-solve(Qstay5+R6%*%Qbackward5)
inv5
R5<- -(Qforward5%*%(inv5))
R5
#level=4
a44 <- c(5*mu1,1*mu2,0,0,0,0)
b44 <- c(0,4*mu1,2*mu2,0,0,0)
c44 <- c(0,0,3*mu1,3*mu2,0,0)
d44 <- c(0,0,0,2*mu1,4*mu2,0)
e44 <- c(0,0,0,0,1*mu1,5*mu2)
Qbackward4 <- matrix(c(a44, b44, c44, d44, e44), nrow = length(a44), byrow = FALSE)
Qbackward4
Qstay4<-diag(c(-(lambda+4*mu1+0*mu2), -(lambda+3*mu1+1*mu2), -(lambda+2*mu1+2*mu2), -(lambda+1*mu1+3*mu2), -(lambda+0*mu1+4*mu2)))
Qstay4
a4 <- c(p*lambda,0,0,0)
b4 <- c((1-p)*lambda,p*lambda,0,0)
c4 <- c(0,(1-p)*lambda,p*lambda,0)
d4 <- c(0,0,(1-p)*lambda,p*lambda)
e4 <- c(0,0,0,(1-p)*lambda)
Qforward4 <- matrix(c(a4, b4, c4, d4, e4), nrow = length(a4), byrow = FALSE)
Qforward4
inv4<-solve(Qstay4+R5%*%Qbackward4)
inv4
R4<- -(Qforward4%*%(inv4))
R4
#level=3
a33 <- c(4*mu1,1*mu2,0,0,0)
b33 <- c(0,3*mu1,2*mu2,0,0)
c33 <- c(0,0,2*mu1,3*mu2,0)
d33 <- c(0,0,0,1*mu1,4*mu2)
Qbackward3 <- matrix(c(a33, b33, c33, d33), nrow = length(a33), byrow = FALSE)
Qbackward3
Qstay3<-diag(c(-(lambda+3*mu1+0*mu2), -(lambda+2*mu1+1*mu2), -(lambda+1*mu1+2*mu2), -(lambda+0*mu1+3*mu2)))
Qstay3
a3 <- c(p*lambda,0,0)
b3 <- c((1-p)*lambda,p*lambda,0)
c3 <- c(0,(1-p)*lambda,p*lambda)
d3 <- c(0,0,(1-p)*lambda)
Qforward3 <- matrix(c(a3, b3, c3, d3), nrow = length(a3), byrow = FALSE)
Qforward3
inv3<-solve(Qstay3+R4%*%Qbackward3)
inv3
R3<- -(Qforward3%*%(inv3))
R3
#level=2
a22 <- c(3*mu1,1*mu2,0,0)
b22 <- c(0,2*mu1,2*mu2,0)
c22 <- c(0,0,1*mu1,3*mu2)
Qbackward2 <- matrix(c(a22, b22, c22), nrow = length(a22), byrow = FALSE)
Qbackward2
Qstay2<-diag(c(-(lambda+2*mu1+0*mu2), -(lambda+1*mu1+1*mu2), -(lambda+0*mu1+2*mu2)))
Qstay2
a2 <- c(p*lambda,0)
b2 <- c((1-p)*lambda,p*lambda)
c2 <- c(0,(1-p)*lambda)
Qforward2 <- matrix(c(a2, b2, c2), nrow = length(a2), byrow = FALSE)
Qforward2
inv2<-inv(Qstay2+R3%*%Qbackward2)
inv2
R2<- -(Qforward2%*%(inv2))
R2
#level=1
a11 <- c(2*mu1,1*mu2,0)
b11 <- c(0,1*mu1,2*mu2)
Qbackward1 <- matrix(c(a11, b11), nrow = length(a11), byrow = FALSE)
Qbackward1
Qstay1<-diag(c(-(lambda+1*mu1+0*mu2), -(lambda+0*mu1+1*mu2)))
Qstay1
a1 <- c(p*lambda)
b1 <- c((1-p)*lambda)
Qforward1 <- matrix(c(a1, b1), nrow = length(a1), byrow = FALSE)
Qforward1
inv1<-solve(Qstay1+R2%*%Qbackward1)
inv1
R1<- -(Qforward1%*%(inv1))
R1
identity_matrix <- diag(11)
sum1<-1+sum(R1%*%R2%*%R3%*%R4%*%R5%*%R6%*%R7%*%R8%*%R9%*%R10%*%solve(identity_matrix-R_at_i))+sum(R1)+sum(R1%*%R2)+sum(R1%*%R2%*%R3)+sum(R1%*%R2%*%R3%*%R4)+sum(R1%*%R2%*%R3%*%R4%*%R5)+sum(R1%*%R2%*%R3%*%R4%*%R5%*%R6)+sum(R1%*%R2%*%R3%*%R4%*%R5%*%R6%*%R7)+sum(R1%*%R2%*%R3%*%R4%*%R5%*%R6%*%R7%*%R8)+sum(R1%*%R2%*%R3%*%R4%*%R5%*%R6%*%R7%*%R8%*%R9)
sum1
pi000<-1/sum1
pi000
pi111<-pi000%*%R1
pi111
pi222<-pi111%*%R2
pi222
pi333<-pi222%*%R3
pi333
pi444<-pi333%*%R4
pi444
pi555<-pi444%*%R5
pi555
pi666<-pi555%*%R6
pi666
pi777<-pi666%*%R7
pi777
pi888<-pi777%*%R8
pi888
pi999<-pi888%*%R9
pi999
pi101010<-pi999%*%R10
pi101010
#general model
rho10<-(lambda)/(10*mu1)
k<-10
sum_krho10<- function(k, rho10) {
  sum(sapply(0:(k-1), function(i) (k*rho10)^i / factorial(i)))
}
result<-sum_krho10(k, rho10)
result
result1<-((k*rho10)^k)/(factorial(k)*(1-rho10))
result1
pi0_MM10<-1/(result+result1)
pi0_MM10