#setting initial value MM5
p=0.999
lambda=5/2
mu1=1
mu2=0.1
#get those matrices  
a <- c(5*p*mu1,p*mu2,0,0,0,0)
b <- c(5*(1-p)*mu1,4*p*mu1+1*(1-p)*mu2,2*p*mu2,0,0,0)
c <- c(0,4*(1-p)*mu1,3*p*mu1+2*(1-p)*mu2,3*p*mu2,0,0)
d <- c(0,0,3*(1-p)*mu1,2*p*mu1+3*(1-p)*mu2,4*p*mu2,0)
e <- c(0,0,0,2*(1-p)*mu1,1*p*mu1+4*(1-p)*mu2,5*p*mu2)
f <- c(0,0,0,0,1*(1-p)*mu1,5*(1-p)*mu2)
Qbackward <- matrix(c(a, b, c, d, e, f), nrow = length(a), byrow = FALSE)
Qbackward
Qstay<-diag(c(-(lambda+5*mu1+0*mu2), -(lambda+4*mu1+1*mu2), -(lambda+3*mu1+2*mu2), -(lambda+2*mu1+3*mu2), -(lambda+1*mu1+4*mu2), -(lambda+0*mu1+5*mu2)))
Qstay
Qforward<-diag(c(lambda, lambda, lambda, lambda, lambda, lambda))
Qforward
#conduct the iteration
library(matlib)
QstayI  <- solve(Qstay)
QstayI
V<- -(Qforward%*%QstayI)
V
W<- -(Qbackward%*%QstayI)
W
# Function to solve the element-wise recurrence relation with matrices
solve_recurrence_matrix <- function(V, W, num_iterations) {
  # Initialize matrices for R
  R <- array(0, dim = c(num_iterations, 6, 6))
  R[1,,] <- matrix(0, nrow = 6, ncol = 6)
  
  # Iterate to calculate R(k+1) element-wise
  for (k in 1:(num_iterations - 1)) {
    R[k + 1,,] <- V + (R[k,,])%*%(R[k,,]) %*% W
  }
  
  # Return a list of matrices representing each iteration
  result_list <- lapply(1:num_iterations, function(i) R[i,,])
  return(result_list)
}

num_iterations <- 20  # Replace with the desired number of iterations

result_list <- solve_recurrence_matrix(V, W, num_iterations)

# Access each iteration
for (i in 1:num_iterations) {
  print(paste("Iteration", i))
  print(result_list[[i]])
}

i <- 20
R_at_i <- result_list[[i]]
print(paste("R matrix at iteration", i))
print(R_at_i)
#With the R is known, we can input the previous level function to solve
#level=5
Qbackward5<- Qbackward
Qstay5<-Qstay
a5 <- c(p*lambda,0,0,0,0)
b5 <- c((1-p)*lambda,p*lambda,0,0,0)
c5 <- c(0,(1-p)*lambda,p*lambda,0,0)
d5 <- c(0,0,(1-p)*lambda,p*lambda,0)
e5 <- c(0,0,0,(1-p)*lambda,p*lambda)
f5 <- c(0,0,0,0,(1-p)*lambda)
Qforward5 <- matrix(c(a5, b5, c5, d5, e5,f5), nrow = length(a5), byrow = FALSE)
Qforward5
inv5<-solve(Qstay5+R_at_i%*%Qbackward5)
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
#start to find the distribution
identity_matrix <- diag(6)
sum<-1+sum(R1%*%R2%*%R3%*%R4%*%R5%*%solve(identity_matrix-R_at_i))+sum(R1)+sum(R1%*%R2)+sum(R1%*%R2%*%R3)+sum(R1%*%R2%*%R3%*%R4)
sum
pi00<-1/sum
pi00
pi11<-pi00%*%R1
pi11
pi22<-pi11%*%R2
pi22
pi33<-pi22%*%R3
pi33
pi44<-pi33%*%R4
pi44
pi55<-pi44%*%R5
pi55
pi66<-pi55%*%R_at_i
pi66
pi77<-pi66%*%R_at_i
pi77
pi88<-pi77%*%R_at_i
pi88
pi99<-pi88%*%R_at_i
pi99
pi1010<-pi99%*%R_at_i
pi1010
#result
pi00+sum(pi11)+sum(pi22)+sum(pi33)+sum(pi44)+sum(pi55)+sum(pi66)+sum(pi77)+sum(pi88)+sum(pi99)+sum(pi1010)
execpted_value_of_jobs_in_system<-0*pi00+1*sum(pi11)+2*sum(pi22)+3*sum(pi33)+4*sum(pi44)+5*sum(pi55)+6*sum(pi66)+7*sum(pi77)+8*sum(pi88)+9*sum(pi99)+10*sum(pi1010)
execpted_value_of_jobs_in_system

#For general M/M/5(1)
#For general M/M/5(2)
rho=lambda/(5*mu1)
rho
pi0_MM5<- 1/(((5*rho)^5)/(120*(1-rho))+5*rho+1+((5*rho)^2)/2+((5*rho)^3)/6+((5*rho)^4)/24)
pi0_MM5
pi1_MM5<-pi0_MM5*5*rho
pi1_MM5
pi2_MM5<-pi0_MM5*((5*rho)^2/2)
pi2_MM5
pi3_MM5<-pi0_MM5*((5*rho)^3/6)
pi3_MM5
pi4_MM5<-pi0_MM5*((5*rho)^4/24)
pi4_MM5
pi5_MM5<-pi0_MM5*((5*rho)^5/120)
pi5_MM5
P_job_wait_MM5<-1-(pi0_MM5+pi1_MM5+pi2_MM5+pi3_MM5+pi4_MM5+pi5_MM5)
P_job_wait_MM5
execpted_value_of_jobs_in_system1<-0*pi0_MM5+1*pi1_MM5+2*pi2_MM5
execpted_value_of_jobs_in_system1
P_job_wait<- 1-(pi00+sum(pi11)+sum(pi22)+sum(pi33)+sum(pi44)+sum(pi55))
P_job_wait