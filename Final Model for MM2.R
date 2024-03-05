#setting inital value

p=0.99
lambda=1
mu1=1
mu2=0.01

#get those matrices  
a <- c(2*p*mu1,p*mu2,0)
b <- c(2*(1-p)*mu1,1*p*mu1+1*(1-p)*mu2,2*p*mu2)
c <- c(0,1*(1-p)*mu1,2*(1-p)*mu2)
Qbackward <- matrix(c(a, b, c), nrow = length(a), byrow = FALSE)
Qbackward
Qstay<- matrix(c(-(lambda+2*mu1),0,0,0,-(lambda+mu1+mu2),0,0,0,-(lambda+2*mu2)),nrow = 3, ncol = 3, byrow = TRUE)
Qstay
Qforward<-matrix(c(lambda,0,0,0,lambda,0,0,0,lambda),nrow = 3, ncol = 3, byrow = TRUE)
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
  R <- array(0, dim = c(num_iterations, 3, 3))
  R[1,,] <- matrix(0, nrow = 3, ncol = 3)
  
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
#level=2
a2 <- c(p*lambda,0)
b2 <- c((1-p)*lambda,p*lambda)
c2 <- c(0,(1-p)*lambda)
Qforward2 <- matrix(c(a2, b2, c2), nrow = length(a2), byrow = FALSE)
Qforward2
x2<-R_at_i%*%Qbackward
x2
Qstay+x2
inv2<-solve(Qstay+x2)
inv2
R2<- -(Qforward2%*%(inv2))
R2

#level=1
Qbackward1<- matrix(c(2*mu1,0,mu2,mu1,0,2*mu2),nrow = 3, ncol = 2, byrow = TRUE)
Qbackward1
Qstay1<- matrix(c(-(lambda+mu1),0,0,-(lambda+mu2)),nrow = 2, ncol = 2, byrow = TRUE)
Qstay1
Qforward1<-matrix(c(p*lambda,(1-p)*lambda),nrow = 1, ncol = 2, byrow = TRUE)
Qforward1
inv1<-solve(Qstay1+R2%*%Qbackward1)
inv1
R1<- -(Qforward1%*%(inv1))
R1
# Combine results into a list
result_list <- list(R1 = R1,
                    R2 = R2,
                    R=R_at_i)

# Print the results
print("Result List:")
print(result_list)
cat("Result List:\n")
for (name in names(result_list)) {
  cat(name, ":\n")
  print(result_list[[name]], quote = FALSE)
  cat("\n")
}
#start to find the distribution
Qb<-matrix(c(mu1,mu2),nrow = 2, ncol = 1, byrow = TRUE)
Qb 
Qa<-matrix(c(p*lambda,(1-p)*lambda),nrow = 1, ncol = 2, byrow = TRUE)
Qa
Qstay1+R2%*%Qbackward1
xxx<-solve(Qstay1+R2%*%Qbackward1)
pi1n<-(-Qa%*%xxx)
pi1n
pi1n%*%Qb-lambda
identity_matrix <- diag(3)
RI<-identity_matrix-R_at_i
RI
RII<-solve(RI)
RII
pi0<-1/(1+sum(pi1n)+sum(pi1n%*%R2%*%RII))
pi0
pi01<-1/(1+sum(R1)+sum(R1%*%R2%*%RII))
pi01
pi1<-pi0%*%R1
pi1
pi2<-pi1%*%R2
pi2
pi3<-pi2%*%R_at_i
pi3
pi4<-pi3%*%R_at_i
pi4
pi5<-pi4%*%R_at_i
pi5
pi6<-pi5%*%R_at_i
pi6
pi7<-pi6%*%R_at_i
pi7
pi8<-pi7%*%R_at_i
pi8
pi9<-pi8%*%R_at_i
pi9
pi10<-pi9%*%R_at_i
pi10
pi0+sum(pi1)+sum(pi2)+sum(pi3)+sum(pi4)+sum(pi5)+sum(pi6)+sum(pi7)+sum(pi8)+sum(pi9)+sum(pi10)
execpted_value_of_jobs_in_system<-0*pi0+1*sum(pi1)+2*sum(pi2)+3*sum(pi3)+4*sum(pi4)+5*sum(pi5)+6*sum(pi6)+7*sum(pi7)+8*sum(pi8)+9*sum(pi9)+10*sum(pi10)
execpted_value_of_jobs_in_system

#For general M/M/2
rho=lambda/(2*mu1)
rho
pi0_MM2<- 1/(((2*rho)^2)/(2*(1-rho))+1+2*rho)
pi0_MM2
pi1_MM2<-pi0_MM2*2*rho
pi1_MM2
pi2_MM2<-pi0_MM2*((2*rho)^2/2)
pi2_MM2
pi3_MM2<-pi0_MM2*(2*rho^3)
pi3_MM2
pi4_MM2<-pi0_MM2*(2*rho^4)
pi4_MM2
pi5_MM2<-pi0_MM2*(2*rho^5)
pi5_MM2
pi6_MM2<-pi0_MM2*(2*rho^6)
pi6_MM2
pi0_MM2+pi1_MM2+pi2_MM2+pi3_MM2+pi4_MM2+pi5_MM2+pi6_MM2
execpted_value_of_jobs_in_system_general<-0*pi0_MM2+1*pi1_MM2+2*pi2_MM2+3*pi3_MM2+4*pi4_MM2+5*pi5_MM2+6*pi6_MM2
execpted_value_of_jobs_in_system_general
P_wait_general_MM2<- 1-(pi0_MM2+pi1_MM2+pi2_MM2)
P_wait_general_MM2
P_wait_modified_MM2<- 1-(pi0+sum(pi1)+sum(pi2))
P_wait_modified_MM2