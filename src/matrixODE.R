## utility function to solve a PK system for an arbitrary coefficient
## matrix

## system ODE is given by x' = A * x
## A <- rbind(c(-ka,          0,    0),
##            c( ka, -k10 - k12,  k21),
##            c(  0,  k12      , -k21))

## solution for the homogeneous system is x(t) = exp(At) * x(0)
## with A= U D U^-1 => exp(At) = U exp(Dt) U^-1

## inhomogeneous solution to x' = A * x + b obtained via
## 0  = A  *  xs + b  <=> xs = - A^-1 * b
## x' = A  * (x - xs)
## x  = xs +  exp(At) * (x0 - xs)

matrixODE <- function(A) {
    Ai <- solve(A)
    es <- eigen(A)

    ev <- es$values
    U  <- es$vectors
    Ui <- solve(U)
    zero <- rep(0, ncol(A))

    state <- function(time, x0=zero, b=zero) {
        xs <- - Ai %*% b
        xs + U %*% diag(exp(ev * time)) %*% Ui %*% (x0 - xs)
    }
    Vectorize(state, "time")
}

## allows to return only a specific component from the model
cmt_select <- function(model, cmt) {
    function(...) {
        model(...)[cmt,]
    }
}

## function to setup 2cmt model matrix
model_2cmt <- function(k10,k12,k21) {
  rbind(c( -k10-k12,  k21)
        ,c(     k12, -k21))
}


