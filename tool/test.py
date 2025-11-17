from quasi_newton import *
import numpy as np

def my_f_simple(x):
    return x[0]**2 + (x[1]-1)**2

def my_grad_exact(x):
    return np.array([2*x[0], 2*(x[1] - 1) ])
# Specify convergence criteria
eps1 = 1E-8
eps2 = 1E-4

# Create a Lambda (anonymous) function for gradient calculation
# calc_grad = lambda x : my_grad_approx(x,my_f_simple,1E-6,verbose=False)
calc_grad = lambda x : my_grad_exact(x)

# Specify starting point
x0 = np.array([-0.1, 0.2])

# Call optimization routine
x,f,p,B = alg1_sr1(x0,my_f_simple,calc_grad,eps1,eps2,verbose=True,max_iter=50)

# SR1 Hessian approximation
print("\nSR1 Hessian approximation. B[k] =")
analyze_hes(B[-1])

# Actual Hessian
print("True Hessian at x*. B =")
analyze_hes(my_hes_approx(x[-1],calc_grad,1E-6))


