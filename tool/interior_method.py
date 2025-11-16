import numpy as np
import matplotlib.pyplot as plt


#  cost function
#  y(5 + x)
def f(x):
    return (5 + x[0]) * x[1]
def grad_f(x):
    return np.array([x[1], (5 + x[0])], dtype=float)
def hess_f(x):
    return np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=float)
def constraint1_func(x):
    x[0]*x[1] -5 -x[3]
def constraint1_grad(x):
    return np.array([x[1], x[19], -1, 0], dtype=float)
def constraint1_hess(x):
    return np.array([[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=float)
def constraint2_func(x):
    return x[0]**2 + x[1]**2 - 20 - x[4]
def constraint2_grad(x):
    return np.array([2*x[0], 2*x[1],  0, -1], dtype=float)
def constraint2_hess(x):
    return np.array([[-2, 0, 0, 0], [0, -2, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=float)

def constraint_grad(x):
    return np.array([constraint1_func(x), constraint2_func(x)])
# initial guess value
x0 = np.array([0, 0, 0, 0, 0])
# initial constraint value
c0 = np.array([constraint1_func(x0), constraint2_func(x0)])



