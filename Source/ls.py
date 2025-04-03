import numpy as np
import math
import sympy as sp

dx = 1  # Replace with your value
dy = 1  # Replace with your value

matrix = np.array([
    [1, -dx, -dy,  dx*dy, 0.5*dx*dx, 0.5*dy*dy],
    [1,   0, -dy,      0,         0, 0.5*dy*dy],
    [1,  dx, -dy, -dx*dy, 0.5*dx*dx, 0.5*dy*dy],
    [1, -dx,   0,      0, 0.5*dx*dx,         0],
    [1,   0,   0,      0,         0,         0],
    [1,  dx,   0,      0, 0.5*dx*dx,         0],
    [1, -dx,  dy, -dx*dy, 0.5*dx*dx, 0.5*dy*dy],
    [1,   0,  dy,      0,         0, 0.5*dy*dy],
    [1,  dx,  dy,  dx*dy, 0.5*dx*dx, 0.5*dy*dy]
])

print(matrix)

dis = math.sqrt(2)/2
rhs = [-dis, 0, -dis, 0, dis, 0, -dis, 0, -dis]
print(rhs)

# Let us say wanna use the least square method to solver matrix * x = rhs, in whcih x is a 6 by 1 array.

# Solve using least squares
x, residuals, rank, s = np.linalg.lstsq(matrix, rhs, rcond=None)

print("\nSolution (x):")
print(x)

# If you're interested in the residuals
print("\nResiduals:")
print(residuals)

# Enable pretty printing
sp.init_printing(use_unicode=True, wrap_line=False)

# Define symbolic variables
dx, dy = sp.symbols('dx dy')
half = sp.Rational(1,2)

matrixsp = sp.Matrix([
    [1, -dx, -dy,  dx*dy, half*dx*dx, half*dy*dy],
    [1,   0, -dy,      0,          0, half*dy*dy],
    [1,  dx, -dy, -dx*dy, half*dx*dx, half*dy*dy],
    [1, -dx,   0,      0, half*dx*dx,          0],
    [1,   0,   0,      0,          0,          0],
    [1,  dx,   0,      0, half*dx*dx,          0],
    [1, -dx,  dy, -dx*dy, half*dx*dx, half*dy*dy],
    [1,   0,  dy,      0,          0, half*dy*dy],
    [1,  dx,  dy,  dx*dy, half*dx*dx, half*dy*dy]
])

result = matrixsp.transpose() * matrixsp
# sp.pprint(result)

# Compute the inverse of the result matrix
inverse_result = result.inv()

# Display the inverse matrix
# sp.pprint(inverse_result)

# Compute the product of the inverse of the result matrix and matrixsp
product = result.inv() * matrixsp.transpose()
sp.pprint(product)

# Substitute in the values for dx and dy
dx_val = 1  # replace with your value if different
dy_val = 1  # replace with your value if different
result_num = result.subs({dx: dx_val, dy: dy_val})

# Convert the SymPy matrix to a numpy array
result_np = np.array(result_num.tolist(), dtype=float)
print(result_np)

x_cal = np.linalg.inv(result_np) @ matrix.T @ rhs
print(x_cal)
