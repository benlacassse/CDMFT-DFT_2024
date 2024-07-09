import numpy as np

def generate_diagonal_pattern_matrix(n):
    if n == 2:
        liste = [0, 1]
    elif n == 3:
        liste = [0, 2, 3]
    elif n == 4:
        liste = [0, 2, 4, 5]
    elif n == 5:
        liste = [0, 3, 5, 7, 8]
        
    # Initialize an n x n matrix with zeros
    matrix = np.zeros((n, n), dtype=int)

    for i in range(n):
        for j in range(n):
            diff = abs(i - j)
            # Determine the value for matrix[i][j]
            diag_num = min(i, j, n - 1 - i, n - 1 - j) + liste[diff]
            matrix[i, j] = diag_num

    return matrix

# Example usage
n = 4

# Generate the template matrix
template = np.array([[1, 2, 2, 3], [2, 1, 3, 2], [2, 3, 1, 2], [3, 2, 2, 1]])

# Generate the coefficients matrix
coefficients_matrix = generate_diagonal_pattern_matrix(n)

# Reshape template and coefficients matrix for broadcasting
template = template.reshape(1, 1, n, n)
coefficients_matrix = (3 * coefficients_matrix).reshape(n, n, 1, 1)

# Create the 4D matrix by adding 3 * x to each element of the template
matrix_4d = template + coefficients_matrix

print(matrix_4d)
print(matrix_4d.shape)
