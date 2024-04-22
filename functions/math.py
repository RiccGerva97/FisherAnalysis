import numpy as np

def matrix_inversion(matrix):
    """Calculates the correlation matrix of matrix m_ij."""

    # Check if the matrix is square
    if len(matrix) != len(matrix[0]):
        raise ValueError("Matrix must be square for inversion.")

    # Get the size of the matrix
    n = len(matrix)

    # Create an identity matrix of the same size
    identity = [[0] * n for _ in range(n)]
    for i in range(n):
        identity[i][i] = 1

    # Perform Gaussian elimination to transform the matrix to the identity matrix
    for i in range(n):
        # Scale the current row so that the diagonal element becomes 1
        factor = matrix[i][i]
        for j in range(n):
            matrix[i][j] /= factor
            identity[i][j] /= factor

        # Subtract the scaled current row from other rows to make other elements in the column 0
        for k in range(n):
            if k != i:
                factor = matrix[k][i]
                for j in range(n):
                    matrix[k][j] -= factor * matrix[i][j]
                    identity[k][j] -= factor * identity[i][j]

    # return identiry                    
    res = []
    for i in range(len(identity)):
        res.append(np.array(identity[i]))
    res = np.array(res)

    return res 


def Hartlap(mat, Nr = 350):
    """Calculates inverse matrix using Hartlap correction.
    Arguments:
    - `mat`: input matrix to invert
    - `Nr`: nuber of realization used o calculated the matrix
    """
    # return (Nr-len(mat)-2)/(Nr-1)*np.linalg.inv(mat)
    return (Nr-len(mat)-2)/(Nr-1)*matrix_inversion(mat)


def gaussian(x, mu, sigma):
    """
    Compute the PDF of a Gaussian distribution.

    Parameters:
    - x: Input values.
    - mu: Mean of the distribution.
    - sigma: Standard deviation of the distribution.

    Returns:
    - PDF values for the input x.
    """
    return 1 / (np.sqrt(2 * sigma**2 * np.pi)) * np.exp(-(x - mu)**2 / (2 * sigma**2))

