import numpy as np

def sub(A: np.ndarray, b: np.ndarray, type: str = 'lower') -> np.ndarray:
    """
    Solve Ax = b for triangular A with forward/backward substitution
    
    Parameters
    ----------
    A : np.ndarray, shape (M, M)
    b : np.ndarray, shape (M,)
    type : str, default is 'lower'
        'lower' -> if A is lower triangular (forward substitution)
        'upper' -> if A is upper triangular (backward substitution)

    Returns
    -------
    x : np.ndarray, shape (M,)
    """
    x = np.zeros_like(b)

    for i in range(len(b)) if type == 'lower' else reversed(range(len(b))):
        x[i] = (b[i] - np.dot(A[i, :], x[:])) / A[i, i]

    return x
