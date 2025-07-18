#pragma once

#include "matrix2d.h"

/**************************************************************************************************
 * LU Solver Module
 *
 * Implements LU decomposition with implicit pivoting and related routines to:
 *   - Decompose a matrix A into L (lower), U (upper), and permutation matrix P
 *   - Compute the determinant via LU
 *   - Solve systems of equations A*x = b using LU decomposition
 *
 * Dependencies:
 *   - Matrix2D<T>
 *   - vector_utils.h for ERROR_MESSAGE and SWAP
 *
 * Limitations:
 *   - Only works for square matrices
 *   - Does not handle complex numbers (extend as needed)
 **************************************************************************************************/

template <class T> void LUdecomp(const Matrix2D<T>& A, Matrix2D<T>& L, Matrix2D<T>& U, Matrix2D<T>& P);
template <class T> void LUdecomp(const Matrix2D<T>& A, Matrix2D<T>& LU, vector<int>& indx, T& d);
template <class T> T DeterminantLU(const Matrix2D<T>& A);
template <class T> vector<T> solveLU(const Matrix2D<T>& A_matrix, const vector<T>& b_vector);
template <class T> Matrix2D<T> solveLU(const Matrix2D<T>& A_matrix, const Matrix2D<T>& b_matrix);


/***********************************************************************************************
 * LUdecomp (Full form)
 * Performs LU decomposition with implicit pivoting: P*A = L*U
 *
 * @param A   Input matrix (must be square)
 * @param L   Output lower triangular matrix
 * @param U   Output upper triangular matrix
 * @param P   Output permutation matrix (row interchanges)
 ***********************************************************************************************/
template <class T>
void LUdecomp(const Matrix2D<T>& A, Matrix2D<T>& L, Matrix2D<T>& U, Matrix2D<T>& P) {
#ifdef _CHECKSIZE_
    if (A.getCols() != A.getRows()) {
        ERROR_MESSAGE("LUdecomp() --> Matrix A must be square.");
    }
#endif

    Matrix2D<T> LU = A;
    int n = LU.getRows();
    L.Resize(n, n);
    U.Resize(n, n);
    P.Resize(n, n);
    P.SetToIdentity();

    const T TINY = static_cast<T>(1.0e-40);
    vector<double> scale(n);
    for (int i = 0; i < n; i++) {
        double maxVal = 0.0;
        for (int j = 0; j < n; j++)
            maxVal = Max(maxVal, abs(LU[i][j]));
        if (maxVal == 0.0) ERROR_MESSAGE("LUdecomp() --> Singular matrix.");
        scale[i] = 1.0 / maxVal;
    }

    for (int k = 0; k < n; k++) {
        int imax = k;
        double maxPivot = 0.0;
        for (int i = k; i < n; i++) {
            double temp = scale[i] * abs(LU[i][k]);
            if (temp > maxPivot) {
                maxPivot = temp;
                imax = i;
            }
        }

        if (k != imax) {
            for (int j = 0; j < n; j++) {
                SWAP(LU[imax][j], LU[k][j]);
                SWAP(P[imax][j], P[k][j]);
            }
            scale[imax] = scale[k];
        }

        if (LU[k][k] == static_cast<T>(0.0)) LU[k][k] = TINY;

        for (int i = k + 1; i < n; i++) {
            T factor = LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < n; j++)
                LU[i][j] -= factor * LU[k][j];
        }
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                L[i][j] = LU[i][j];
                U[i][j] = static_cast<T>(0.0);
            }
            else if (i == j) {
                L[i][j] = static_cast<T>(1.0);
                U[i][j] = LU[i][j];
            }
            else {
                L[i][j] = static_cast<T>(0.0);
                U[i][j] = LU[i][j];
            }
        }
    }
}


/***********************************************************************************************
 * LUdecomp (Compact form)
 * Performs LU decomposition and stores both L and U in a single matrix (in-place).
 * Suitable for solving systems.
 *
 * @param A     Input matrix (must be square)
 * @param LU    Output matrix with L and U combined
 * @param indx  Row permutations
 * @param d     +1 or -1 based on number of row interchanges
 ***********************************************************************************************/
template <class T>
void LUdecomp(const Matrix2D<T>& A, Matrix2D<T>& LU, vector<int>& indx, T& d) {
#ifdef _CHECKSIZE_
    if (A.getCols() != A.getRows()) {
        ERROR_MESSAGE("LUdecomp() --> Matrix A must be square.");
    }
#endif

    LU = A;
    int n = LU.getRows();
    const T TINY = static_cast<T>(1.0e-40);
    d = static_cast<T>(1.0);
    indx.resize(n);
    vector<double> scale(n);

    for (int i = 0; i < n; i++) {
        double maxVal = 0.0;
        for (int j = 0; j < n; j++)
            maxVal = Max(maxVal, abs(LU[i][j]));
        if (maxVal == 0.0) ERROR_MESSAGE("LUdecomp() --> Singular matrix.");
        scale[i] = 1.0 / maxVal;
    }

    for (int k = 0; k < n; k++) {
        int imax = k;
        double maxPivot = 0.0;
        for (int i = k; i < n; i++) {
            double temp = scale[i] * abs(LU[i][k]);
            if (temp > maxPivot) {
                maxPivot = temp;
                imax = i;
            }
        }

        if (k != imax) {
            for (int j = 0; j < n; j++)
                Swap(LU[k][j], LU[imax][j]);
            d = -d;
            scale[imax] = scale[k];
        }

        indx[k] = imax;
        if (LU[k][k] == static_cast<T>(0.0)) LU[k][k] = TINY;

        for (int i = k + 1; i < n; i++) {
            T factor = LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < n; j++)
                LU[i][j] -= factor * LU[k][j];
        }
    }
}


/***********************************************************************************************
 * DeterminantLU
 * Computes the determinant of a square matrix using LU decomposition.
 *
 * @param A   Input matrix (must be square)
 * @return    Determinant of matrix A
 ***********************************************************************************************/
template <class T>
T DeterminantLU(const Matrix2D<T>& A) {
    Matrix2D<T> LU;
    vector<int> indx;
    T d;
    LUdecomp(A, LU, indx, d);
    for (int i = 0; i < LU.getRows(); i++)
        d *= LU[i][i];
    return d;
}


/**************************************************************************************************

    solveLU

    Solves a system of linear equations A * x = b using LU decomposition.
    This version handles the case where b is a 1D vector.

    NOTE: LU decomposition is generally the preferred method when solving Ax = b with a square A.

    *** INPUTS ***

    A_matrix        Matrix2D<T>&       Coefficient matrix A (must be square)
    b_vector        vector<T>&         Right-hand-side vector b

    *** OUTPUTS ***

    x_result        vector<T>          Solution vector x

**************************************************************************************************/
template <class T>
vector<T> solveLU(const Matrix2D<T>& A_matrix, const vector<T>& b_vector) {
#ifdef _CHECKSIZE_
    if (A_matrix.getCols() != A_matrix.getRows())
        ERROR_MESSAGE("solveLU() --> Matrix A must be square.");

    if (A_matrix.getCols() != static_cast<int>(b_vector.size()))
        ERROR_MESSAGE("solveLU() --> Number of columns in matrix A must match size of vector b.");
#endif

    vector<T> x_result(b_vector);
    Matrix2D<T> LU;
    vector<int> indx;
    T d;

    LUdecomp(A_matrix, LU, indx, d);

    int n = LU.getRows();
    int ii = 0;

    for (int i = 0; i < n; ++i) {
        int ip = indx[i];
        T sum = x_result[ip];
        x_result[ip] = x_result[i];

        if (ii != 0) {
            for (int j = ii - 1; j < i; ++j)
                sum -= LU[i][j] * x_result[j];
        }
        else if (sum != static_cast<T>(0.0)) {
            ii = i + 1;
        }

        x_result[i] = sum;
    }

    for (int i = n - 1; i >= 0; --i) {
        T sum = x_result[i];
        for (int j = i + 1; j < n; ++j)
            sum -= LU[i][j] * x_result[j];

        x_result[i] = sum / LU[i][i];
    }

    return x_result;
}

/**************************************************************************************************

    solveLU

    Solves a system of matrix equations A * X = B using LU decomposition.
    This version handles the case where B is a matrix with multiple right-hand sides.

    NOTE: LU decomposition is preferred when B has many columns (e.g., multiple right-hand-sides),
    especially if m (columns of B) << n.

    *** INPUTS ***

    A_matrix        Matrix2D<T>&       Coefficient matrix A (must be square)
    b_matrix        Matrix2D<T>&       Right-hand-side matrix B

    *** OUTPUTS ***

    x_result        Matrix2D<T>        Solution matrix X

**************************************************************************************************/
template <class T>
Matrix2D<T> solveLU(const Matrix2D<T>& A_matrix, const Matrix2D<T>& b_matrix) {
#ifdef _CHECKSIZE_
    if (A_matrix.getCols() != A_matrix.getRows())
        ERROR_MESSAGE("solveLU() --> Matrix A must be square.");

    if (A_matrix.getCols() != b_matrix.getRows())
        ERROR_MESSAGE("solveLU() --> Number of columns in A must equal number of rows in B.");
#endif

    Matrix2D<T> x_result(b_matrix);
    Matrix2D<T> LU;
    vector<int> indx;
    T d;

    LUdecomp(A_matrix, LU, indx, d);

    int n = LU.getRows();
    int m = b_matrix.getCols();

    for (int j = 0; j < m; ++j) {
        int ii = 0;

        for (int i = 0; i < n; ++i) {
            int ip = indx[i];
            T sum = x_result[ip][j];
            x_result[ip][j] = x_result[i][j];

            if (ii != 0) {
                for (int k = ii - 1; k < i; ++k)
                    sum -= LU[i][k] * x_result[k][j];
            }
            else if (sum != static_cast<T>(0.0)) {
                ii = i + 1;
            }

            x_result[i][j] = sum;
        }

        for (int i = n - 1; i >= 0; --i) {
            T sum = x_result[i][j];
            for (int k = i + 1; k < n; ++k)
                sum -= LU[i][k] * x_result[k][j];

            x_result[i][j] = sum / LU[i][i];
        }
    }

    return x_result;
}
