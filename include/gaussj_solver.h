#pragma once

#include "matrix2d.h"  // Uses Matrix2D<T> and _CHECKSIZE_

/**************************************************************************************************
 * @brief Solves a linear system A·x = b using Gauss-Jordan elimination with full pivoting.
 *
 * This overload accepts a right-hand side vector `b` and returns a vector of solutions `x`.
 * This method modifies the matrix and RHS internally, so copies are made to preserve originals.
 *
 * @tparam T Numeric type (float, double, complex<T>, etc.)
 * @param A_matrix Coefficient square matrix A (must be NxN)
 * @param b_vector Right-hand side vector b (size N)
 * @return vector<T> Solution vector x
 *
 * @note LU decomposition is preferred when solving multiple systems with the same A.
 * @ERROR_MESSAGE If A is not square or dimensions mismatch (when _CHECKSIZE_ is defined)
 **************************************************************************************************/
template <class T>
vector<T> solveGaussJ(const Matrix2D<T>& A_matrix, const vector<T>& b_vector) {
#ifdef _CHECKSIZE_
	if (A_matrix.getCols() != A_matrix.getRows())
		ERROR_MESSAGE("solveGaussJ() --> Matrix A must be square.");
	if (A_matrix.getCols() != int(b_vector.size()))
		ERROR_MESSAGE("solveGaussJ() --> Dimension mismatch: A.cols must equal b.size.");
#endif

	vector<T> x_result(b_vector);
	Matrix2D<T> A(A_matrix);
	int n = A.getRows();
	vector<int> indxc(n), indxr(n), ipiv(n, 0);
	int irow = 0, icol = 0;

	for (int i = 0; i < n; ++i) {
		// Find pivot
		double big = 0.0;
		for (int j = 0; j < n; ++j)
			if (!ipiv[j])
				for (int k = 0; k < n; ++k)
					if (!ipiv[k] && std::abs(A[j][k]) >= big) {
						big = std::abs(A[j][k]);
						irow = j;
						icol = k;
					}
		++ipiv[icol];

		// Swap rows if needed
		if (irow != icol) {
			for (int l = 0; l < n; ++l) Swap(A[irow][l], A[icol][l]);
			Swap(x_result[irow], x_result[icol]);
		}

		indxr[i] = irow;
		indxc[i] = icol;

		if (A[icol][icol] == static_cast<T>(0.0))
			ERROR_MESSAGE("solveGaussJ() --> Singular matrix encountered.");

		T pivinv = static_cast<T>(1.0) / A[icol][icol];
		A[icol][icol] = static_cast<T>(1.0);
		for (int l = 0; l < n; ++l) A[icol][l] *= pivinv;
		x_result[icol] *= pivinv;

		// Eliminate all other rows
		for (int ll = 0; ll < n; ++ll) {
			if (ll != icol) {
				T dum = A[ll][icol];
				A[ll][icol] = static_cast<T>(0.0);
				for (int l = 0; l < n; ++l) A[ll][l] -= A[icol][l] * dum;
				x_result[ll] -= x_result[icol] * dum;
			}
		}
	}

	return x_result;
}

/**************************************************************************************************
 * @brief Solves a system A·X = B using Gauss-Jordan elimination (multiple RHS).
 *
 * This overload accepts a right-hand side matrix B (N×M), and solves for the matrix X (N×M)
 * such that A·X = B. This is useful when solving for multiple right-hand sides simultaneously.
 *
 * @tparam T Numeric type
 * @param A_matrix Square coefficient matrix A (NxN)
 * @param b_matrix Right-hand side matrix B (NxM)
 * @return Matrix2D<T> Solution matrix X (NxM)
 *
 * @note LU decomposition is more efficient when M << N or for repeated systems.
 * @ERROR_MESSAGE If dimensions do not match or matrix is not square (when _CHECKSIZE_ is defined)
 **************************************************************************************************/
template <class T>
Matrix2D<T> solveGaussJ(const Matrix2D<T>& A_matrix, const Matrix2D<T>& b_matrix) {
#ifdef _CHECKSIZE_
	if (A_matrix.getCols() != A_matrix.getRows())
		ERROR_MESSAGE("solveGaussJ() --> Matrix A must be square.");
	if (A_matrix.getCols() != b_matrix.getRows())
		ERROR_MESSAGE("solveGaussJ() --> Dimension mismatch: A.cols must equal B.rows.");
#endif

	Matrix2D<T> x_result(b_matrix);
	Matrix2D<T> A(A_matrix);
	int n = A.getRows();
	int m = x_result.getCols();

	vector<int> indxc(n), indxr(n), ipiv(n, 0);
	int irow = 0, icol = 0;

	for (int i = 0; i < n; ++i) {
		// Find pivot
		double big = 0.0;
		for (int j = 0; j < n; ++j)
			if (!ipiv[j])
				for (int k = 0; k < n; ++k)
					if (!ipiv[k] && std::abs(A[j][k]) >= big) {
						big = std::abs(A[j][k]);
						irow = j;
						icol = k;
					}
		++ipiv[icol];

		if (irow != icol) {
			for (int l = 0; l < n; ++l) Swap(A[irow][l], A[icol][l]);
			for (int l = 0; l < m; ++l) Swap(x_result[irow][l], x_result[icol][l]);
		}

		indxr[i] = irow;
		indxc[i] = icol;

		if (A[icol][icol] == static_cast<T>(0.0))
			ERROR_MESSAGE("solveGaussJ() --> Singular matrix encountered.");

		T pivinv = static_cast<T>(1.0) / A[icol][icol];
		A[icol][icol] = static_cast<T>(1.0);
		for (int l = 0; l < n; ++l) A[icol][l] *= pivinv;
		for (int l = 0; l < m; ++l) x_result[icol][l] *= pivinv;

		for (int ll = 0; ll < n; ++ll) {
			if (ll != icol) {
				T dum = A[ll][icol];
				A[ll][icol] = static_cast<T>(0.0);
				for (int l = 0; l < n; ++l) A[ll][l] -= A[icol][l] * dum;
				for (int l = 0; l < m; ++l) x_result[ll][l] -= x_result[icol][l] * dum;
			}
		}
	}

	return x_result;
}
