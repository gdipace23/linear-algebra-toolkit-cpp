#pragma once

#include "matrix2d.h"     // Matrix2D<T> definition
#include "vector_utils.h" // Includes utility macros (SWAP, ERROR_MESSAGE)

/**************************************************************************************************
 * @brief Solves overdetermined linear systems A·x ≈ b using Least Squares approximation.
 *
 * Solves min ||Ax - b||² using the normal equations:
 *      x = (AᵗA)⁻¹ Aᵗ b
 *
 * @tparam T Numeric type (float, double, etc.)
 * @param A_matrix Rectangular matrix (more rows than columns)
 * @param b_vector Right-hand side vector (same number of rows as A)
 * @return vector<T> Least squares solution vector
 *
 * @note For better numerical stability with large problems, prefer QR decomposition.
 * @ERROR_MESSAGE If dimensions mismatch or matrix is underdetermined (when _CHECKSIZE_ is defined)
 **************************************************************************************************/
template <typename T>
vector<T> LSQ(const Matrix2D<T>& A_matrix, const vector<T>& b_vector) {
#ifdef _CHECKSIZE_
	if (A_matrix.getRows() != int(b_vector.size()))
		ERROR_MESSAGE("LSQ() --> Dimension mismatch: A.rows must equal b.size().");
	if (A_matrix.getRows() < A_matrix.getCols())
		ERROR_MESSAGE("LSQ() --> Underdetermined system: rows(A) < cols(A).");
#endif

	Matrix2D<T> At = A_matrix.Transpose();
	Matrix2D<T> AtA = At * A_matrix;
	Matrix2D<T> AtA_inv = AtA.InverseGaussJ();
	Matrix2D<T> AtA_inv_At = AtA_inv * At;
	return AtA_inv_At * b_vector;
}


/**************************************************************************************************
 * @brief Solves complex overdetermined systems A·x ≈ b using Least Squares.
 *
 * Uses conjugate transpose: x = (AᴴA)⁻¹ Aᴴ b
 *
 * @tparam T Base real type (e.g., double)
 * @param A_matrix Complex coefficient matrix
 * @param b_vector Complex right-hand side vector
 * @return vector<complex<T>> Least squares complex solution vector
 **************************************************************************************************/
template <typename T>
vector<complex<T>> LSQComplex(const Matrix2D<complex<T>>& A_matrix, const vector<complex<T>>& b_vector) {
#ifdef _CHECKSIZE_
	if (A_matrix.getRows() != int(b_vector.size()))
		ERROR_MESSAGE("LSQComplex() --> Dimension mismatch: A.rows must equal b.size().");
	if (A_matrix.getRows() < A_matrix.getCols())
		ERROR_MESSAGE("LSQComplex() --> Underdetermined system: rows(A) < cols(A).");
#endif

	Matrix2D<complex<T>> At = A_matrix.ConjugateTranspose();
	Matrix2D<complex<T>> AtA = At * A_matrix;
	Matrix2D<complex<T>> AtA_inv = AtA.InverseGaussJ();
	Matrix2D<complex<T>> AtA_inv_At = AtA_inv * At;
	return AtA_inv_At * b_vector;
}


/**************************************************************************************************
 * @brief Solves multiple overdetermined systems A·X ≈ B (batched LSQ with multiple RHS).
 *
 * @tparam T Numeric type
 * @param A_matrix Coefficient matrix (MxN)
 * @param b_matrix Right-hand side matrix (MxK)
 * @return Matrix2D<T> Least squares solution matrix (NxK)
 *
 * Each column in `b_matrix` represents a separate LSQ problem.
 **************************************************************************************************/
template <typename T>
Matrix2D<T> LSQ(const Matrix2D<T>& A_matrix, const Matrix2D<T>& b_matrix) {
#ifdef _CHECKSIZE_
	if (A_matrix.getRows() != b_matrix.getRows())
		ERROR_MESSAGE("LSQ() --> A.rows must equal B.rows.");
	if (A_matrix.getRows() < A_matrix.getCols())
		ERROR_MESSAGE("LSQ() --> Underdetermined system: rows(A) < cols(A).");
#endif

	Matrix2D<T> At = A_matrix.Transpose();
	Matrix2D<T> AtA = At * A_matrix;
	Matrix2D<T> AtA_inv = AtA.InverseGaussJ();
	Matrix2D<T> AtA_inv_At = AtA_inv * At;
	return AtA_inv_At * b_matrix;
}


/**************************************************************************************************
 * @brief Solves complex systems A·X ≈ B using least squares with multiple RHS.
 *
 * @tparam T Base real type
 * @param A_matrix Complex coefficient matrix (MxN)
 * @param b_matrix Complex right-hand side matrix (MxK)
 * @return Matrix2D<complex<T>> Least squares complex solution matrix (NxK)
 **************************************************************************************************/
template <typename T>
Matrix2D<complex<T>> LSQComplex(const Matrix2D<complex<T>>& A_matrix, const Matrix2D<complex<T>>& b_matrix) {
#ifdef _CHECKSIZE_
	if (A_matrix.getRows() != b_matrix.getRows())
		ERROR_MESSAGE("LSQComplex() --> A.rows must equal B.rows.");
	if (A_matrix.getRows() < A_matrix.getCols())
		ERROR_MESSAGE("LSQComplex() --> Underdetermined system: rows(A) < cols(A).");
#endif

	Matrix2D<complex<T>> At = A_matrix.ConjugateTranspose();
	Matrix2D<complex<T>> AtA = At * A_matrix;
	Matrix2D<complex<T>> AtA_inv = AtA.InverseGaussJ();
	Matrix2D<complex<T>> AtA_inv_At = AtA_inv * At;
	return AtA_inv_At * b_matrix;
}
