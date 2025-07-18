#pragma once

#include <vector>
#include <complex>
#include <iostream>
#include "vector_utils.h"

using std::vector;
using std::complex;
using std::cout;
using std::endl;

/**************************************************************************************************
 * Matrix2D
 *
 * A generic templated class for handling two-dimensional matrices with dynamic allocation.
 * Supports arithmetic operations, resizing, identity matrices, transposition, inversion,
 * and more.
 *
 * @tparam T  Numeric type (e.g., float, double, complex<double>)
 **************************************************************************************************/
template <class T>
class Matrix2D {
private:
	int rows;
	int cols;
	T** data;

public:
	/* === Constructors === */

	Matrix2D(); // Default constructor
	Matrix2D(int numRows, int numCols); // Create matrix of size [numRows x numCols]
	Matrix2D(int numRows, int numCols, const T& a); // Initialize matrix with constant value
	Matrix2D(int numRows, int numCols, const T* a); // Initialize matrix from array
	Matrix2D(int numRows, int numCols, const vector<T>& a); // Initialize matrix from vector
	Matrix2D(const Matrix2D<T>& rhs); // Copy constructor

	~Matrix2D(); // Destructor

	/* === Getters === */

	inline int getRows() const;  // Returns number of rows
	inline int getCols() const;  // Returns number of columns
	inline bool Empty() const;   // Returns true if matrix is empty

	/* === Configuration Methods === */

	void Resize(int numRows, int numCols);                     // Resize matrix (data not preserved)
	void Assign(int numRows, int numCols, const T& a);         // Resize and fill with constant
	void Assign(int numRows, int numCols, const vector<T>& a); // Resize and fill with vector
	void SetToIdentity();                                      // Set matrix as identity

	/* === Matrix Operations === */

	Matrix2D<T> InverseGaussJ() const;              // Inverse via Gauss-Jordan with full pivoting
	Matrix2D<T> Transpose() const;                  // Returns the transpose of the matrix
	Matrix2D<T> ConjugateTranspose() const;         // Returns the conjugate transpose (for complex types)
	Matrix2D<T> SubMatrix(int Row, int Col) const;  // Extract submatrix by removing row and column
	T Determinant() const;                          // Computes the determinant

	/* === Printing === */

	void printMatrix();                      // Print matrix to stdout
	void printMatrix(int significantDigits); // Print matrix with user-defined precision

	/* === Operator Overloads === */

	// Assignment
	Matrix2D<T>& operator=(const Matrix2D<T>& rhs); // Assignment operator

	// Subscript access
	inline T* operator[](const int i);             // Access row i (non-const)
	inline const T* operator[](const int i) const; // Access row i (const)

	// Matrix addition
	Matrix2D<T> operator+(const Matrix2D<T>& rhs) const;
	template <class U> friend Matrix2D<U> operator+(const U& lhs, const Matrix2D<U>& rhs); // Scalar + Matrix
	template <class U> friend Matrix2D<U> operator+(const Matrix2D<U>& lhs, const U& rhs); // Matrix + Scalar

	// Matrix subtraction
	Matrix2D<T> operator-(const Matrix2D<T>& rhs) const;
	template <class U> friend Matrix2D<U> operator-(const U& lhs, const Matrix2D<U>& rhs); // Scalar - Matrix
	template <class U> friend Matrix2D<U> operator-(const Matrix2D<U>& lhs, const U& rhs); // Matrix - Scalar

	// Matrix multiplication
	Matrix2D<T> operator*(const Matrix2D<T>& rhs) const;
	template <class U> friend Matrix2D<U> operator*(const U& lhs, const Matrix2D<U>& rhs);         // Scalar * Matrix
	template <class U> friend Matrix2D<U> operator*(const Matrix2D<U>& lhs, const U& rhs);         // Matrix * Scalar
	template <class U> friend Matrix2D<U> operator*(const vector<U>& lhs, const Matrix2D<U>& rhs); // Row vector * Matrix
	template <class U> friend vector<U> operator*(const Matrix2D<U>& lhs, const vector<U>& rhs);   // Matrix * Column vector

	// Convert matrix to flattened vector (row-major order)
	template <class U> friend vector<U> ConvertToVector(const Matrix2D<U>& rhs);

	// Complex matrix handling
	template <class U> friend Matrix2D<U> ExpandCmplxMatrix(const Matrix2D<complex<U>>& rhs);   // [real, -imag; imag, real]
	template <class U> friend Matrix2D<complex<U>> ContractCmplxMatrix(const Matrix2D<U>& rhs); // Complex from 2x real matrix
};

// ============================================================================================
// Matrix2D<T> Implementation
// ============================================================================================

/**
 * @brief Default constructor. Creates an empty matrix.
 */
template <class T>
Matrix2D<T>::Matrix2D() : rows(0), cols(0), data(nullptr) {}

/**
 * @brief Constructor: Initializes matrix with given size (uninitialized values).
 *
 * @param numRows Number of rows
 * @param numCols Number of columns
 */
template <class T>
Matrix2D<T>::Matrix2D(int numRows, int numCols) : rows(numRows), cols(numCols) {
	if (rows * cols > 0) {
		data = new T * [rows];
		for (int i = 0; i < rows; ++i)
			data[i] = new T[cols];
	}
	else {
		data = nullptr;
	}
}

/**
 * @brief Constructor: Initializes matrix with a constant value.
 *
 * @param numRows Number of rows
 * @param numCols Number of columns
 * @param a Value to assign to all elements
 */
template <class T>
Matrix2D<T>::Matrix2D(int numRows, int numCols, const T& a) : rows(numRows), cols(numCols) {
	if (rows * cols > 0) {
		data = new T * [rows];
		for (int i = 0; i < rows; ++i)
			data[i] = new T[cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				data[i][j] = a;
	}
	else {
		data = nullptr;
	}
}

/**
 * @brief Constructor: Initializes matrix from a linear array (row-major order).
 *
 * @param numRows Number of rows
 * @param numCols Number of columns
 * @param a Pointer to input array (size must be numRows * numCols)
 */
template <class T>
Matrix2D<T>::Matrix2D(int numRows, int numCols, const T* a) : rows(numRows), cols(numCols) {
	if (rows * cols > 0) {
		data = new T * [rows];
		for (int i = 0; i < rows; ++i)
			data[i] = new T[cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				data[i][j] = a[j + i * cols];
	}
	else {
		data = nullptr;
	}
}

/**
 * @brief Constructor: Initializes matrix from a std::vector (row-major order).
 *
 * @param numRows Number of rows
 * @param numCols Number of columns
 * @param a Input vector of size numRows * numCols
 */
template <class T>
Matrix2D<T>::Matrix2D(int numRows, int numCols, const vector<T>& a) : rows(numRows), cols(numCols) {
	if (rows * cols > 0) {
		data = new T * [rows];
		for (int i = 0; i < rows; ++i)
			data[i] = new T[cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				data[i][j] = a[i * cols + j];
	}
	else {
		data = nullptr;
	}
}

/**
 * @brief Copy constructor.
 *
 * @param rhs Matrix to copy
 */
template <class T>
Matrix2D<T>::Matrix2D(const Matrix2D<T>& rhs) : rows(rhs.rows), cols(rhs.cols) {
	if (rows * cols > 0) {
		data = new T * [rows];
		for (int i = 0; i < rows; ++i)
			data[i] = new T[cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				data[i][j] = rhs.data[i][j];
	}
	else {
		data = nullptr;
	}
}

/**
 * @brief Destructor: Frees allocated memory.
 */
template <class T>
Matrix2D<T>::~Matrix2D() {
	if (data != nullptr) {
		for (int i = 0; i < rows; ++i)
			delete[] data[i];
		delete[] data;
		data = nullptr;
	}
}

/**
 * @brief Returns number of rows.
 */
template <class T>
inline int Matrix2D<T>::getRows() const {
	return rows;
}

/**
 * @brief Returns number of columns.
 */
template <class T>
inline int Matrix2D<T>::getCols() const {
	return cols;
}

/**
 * @brief Checks if matrix is empty (null pointer).
 */
template <class T>
inline bool Matrix2D<T>::Empty() const {
	return data == nullptr;
}

/**
 * @brief Resizes the matrix (contents not preserved).
 *
 * @param numRows New number of rows
 * @param numCols New number of columns
 */
template <class T>
void Matrix2D<T>::Resize(int numRows, int numCols) {
	if (rows != numRows || cols != numCols) {
		this->~Matrix2D(); // Clean up existing memory
		rows = numRows;
		cols = numCols;

		if (rows * cols > 0) {
			data = new T * [rows];
			for (int i = 0; i < rows; ++i)
				data[i] = new T[cols];
		}
		else {
			data = nullptr;
		}
	}
}

/**
 * @brief Resizes the matrix and fills with a constant value.
 */
template <class T>
void Matrix2D<T>::Assign(int numRows, int numCols, const T& a) {
	Resize(numRows, numCols);
	if (data) {
		for (int r = 0; r < rows; ++r)
			for (int c = 0; c < cols; ++c)
				data[r][c] = a;
	}
}

/**
 * @brief Resizes the matrix and fills using a flat vector (row-major order).
 */
template <class T>
void Matrix2D<T>::Assign(int numRows, int numCols, const vector<T>& a) {
	Resize(numRows, numCols);
	if (data) {
		for (int r = 0; r < rows; ++r)
			for (int c = 0; c < cols; ++c)
				data[r][c] = a[r * cols + c];
	}
}

/**
 * @brief Converts the matrix into an identity matrix.
 *
 * @ERROR_MESSAGEs if matrix is not square or uninitialized (if _CHECKSIZE_ is defined).
 */
template <class T>
void Matrix2D<T>::SetToIdentity() {
#ifdef _CHECKSIZE_
	if (rows != cols)
		ERROR_MESSAGE("SetToIdentity() --> Matrix must be square.");
	if (!data)
		ERROR_MESSAGE("SetToIdentity() --> Matrix is empty.");
#endif

	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			data[i][j] = (i == j) ? static_cast<T>(1.0) : static_cast<T>(0.0);
}

/**
 * @brief Computes the inverse of the matrix using Gauss-Jordan elimination with full pivoting.
 *
 * @return Inverse of the current matrix.
 *
 * @ERROR_MESSAGEs If the matrix is not square or is singular (when _CHECKSIZE_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::InverseGaussJ() const {
#ifdef _CHECKSIZE_
	if (rows != cols)
		ERROR_MESSAGE("InverseGaussJ() --> Matrix must be square.");
	if (data == nullptr)
		ERROR_MESSAGE("InverseGaussJ() --> Matrix is empty.");
#endif

	Matrix2D<T> invA(*this);  // Copy of original matrix to operate on
	int n = invA.getRows();
	int irow = 0, icol = 0;
	double big;
	T pivinv, dum;

	vector<int> indxc(n), indxr(n), ipiv(n, 0);  // Pivot bookkeeping vectors

	for (int i = 0; i < n; ++i) {
		// Select the pivot element
		big = 0.0;
		for (int j = 0; j < n; ++j) {
			if (ipiv[j] != 1) {
				for (int k = 0; k < n; ++k) {
					if (ipiv[k] == 0 && std::abs(invA[j][k]) >= big) {
						big = std::abs(invA[j][k]);
						irow = j;
						icol = k;
					}
				}
			}
		}

		++ipiv[icol];

		// Swap rows if necessary to position pivot
		if (irow != icol) {
			for (int l = 0; l < n; ++l)
				Swap(invA[irow][l], invA[icol][l]);
		}

		indxr[i] = irow;
		indxc[i] = icol;

		if (invA[icol][icol] == static_cast<T>(0.0))
			ERROR_MESSAGE("InverseGaussJ() --> Singular matrix.");

		// Normalize the pivot row
		pivinv = static_cast<T>(1.0) / invA[icol][icol];
		invA[icol][icol] = static_cast<T>(1.0);
		for (int l = 0; l < n; ++l)
			invA[icol][l] *= pivinv;

		// Eliminate all other rows
		for (int ll = 0; ll < n; ++ll) {
			if (ll != icol) {
				dum = invA[ll][icol];
				invA[ll][icol] = static_cast<T>(0.0);
				for (int l = 0; l < n; ++l)
					invA[ll][l] -= invA[icol][l] * dum;
			}
		}
	}

	// Reverse the column interchanges
	for (int l = n - 1; l >= 0; --l) {
		if (indxr[l] != indxc[l]) {
			for (int k = 0; k < n; ++k)
				Swap(invA[k][indxr[l]], invA[k][indxc[l]]);
		}
	}

	return invA;
}

/**
 * @brief Returns the transpose of the current matrix.
 *
 * @return Transposed matrix (rows become columns and vice versa).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::Transpose() const {
	Matrix2D<T> transpose(cols, rows);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			transpose.data[j][i] = data[i][j];
	return transpose;
}

/**
 * @brief Returns the conjugate transpose of a complex matrix.
 *
 * @return Conjugate transposed matrix.
 *
 * @ERROR_MESSAGEs If the matrix does not contain complex elements (when _CHECKTYPE_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::ConjugateTranspose() const {
#ifdef _CHECKTYPE_
	if (typeid(T) != typeid(complex<float>) &&
		typeid(T) != typeid(complex<double>) &&
		typeid(T) != typeid(complex<long double>)) {
		ERROR_MESSAGE("ConjugateTranspose() --> Matrix is not complex.");
	}
#endif

	Matrix2D<T> conjTran(cols, rows);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			conjTran.data[j][i] = std::conj(data[i][j]);
	return conjTran;
}

/**
 * @brief Extracts the (rows-1)x(cols-1) submatrix excluding the given row and column.
 *
 * @param Row Row to exclude
 * @param Col Column to exclude
 * @return Submatrix with the specified row and column removed.
 *
 * @ERROR_MESSAGEs If Row or Col are out of bounds (when _CHECKBOUNDS_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::SubMatrix(int Row, int Col) const {
#ifdef _CHECKBOUNDS_
	if (Row < 0 || Row >= rows || Col < 0 || Col >= cols)
		ERROR_MESSAGE("SubMatrix() --> Index out of bounds.");
#endif

	Matrix2D<T> subMatrix(rows - 1, cols - 1);
	int sub_i = 0;
	for (int i = 0; i < rows; ++i) {
		if (i == Row) continue;
		int sub_j = 0;
		for (int j = 0; j < cols; ++j) {
			if (j == Col) continue;
			subMatrix[sub_i][sub_j++] = data[i][j];
		}
		++sub_i;
	}
	return subMatrix;
}

/**
 * @brief Computes the determinant of the matrix recursively.
 *
 * @return Scalar determinant value.
 *
 * @ERROR_MESSAGEs If the matrix is not square (when _CHECKSIZE_ is defined).
 */
template <class T>
T Matrix2D<T>::Determinant() const {
#ifdef _CHECKSIZE_
	if (rows != cols)
		ERROR_MESSAGE("Determinant() --> Matrix must be square.");
#endif

	// Base case: 1x1 matrix
	if (rows == 1)
		return data[0][0];

	// Recursive case: expand along first row
	T det = static_cast<T>(0.0);
	T sign = static_cast<T>(1.0);
	for (int j = 0; j < cols; ++j) {
		Matrix2D<T> sub = SubMatrix(0, j);
		det += sign * data[0][j] * sub.Determinant();
		sign = -sign;
	}
	return det;
}

/**
 * @brief Prints the contents of the matrix to standard output using default precision.
 *
 * @ERROR_MESSAGEs If the matrix is empty.
 */
template <class T>
void Matrix2D<T>::printMatrix() {
	if (data == nullptr)
		ERROR_MESSAGE("printMatrix() --> Matrix is empty.");

	cout.precision(6);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j)
			cout << data[i][j] << "\t";
		cout << "\n";
	}
}

/**
 * @brief Prints the contents of the matrix with a user-defined number of significant digits.
 *
 * @param significantDigits Number of digits to print after the decimal point.
 *
 * @ERROR_MESSAGEs If the matrix is empty.
 */
template <class T>
void Matrix2D<T>::printMatrix(int significantDigits) {
	if (data == nullptr)
		ERROR_MESSAGE("printMatrix() --> Matrix is empty.");

	cout.precision(significantDigits);
	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j)
			cout << data[i][j] << "\t";
		cout << "\n";
	}
}

/**
 * @brief Assignment operator.
 *
 * Performs deep copy from another matrix of the same type.
 *
 * @param rhs Matrix to assign from.
 * @return Reference to this matrix after assignment.
 */
template <class T>
Matrix2D<T>& Matrix2D<T>::operator=(const Matrix2D<T>& rhs) {
	if (this != &rhs) {
		if (rows != rhs.rows || cols != rhs.cols) {
			// Deallocate existing memory
			if (data != nullptr) {
				for (int i = 0; i < rows; ++i)
					delete[] data[i];
				delete[] data;
				data = nullptr;
			}
			rows = rhs.rows;
			cols = rhs.cols;

			if (rows * cols > 0) {
				data = new T * [rows];
				for (int i = 0; i < rows; ++i)
					data[i] = new T[cols];
			}
		}

		// Copy values
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				data[i][j] = rhs[i][j];
	}
	return *this;
}

/**
 * @brief Row access operator (non-const).
 *
 * @param i Row index.
 * @return Pointer to the i-th row.
 *
 * @ERROR_MESSAGEs If index is out of bounds (when _CHECKBOUNDS_ is defined).
 */
template <class T>
inline T* Matrix2D<T>::operator[](const int i) {
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= rows)
		ERROR_MESSAGE("operator[] --> Row index out of bounds.");
#endif
	return data[i];
}

/**
 * @brief Row access operator (const).
 *
 * @param i Row index.
 * @return Const pointer to the i-th row.
 *
 * @ERROR_MESSAGEs If index is out of bounds (when _CHECKBOUNDS_ is defined).
 */
template <class T>
inline const T* Matrix2D<T>::operator[](const int i) const {
#ifdef _CHECKBOUNDS_
	if (i < 0 || i >= rows)
		ERROR_MESSAGE("operator[] --> Row index out of bounds.");
#endif
	return data[i];
}

/**
 * @brief Performs element-wise matrix addition.
 *
 * @param rhs Right-hand side matrix.
 * @return Matrix2D<T> Resulting matrix from the addition.
 * @ERROR_MESSAGEs If matrix dimensions do not match (when _CHECKSIZE_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::operator+(const Matrix2D<T>& rhs)const {
#ifdef _CHECKSIZE_
	if (rows != rhs.rows || cols != rhs.cols) {
		ERROR_MESSAGE("\n     operator+ --> Incorrect dimensions for matrix sum.\n"
			"     Matrix dimensions must agree.");
	}
#endif // _CHECKSIZE_

	Matrix2D<T> result(rows, cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result[i][j] = data[i][j] + rhs[i][j];

	return result;
}

/**
 * @brief Adds a scalar value to each element of the matrix.
 *
 * @param lhs Scalar value.
 * @param rhs Matrix operand.
 * @return Matrix2D<U> Resulting matrix after scalar addition.
 */
template <class U>
Matrix2D<U> operator+ (const U& lhs, const Matrix2D<U>& rhs) {

	Matrix2D<U> result(rhs.rows, rhs.cols);
	for (int i = 0; i < rhs.rows; i++)
		for (int j = 0; j < rhs.cols; j++)
			result[i][j] = lhs + rhs[i][j];

	return result;
}

/**
 * @brief Adds a scalar value to each element of the matrix.
 *
 * @param lhs Matrix operand.
 * @param rhs Scalar value.
 * @return Matrix2D<U> Resulting matrix after scalar addition.
 */
template <class U>
Matrix2D<U> operator+ (const Matrix2D<U>& lhs, const U& rhs) {

	Matrix2D<U> result(lhs.rows, lhs.cols);
	for (int i = 0; i < lhs.rows; i++)
		for (int j = 0; j < lhs.cols; j++)
			result[i][j] = lhs[i][j] + rhs;

	return result;
}

/**
 * @brief Performs element-wise matrix subtraction.
 *
 * @param rhs Right-hand side matrix.
 * @return Matrix2D<T> Resulting matrix from the subtraction.
 * @ERROR_MESSAGEs If matrix dimensions do not match (when _CHECKSIZE_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::operator-(const Matrix2D<T>& rhs)const {
#ifdef _CHECKSIZE_
	if (rows != rhs.rows || cols != rhs.cols) {
		ERROR_MESSAGE("\n     operator- --> Incorrect dimensions for matrix subtraction.\n"
			"     Matrix dimensions must agree.");
	}
#endif // _CHECKSIZE_

	Matrix2D<T> result(rows, cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result[i][j] = data[i][j] - rhs[i][j];

	return result;
}

/**
 * @brief Subtracts each matrix element from a scalar value.
 *
 * @param lhs Scalar value.
 * @param rhs Matrix operand.
 * @return Matrix2D<U> Resulting matrix after scalar-matrix subtraction.
 */
template <class U>
Matrix2D<U> operator- (const U& lhs, const Matrix2D<U>& rhs) {

	Matrix2D<U> result(rhs.rows, rhs.cols);
	for (int i = 0; i < rhs.rows; i++)
		for (int j = 0; j < rhs.cols; j++)
			result[i][j] = lhs - rhs[i][j];

	return result;
}

/**
 * @brief Subtracts a scalar value from each element of the matrix.
 *
 * @param lhs Matrix operand.
 * @param rhs Scalar value.
 * @return Matrix2D<U> Resulting matrix after matrix-scalar subtraction.
 */
template <class U>
Matrix2D<U> operator- (const Matrix2D<U>& lhs, const U& rhs) {

	Matrix2D<U> result(lhs.rows, lhs.cols);
	for (int i = 0; i < lhs.rows; i++)
		for (int j = 0; j < lhs.cols; j++)
			result[i][j] = lhs[i][j] - rhs;

	return result;
}

/**
 * @brief Performs matrix-matrix multiplication.
 *
 * @param rhs Right-hand side matrix.
 * @return Matrix2D<T> Result of the matrix product.
 * @ERROR_MESSAGEs If lhs.cols != rhs.rows (when _CHECKSIZE_ is defined).
 */
template <class T>
Matrix2D<T> Matrix2D<T>::operator*(const Matrix2D<T>& rhs)const {
#ifdef _CHECKSIZE_
	if (cols != rhs.rows) {
		ERROR_MESSAGE("\n     operator* --> Incorrect dimensions for matrix multiplication.\n"
			"     Check that the number of columns in the first matrix matches the number of rows in the second matrix");
	}
#endif // _CHECKSIZE_

	Matrix2D<T> result(rows, rhs.cols);
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < rhs.cols; j++) {
			result[i][j] = static_cast<T>(0.0);
			for (int k = 0; k < rhs.rows; k++)
				result[i][j] += data[i][k] * rhs[k][j];
		}
	}
	return result;
}

/**
 * @brief Multiplies each element of the matrix by a scalar from the left.
 *
 * @param lhs Scalar multiplier.
 * @param rhs Matrix operand.
 * @return Matrix2D<U> Resulting matrix.
 */
template <class U>
Matrix2D<U> operator* (const U& lhs, const Matrix2D<U>& rhs) {

	Matrix2D<U> result(rhs.rows, rhs.cols);
	for (int i = 0; i < rhs.rows; i++)
		for (int j = 0; j < rhs.cols; j++)
			result[i][j] = lhs * rhs[i][j];

	return result;
}

/**
 * @brief Performs row-vector and matrix multiplication.
 *
 * The vector is treated as a 1-row matrix.
 *
 * @param lhs Row vector (size must match rhs.rows).
 * @param rhs Matrix operand (must have 1 row).
 * @return Matrix2D<U> Resulting matrix.
 * @ERROR_MESSAGEs If rhs.rows != 1 (when _CHECKSIZE_ is defined).
 */
template <class U>
Matrix2D<U> operator* (const Matrix2D<U>& lhs, const U& rhs) {

	Matrix2D<U> result(lhs.rows, lhs.cols);
	for (int i = 0; i < lhs.rows; i++)
		for (int j = 0; j < lhs.cols; j++)
			result[i][j] = lhs[i][j] * rhs;

	return result;
}

/**
 * @brief Performs row-vector and matrix multiplication.
 *
 * The vector is treated as a 1-row matrix.
 *
 * @param lhs Row vector (size must match rhs.rows).
 * @param rhs Matrix operand (must have 1 row).
 * @return Matrix2D<U> Resulting matrix.
 * @ERROR_MESSAGEs If rhs.rows != 1 (when _CHECKSIZE_ is defined).
 */
template <class U>
Matrix2D<U> operator* (const vector<U>& lhs, const Matrix2D<U>& rhs) {
#ifdef _CHECKSIZE_
	if (rhs.rows != 1) {
		ERROR_MESSAGE("\n     operator* --> Incorrect dimensions for matrix multiplication.\n"
			"     Check that the number of columns in the first matrix matches the number of rows in the second matrix");
	}
#endif // _CHECKSIZE_

	Matrix2D<U> result(int(lhs.size()), rhs.cols);
	for (int i = 0; i < result.rows; i++)
		for (int j = 0; j < result.cols; j++)
			result[i][j] = lhs[i] * rhs[0][j];

	return result;
}

/**
 * @brief Performs matrix and column-vector multiplication.
 *
 * The vector is treated as a column vector.
 *
 * @param lhs Matrix operand.
 * @param rhs Column vector (size must match lhs.cols).
 * @return std::vector<U> Resulting vector.
 * @ERROR_MESSAGEs If lhs.cols != rhs.size() (when _CHECKSIZE_ is defined).
 */
template <class U>
vector<U> operator* (const Matrix2D<U>& lhs, const vector<U>& rhs) {
#ifdef _CHECKSIZE_
	if (lhs.cols != int(rhs.size())) {
		ERROR_MESSAGE("\n     operator* --> Incorrect dimensions for matrix multiplication.\n"
			"     Check that the number of columns in the first matrix matches the number of rows in the second matrix");
	}
#endif // _CHECKSIZE_

	vector<U> result(lhs.rows);
	for (int i = 0; i < lhs.rows; i++) {
		result[i] = static_cast<U>(0.0);
		for (int j = 0; j < lhs.cols; j++)
			result[i] += lhs[i][j] * rhs[j];
	}

	return result;
}

/**
 * @brief Flattens a matrix into a row-major vector.
 *
 * @param rhs Input matrix.
 * @return Flattened std::vector<T> containing all elements.
 */
template <class U>
vector<U> ConvertToVector(const Matrix2D<U>& rhs) {
	int sz = rhs.rows * rhs.cols;
	vector<U> result(sz);
	for (int i = 0; i < rhs.rows; ++i)
		for (int j = 0; j < rhs.cols; ++j)
			result[i * rhs.cols + j] = rhs[i][j];
	return result;
}

/**
 * @brief Expands a complex matrix to a 2N x 2M real matrix.
 *
 * Format:
 * [ real(A)   -imag(A) ]
 * [ imag(A)    real(A) ]
 *
 * @param rhs Complex matrix.
 * @return Expanded real matrix.
 */
template <class U>
Matrix2D<U> ExpandCmplxMatrix(const Matrix2D<complex<U>>& rhs) {
	Matrix2D<U> result(2 * rhs.rows, 2 * rhs.cols);
	for (int i = 0; i < rhs.rows; ++i) {
		for (int j = 0; j < rhs.cols; ++j) {
			result[i][j] = real(rhs[i][j]);
			result[i][j + rhs.cols] = -imag(rhs[i][j]);
			result[i + rhs.rows][j] = imag(rhs[i][j]);
			result[i + rhs.rows][j + rhs.cols] = real(rhs[i][j]);
		}
	}
	return result;
}

/**
 * @brief Contracts a real-valued 2N x 2M matrix into a complex N x M matrix.
 *
 * Format:
 * Input assumed to be:
 * [ real(A)   -imag(A) ]
 * [ imag(A)    real(A) ]
 *
 * @param rhs Real matrix of shape 2N x 2M.
 * @return Contracted complex matrix.
 *
 * @ERROR_MESSAGEs If dimensions are not even (when _CHECKSIZE_ is defined).
 */
template <class U>
Matrix2D<complex<U>> ContractCmplxMatrix(const Matrix2D<U>& rhs) {
#ifdef _CHECKSIZE_
	if (rhs.getRows() % 2 != 0 || rhs.getCols() % 2 != 0)
		ERROR_MESSAGE("ContractCmplxMatrix() --> Matrix dimensions must be even.");
#endif
	int rows = rhs.getRows() / 2;
	int cols = rhs.getCols() / 2;
	Matrix2D<complex<U>> result(rows, cols);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			result[i][j] = complex<U>(rhs[i][j], rhs[i + rows][j]);
	return result;
}
