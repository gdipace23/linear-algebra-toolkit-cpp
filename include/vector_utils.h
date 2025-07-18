#pragma once

#define _CHECKSIZE_   // Comment this line if you don't want to check the SIZE
#define _CHECKTYPE_   // Comment this line if you don't want to check the TYPE
#define _CHECKBOUNDS_ // Comment this line if you don't want to check the BOUNDS

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <limits>

using std::vector;
using std::complex;
using std::cout;
using std::endl;
using std::numeric_limits;

/************************************************************************************
 * ERROR_MESSAGE Macro
 *
 * This macro emulates basic exception handling in a minimalistic way.
 * It prints a formatted error message including the file and line number
 * where the macro was triggered. The program then waits for user input
 * before terminating execution.
 *
 * This is useful for debugging and simple applications where full exception
 * handling is not necessary or desired.
 *
 * NOTE: In production code, it is recommended to use standard C++ exception
 *       mechanisms (e.g., try-catch blocks) for more flexible error handling.
 ************************************************************************************/
#define ERROR_MESSAGE(msg)                                      \
{                                                               \
    std::cerr << "ERROR: " << msg << "\n"                       \
              << "       in file " << __FILE__                  \
              << ", line " << __LINE__ << "\n";                 \
    std::cerr.flush();                                          \
    std::cerr << "\nPress Enter to exit...";                    \
    std::cin.get();                                             \
    std::exit(EXIT_FAILURE);                                    \
}

 /************************************************************************************
  * Constant Representation of NaN
  *
  * Provides a generic way to obtain a quiet NaN (Not a Number) value for any
  * numeric type T.
  *
  * Example usage:
  *    double x = NaN<double>;
  ************************************************************************************/
template <typename T> static const T NaN = std::numeric_limits<T>::quiet_NaN();

/* === VECTOR OPERATIONS === */

// Overloaded operators
template <typename T> vector<T> operator+(const vector<T>& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator+(const T& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator+(const vector<T>& lhs, const T& rhs);

template <typename T> vector<T> operator-(const vector<T>& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator-(const T& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator-(const vector<T>& lhs, const T& rhs);

template <typename T> vector<T> operator*(const T& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator*(const vector<T>& lhs, const T& rhs);

template <typename T> vector<T> operator/(const T& lhs, const vector<T>& rhs);
template <typename T> vector<T> operator/(const vector<T>& lhs, const T& rhs);

// Vector math
template <typename T> T Dot(const vector<T>& rhs, const vector<T>& lhs);
template <typename T> T Norm(const vector<T>& rhs);
template <typename T> vector<T> Normalize(const vector<T>& rhs);

// I/O
template <typename T> void printVector(const vector<T>& rhs);
template <typename T> void printVector(const vector<T>& rhs, int significantDigits);

// Complex conversion
template <typename T> vector<T> ExpandCmplxVector(const vector<complex<T>>& rhs);
template <typename T> vector<complex<T>> ContractCmplxVector(const vector<T>& rhs);

/************************************************************************************
 * Inline Utility Functions
 *
 * These template functions provide general-purpose utilities similar to common
 * macros in C-style programming, but implemented safely as inline functions.
 *
 * @tparam T  Numeric type (e.g., int, float, double, etc.)
 ************************************************************************************/

 /**
  * @brief Computes the square of a value.
  *
  * @param x Input value.
  * @return x squared.
  */
template <typename T> inline T Square(const T& x) { return x * x; }

/**
 * @brief Returns the maximum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The greater of a and b.
 */
template <typename T> inline T Max(const T& a, const T& b) { return (b > a) ? b : a; }

/**
 * @brief Returns the minimum of two values.
 *
 * @param a First value.
 * @param b Second value.
 * @return The smaller of a and b.
 */
template <typename T> inline T Min(const T& a, const T& b) { return (b < a) ? b : a; }

/**
 * @brief Returns a value with the magnitude of a and the sign of b.
 *
 * @param a Magnitude source.
 * @param b Sign source.
 * @return a with the sign of b.
 */
template <typename T> inline T Sign(const T& a, const T& b) { return (b >= 0) ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a); }

/**
 * @brief Swaps the values of two variables.
 *
 * @param a Reference to the first variable.
 * @param b Reference to the second variable.
 */
template <typename T> inline void Swap(T& a, T& b) { T temp = a; a = b; b = temp; }

/************************************************************************************
 * Arithmetic Operator Overloads for std::vector<T>
 *
 * These overloaded operators allow intuitive element-wise arithmetic between:
 *   - two vectors of the same size
 *   - a scalar and a vector
 *
 * Size checking is enabled if _CHECKSIZE_ is defined.
 ************************************************************************************/

 // vector + vector
template <typename T>
vector<T> operator+(const vector<T>& lhs, const vector<T>& rhs) {
#ifdef _CHECKSIZE_
    if (lhs.size() != rhs.size()) {
        ERROR_MESSAGE("operator+ --> Vector dimensions must agree.");
    }
#endif
    vector<T> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] + rhs[i];
    return result;
}

// scalar + vector
template <typename T>
vector<T> operator+(const T& scalar, const vector<T>& vec) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = scalar + vec[i];
    return result;
}

// vector + scalar
template <typename T>
vector<T> operator+(const vector<T>& vec, const T& scalar) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i] + scalar;
    return result;
}

// vector - vector
template <typename T>
vector<T> operator-(const vector<T>& lhs, const vector<T>& rhs) {
#ifdef _CHECKSIZE_
    if (lhs.size() != rhs.size()) {
        ERROR_MESSAGE("operator- --> Vector dimensions must agree.");
    }
#endif
    vector<T> result(lhs.size());
    for (size_t i = 0; i < lhs.size(); ++i)
        result[i] = lhs[i] - rhs[i];
    return result;
}

// scalar - vector
template <typename T>
vector<T> operator-(const T& scalar, const vector<T>& vec) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = scalar - vec[i];
    return result;
}

// vector - scalar
template <typename T>
vector<T> operator-(const vector<T>& vec, const T& scalar) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i] - scalar;
    return result;
}

// scalar * vector
template <typename T>
vector<T> operator*(const T& scalar, const vector<T>& vec) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = scalar * vec[i];
    return result;
}

// vector * scalar
template <typename T>
vector<T> operator*(const vector<T>& vec, const T& scalar) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i] * scalar;
    return result;
}

// scalar / vector
template <typename T>
vector<T> operator/(const T& scalar, const vector<T>& vec) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = scalar / vec[i];
    return result;
}

// vector / scalar
template <typename T>
vector<T> operator/(const vector<T>& vec, const T& scalar) {
    vector<T> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i)
        result[i] = vec[i] / scalar;
    return result;
}

/******************************************************************************************
 * Dot
 * Computes the scalar (dot) product of two vectors.
 *
 * @param a  First vector operand
 * @param b  Second vector operand
 * @return   Scalar dot product a·b
 *
 * Throws if input vectors have different sizes (when _CHECKSIZE_ is defined).
 ******************************************************************************************/
template <typename T>
T Dot(const vector<T>& a, const vector<T>& b) {
#ifdef _CHECKSIZE_
    if (a.size() != b.size()) {
        ERROR_MESSAGE("Dot() --> Vector dimensions must agree.");
    }
#endif
    T result = static_cast<T>(0.0);
    for (size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

/******************************************************************************************
 * Norm
 * Computes the Euclidean (L2) norm of a vector.
 *
 * @param vec  Input vector
 * @return     Euclidean norm: ||vec||
 ******************************************************************************************/
template <typename T>
T Norm(const vector<T>& vec) {
    T result = static_cast<T>(0.0);
    for (size_t i = 0; i < vec.size(); ++i)
        result += std::abs(vec[i] * vec[i]);
    return std::sqrt(result);
}

/******************************************************************************************
 * Normalize
 * Normalizes a vector to unit Euclidean norm.
 *
 * @param vec  Input vector
 * @return     Normalized vector: vec / ||vec||
 ******************************************************************************************/
template <typename T>
vector<T> Normalize(const vector<T>& vec) {
    T norm = Norm(vec);
    return vec / norm;
}

/******************************************************************************************
 * printVector
 * Prints the contents of a vector to standard output (stdout) with default precision.
 *
 * @param vec  Vector to print
 *
 * Throws if vector is empty.
 ******************************************************************************************/
template <typename T>
void printVector(const vector<T>& vec) {
    if (vec.empty())
        ERROR_MESSAGE("printVector() --> Empty vector.");

    cout.precision(6);
    for (const auto& val : vec)
        cout << val << "\n";
}

/******************************************************************************************
 * printVector (with precision)
 * Prints the contents of a vector with user-defined precision.
 *
 * @param vec                Vector to print
 * @param significantDigits  Number of digits to show after decimal point
 *
 * Throws if vector is empty.
 ******************************************************************************************/
template <typename T>
void printVector(const vector<T>& vec, int significantDigits) {
    if (vec.empty())
        ERROR_MESSAGE("printVector() --> Empty vector.");

    cout.precision(significantDigits);
    for (const auto& val : vec)
        cout << val << "\n";
}

/******************************************************************************************
 * ExpandCmplxVector
 * Expands a complex vector into a real-valued vector by separating real and imaginary parts.
 *
 * Format: [real_0, ..., real_N-1, imag_0, ..., imag_N-1]
 *
 * @param complexVec  Vector of std::complex<T>
 * @return            Real-valued vector of size 2*N
 ******************************************************************************************/
template <typename T>
vector<T> ExpandCmplxVector(const vector<complex<T>>& complexVec) {
    size_t N = complexVec.size();
    vector<T> result(2 * N);

    for (size_t i = 0; i < N; ++i) {
        result[i] = std::real(complexVec[i]);
        result[i + N] = std::imag(complexVec[i]);
    }

    return result;
}

/******************************************************************************************
 * ContractCmplxVector
 * Contracts a real-valued vector into a complex vector, assuming format:
 * [real_0, ..., real_N-1, imag_0, ..., imag_N-1]
 *
 * @param realVec  Real-valued vector of size 2*N
 * @return         Vector of std::complex<T> of size N
 *
 * Throws if input size is not even (when _CHECKSIZE_ is defined).
 ******************************************************************************************/
template <typename T>
vector<complex<T>> ContractCmplxVector(const vector<T>& realVec) {
#ifdef _CHECKSIZE_
    if (realVec.size() % 2 != 0)
        ERROR_MESSAGE("ContractCmplxVector() --> Input vector must have even size.");
#endif

    size_t N = realVec.size() / 2;
    vector<complex<T>> result(N);

    for (size_t i = 0; i < N; ++i)
        result[i] = complex<T>(realVec[i], realVec[i + N]);

    return result;
}