#include "matrix2d.h"
#include "matrix2d_impl.h"
#include "gaussj_solver.h"
#include "lu_solver.h"
#include "lsq_solver.h"

using namespace std;

int main() {
    cout << "=== TESTING SOLVERS ===\n" << endl;

    // ================================
    // Test: Gauss-Jordan (vector)
    // ================================
    Matrix2D<double> A1(3, 3, {
        2, -1, 0,
        -1, 2, -1,
        0, -1, 2
        });

    vector<double> b1 = { 1, 0, 1 };

    vector<double> x1_GJ = solveGaussJ(A1, b1);
    cout << "[Gauss-Jordan - vector]: Solution x = ";
    for (double val : x1_GJ) cout << val << " ";
    cout << "\n" << endl;

    // ================================
    // Test: Gauss-Jordan (matrix)
    // ================================
    Matrix2D<double> B1(3, 2, {
        1, 2,
        0, 1,
        1, 0
        });

    Matrix2D<double> X1_GJ = solveGaussJ(A1, B1);
    cout << "[Gauss-Jordan - matrix]: Solution X = " << endl;
    X1_GJ.printMatrix();

    // ================================
    // Test: LU (vector)
    // ================================
    vector<double> x1_LU = solveLU(A1, b1);
    cout << "\n[LU - vector]: Solution x = ";
    for (double val : x1_LU) cout << val << " ";
    cout << "\n" << endl;

    // ================================
    // Test: LU (matrix)
    // ================================
    Matrix2D<double> X1_LU = solveLU(A1, B1);
    cout << "[LU - matrix]: Solution X = " << endl;
    X1_LU.printMatrix();

    // ================================
    // Test: Least Squares (vector)
    // Overdetermined system: 4 equations, 2 unknowns
    // ================================
    Matrix2D<double> A2(4, 2, {
        1, 1,
        2, 1,
        3, 1,
        4, 1
        });

    vector<double> b2 = { 6, 5, 7, 10 };

    vector<double> x2_LSQ = LSQ(A2, b2);
    cout << "\n[Least Squares - vector]: Solution x = ";
    for (double val : x2_LSQ) cout << val << " ";
    cout << "\n" << endl;

    // ================================
    // Test: Least Squares (matrix)
    // ================================
    Matrix2D<double> B2(4, 2, {
        6, 1,
        5, 2,
        7, 3,
        10, 4
        });

    Matrix2D<double> X2_LSQ = LSQ(A2, B2);
    cout << "[Least Squares - matrix]: Solution X = " << endl;
    X2_LSQ.printMatrix();

    // ================================
    // Test: Least Squares Complex
    // ================================
    using cdouble = complex<double>;

    Matrix2D<cdouble> A3(3, 2, {
        {1, 0}, {0, 1},
        {2, 1}, {1, -1},
        {1, -1}, {2, 0}
        });

    vector<cdouble> b3 = {
        {1, 2},
        {2, 1},
        {3, 0}
    };

    vector<cdouble> x3_LSQ = LSQComplex(A3, b3);
    cout << "\n[Least Squares - complex vector]: Solution x = ";
    for (auto z : x3_LSQ) cout << z << " ";
    cout << "\n" << endl;

    cout << "=== END TESTS ===" << endl;

    // Keep the window open until Enter is pressed
    cout << "\nPress Enter to exit...";
    std::cin.get();

    return 0;
}
