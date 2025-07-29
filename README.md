# 🧮 Linear Algebra Toolkit (C++)

**Linear Algebra Toolkit** is a modular C++ library for performing numerical operations on vectors and matrices.  
It is designed to be lightweight, extensible, and ideal for educational, scientific, and engineering applications.

The toolkit includes efficient routines for solving linear systems using Gauss-Jordan elimination, LU decomposition, and Least Squares methods — all with full support for complex numbers.

---

## 📌 Key Features

✅ Templated matrix class `Matrix2D<T>` supporting:
- Dynamic allocation and resizing
- Arithmetic operations and transposition
- Determinant and matrix inversion
- Identity and conjugate transpose

✅ Solvers:
- Gauss-Jordan elimination with full pivoting
- LU decomposition with implicit pivoting
- Least Squares (real and complex systems)

✅ Clean and portable code with full complex support  
✅ Easy to extend for your custom numerical needs  
✅ Practical usage examples provided in the `examples/` folder  

---

## 📂 Project Structure

```
LinearAlgebraToolkit/
│
├── include/
│ ├── matrix2d.h # Matrix2D class
│ ├── vector_utils.h # Inline utilities (e.g., MAX, MIN, SWAP, etc.)
│ ├── gaussj_solver.h # Gauss-Jordan elimination
│ ├── lu_solver.h # LU decomposition and solver
│ └── lsq_solver.h # Least squares solver
│
├── examples/
│ └── example_solvers.cpp # Demonstration of all solver methods
│
├── README.md
└── LICENSE (MIT)

```

---

## ⚙️ Requirements

- C++17 or later  
- No external dependencies

---

### 🔁 How to Use

1. Include the necessary headers in your source file:
```cpp
#include "matrix2d.h"
#include "gaussj_solver.h"
#include "lu_solver.h"
#include "lsq_solver.h"
```

2. Compile the example or your own code:

```bash
g++ -std=c++17 examples/example_solvers.cpp -o solvers
```

3. Run the executable:

```bash
./solvers
```

---

## 📌 The example will pause at the end and wait for ENTER to be pressed, so you can review the output.

---

## 📘 Module Overview

### ✔️ matrix2d.h 
- Matrix2D<T> class for general matrix operations

- Methods: Resize(), Assign(), SetToIdentity(), Transpose(), InverseGaussJ(), Determinant()

- Operator overloads: +, -, *, []

- Compatible with std::vector<T>

### ✔️ vector_utils.h
- Inline utility macros and functions:

- Square(a) – square of a value

- Max(a,b) / Min(a,b) – elementwise max/min

- Sign(a,b) – applies the sign of b to a

- Swap(a,b) – in-place swap

### ✔️ gaussj_solver.h
- solveGaussJ(A, b) – Solves A * x = b using Gauss-Jordan elimination with full pivoting

- Supports both vector and matrix right-hand sides

- Robust against poorly conditioned matrices

### ✔️ lu_solver.h
- LUdecomp() – Performs LU decomposition with implicit pivoting

- solveLU(A, b) – Solves A * x = b using LU decomposition

- DeterminantLU() – Computes the determinant via LU

- Supports both scalar and matrix systems

### ✔️ lsq_solver.h
- LSQ(A, b) – Least Squares solver for overdetermined systems

- LSQComplex() – Same but for complex-valued matrices/vectors

- Internally uses normal equations: x = (AᵗA)⁻¹ Aᵗb

---

## 💡 Included Examples

| File                    | Description                                           |
|-----------------------------|-------------------------------------------------------|
| `example_solvers.cpp`  | Demonstrates the use of Gauss-Jordan, LU, and Least Squares solvers with both vector and matrix inputs, real and complex types         |

---

## 🚀 Suggestions
Use this library as a numerical backend for:

Signal processing pipelines

System identification and modeling

Scientific and academic computation

Finite element methods and numerical simulations

All modules are templated and extensible — you can easily add support for additional solvers or operations.

---

## 👤 Author

Developed by **Giuseppe Dipace**.

---

## 📄 License

This project is licensed under the MIT License.  
See `LICENSE.md` for details.

---

## ⭐ Contributing

Pull requests and suggestions are welcome!  
Feel free to open an [issue](https://github.com/gdipace23/linear-algebra-toolkit-cpp/issues) or create a PR.
---
