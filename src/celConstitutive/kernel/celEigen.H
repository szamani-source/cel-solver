/*---------------------------------------------------------------------------*\
Function
    cel::constitutive::symmetricEigen

Description
    Eigenvalues and eigenvectors of a symmetric 3x3 matrix by cyclic Jacobi
    rotations (WP 2.5b). Plain C++, no OpenFOAM, no external library.

    Jacobi is used rather than the closed-form (trigonometric) solution
    because it stays accurate to round-off relative to the matrix norm when
    eigenvalues are repeated or nearly repeated (hydrostatic and triaxial
    stress states), where the closed form loses digits. The eigenvectors it
    returns are always orthonormal; for repeated eigenvalues any orthonormal
    basis of the eigenspace is returned, so only V diag(w) V^T is unique.

SourceFiles
    celEigen.C
\*---------------------------------------------------------------------------*/

#ifndef celEigen_H
#define celEigen_H

namespace cel
{
namespace constitutive
{

//- Eigen-decomposition of the symmetric matrix A (only A[i][j], j >= i, is
//  read). Eigenvalues w in descending order, w[0] >= w[1] >= w[2]; column k
//  of V (V[i][k]) is the unit eigenvector of w[k].
void symmetricEigen(const double A[3][3], double w[3], double V[3][3]);

//- Symmetric 3x3 matrix from Voigt stress (xx, yy, zz, xy, yz, zx)
void voigtToMatrix(const double s[6], double A[3][3]);

//- V diag(w) V^T as Voigt stress (xx, yy, zz, xy, yz, zx)
void spectralToVoigt(const double w[3], const double V[3][3], double s[6]);

} // End namespace constitutive
} // End namespace cel

#endif
