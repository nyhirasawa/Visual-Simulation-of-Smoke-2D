#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H
#include "sparse_matrix.h"
#include <vector>
namespace linear_algebra{
//CG法Denseバージョン
void conjugate_gradient(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, int n);
//CG法sparseバージョン
void conjugate_gradient(const sparse_matrix &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr, double eps);

//不完全コレスキー分解
void incomplete_cholesky_decomposition(const sparse_matrix &A, sparse_matrix &L, std::vector<double> &d, const int n);
//不完全コレスキー分解によって得られた下三角行列Lと対角行列(の対角成分を保持したベクトルd)から
//(LDL^{T})^{-1}rを計算する関数. 結果はresultに格納される
void calc_LDLt_inv_r(const sparse_matrix& L, const std::vector<double>& d, const std::vector<double>& r, std::vector<double>& result, int n);
//ICCG法sparseバージョン
void incomplete_cholesky_conjugate_gradient(const sparse_matrix &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr, double eps);


//gauss_seidel法sparseバージョン
void gauss_seidel(const sparse_matrix_with_diagonal_element &A, const std::vector<double> &b, std::vector<double> &x, int n, int max_itr);
}//namespace linear_algebra
#endif
