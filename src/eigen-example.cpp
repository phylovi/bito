#include <Eigen/Core>
#include <iostream>
using namespace Eigen;
using namespace std;
template <typename Derived>
Eigen::Block<Derived> topLeftCorner(MatrixBase<Derived>& m, int rows,
                                    int cols) {
  return Eigen::Block<Derived>(m.derived(), 0, 0, rows, cols);
}
template <typename Derived>
const Eigen::Block<const Derived> topLeftCorner(const MatrixBase<Derived>& m,
                                                int rows, int cols) {
  return Eigen::Block<const Derived>(m.derived(), 0, 0, rows, cols);
}

using EigenMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

int main(int, char**) {
  Matrix4d m = Matrix4d::Identity();
  cout << topLeftCorner(4 * m, 2, 3) << endl;  // calls the const version
  topLeftCorner(m, 2, 3) *= 5;                 // calls the non-const version
  cout << "Now the matrix m is:" << endl << m << endl;

  EigenMatrixXd x(3, 0);
  cout << "Now the matrix x is:" << endl << x << endl;
  cout << "it has:" << x.rows() << "and" << x.cols() << endl;
  return 0;
}
