#include <iostream>
#include <cmath>
#include <Eigen/Dense>
void set_as_hilbert(Eigen::MatrixXd & M);
double condition_number(Eigen::MatrixXd & M);

int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  for (int N = 1; N <= 4096; N *= 2) {
    Eigen::MatrixXd M(N, N);
    set_as_hilbert(M);
    std::cout << N << "  " << condition_number(M) << std::endl;
  }
  return 0; }
void set_as_hilbert(Eigen::MatrixXd & M)
{
  M.setZero();
  // Escriba el codigo para iniciar esta matriz como matriz de Hilbert
  for(int n=0;n<M.cols();n++)
    {
      for(int m=0;m<M.cols();m++)
        {
          M(m,n)=1.0/(n+m+1);
        }
    }
}
double condition_number(Eigen::MatrixXd & M)
{
  // Escriba aca el codigo para calcular el condition number.  Puede
  // ser una sola linea!  Se le recomienda buscar como calcular la
  // norma de una matriz de eigen (funcion .norm() ?)
  double k=M.norm()*M.inverse().norm();

     return k;
}
