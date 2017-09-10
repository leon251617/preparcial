
#include <iostream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include<Eigen/Core>
#include <Eigen/Eigenvalues>
void set_H0(Eigen::MatrixXd & M);
void set_X(Eigen::MatrixXd & M);
double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, const double lambda, const int index);
int main(int argc, char **argv)
{
  std::cout.precision(16);
  std::cout.setf(std::ios::scientific);
  double lambda = 0.2;
  for (int N = 2; N <= 1024; N *= 2) 
     {
    Eigen::MatrixXd X(N, N), H0(N, N);
    set_H0(H0);
    set_X(X);
    std::cout << 1.0/N << " " << eigen_energy(H0, X, lambda, 0) << std::endl;
     }
  
}

void set_H0(Eigen::MatrixXd & M)
{
  M.setZero();
  // Escriba aca el codigo que crea la matriz H0
  // Puede guiarse por la funcion set_X
  for(int n=0;n<M.cols();n++)
    {
      for(int m=0;m<M.cols();m++)
        {
          if(n==m)
            {
              M(n,m)+=(n+0.5);
            }
        }
    }
}
void set_X(Eigen::MatrixXd & M)
{
  M.setZero();
  for (int n = 0; n < M.cols(); ++n)
    {
    for (int m = 0; m < M.cols(); ++m)
        {
      if (n == m+1) M(n, m) += std::sqrt((m+1.0)/2.0);
      if (n == m-1) M(n, m) += std::sqrt((m)/2.0);
    }
    
        }
 }
double eigen_energy(Eigen::MatrixXd & H, Eigen::MatrixXd & X, const double lambda, const int index)
{
  // Implemente aca el calculo de los valores propios, usando la libreria eigen
  // - Calculo de Hlambda :
  Eigen::MatrixXd HL=H+lambda*X*X*X*X;
  // - Calculo de los valores propios (y vectores propios) :
  Eigen::EigenSolver<Eigen::MatrixXd> es(HL);
  
  //2
  
  // - Extraer los valores propios al vector Eigen::VectorXd evals, solamente
  Eigen::VectorXd  evals=es.eigenvalues().real();
  //   la parte real (cuando pida los valores propios escriba .eigenvalues.real()) :
  
  // ordenar los valores propios :
   std::sort(evals.data(), evals.data() + evals.size());
  //retornar el valor propio (para el estado base index == 0) : return evals(index);
 
   return evals(index);                                                                                                                      
                                                                                      
 
}
