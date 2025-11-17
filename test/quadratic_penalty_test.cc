#include "../core/quad_penalty.h"

using namespace gons;
using namespace ::utilites::LOG_MSG;

using X = ::Vector<double, 2>;

int main() {
  // equality constrained quadratic program problem
  // 0.5 * x^T P x + q^T x
  // s.t.
  // Ax = b
  Matrix<double, 2, 2> P = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q = {1.0, 1.0};
  Matrix<double, 1, 2> A = {{1.0, 1.0}};
  Vector<double, 1> b = {1.0};
  // construct corresponding solver
  penalty_function::QuadPenaltyFunction_ECQP<double, 2, 1> PenaltyFunc(P, q, A, b);
  
    // run optimization                                                                     
  PenaltyFunc.Optimize();
  X rlt = PenaltyFunc.get_result(); // get result
  std::cout << "result: " << rlt << std::endl;
  // inequality constrained quadratic program problem
  // 0.5 * x^T P x + q^T x
  // s.t.
  // l <= Ax <= u [Ax = b => b <= Ax <= b]
  Matrix<double, 2, 2> P2 = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q2 = {1.0, 1.0};
  Matrix<double, 3, 2> A2 = {{1.0, 1.0}, {1.0, 0.0}, {0.0, 1.0}};
  Vector<double, 3> l2 = {1.0, 0.0, 0.0};
  Vector<double, 3> u2 = {1.0, 0.7, 0.7};
  // construct corresponding solver
  penalty_function::QuadPenaltyFunction_IECQP<double, 2, 3> PenaltyFunc2(P2, q2, A2, l2, u2);
  // run optimization
  PenaltyFunc2.Optimize();
  // Before internal unconstrained solver oscillates best approx solution
  // is: 0.299622 0.700024
  // while osqp gives : 0.29877108 0.701228
  // This is the best i can do....
  std::cout << "result: " << PenaltyFunc2.get_result() << std::endl; // get and print the result
  return 0;
}
