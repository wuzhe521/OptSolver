#ifndef GONS_QUAD_PENALTY_H
#define GONS_QUAD_PENALTY_H

#include "algebra.h"
#include "config.h"
#include "constant.h"
#include "function.h"
#include "gradientsearch.h"
#include "linearsearch.h"
#include "matrix.h"
#include "newton.h"
#include "quasi_newton.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
namespace penalty_function {
using namespace gons::utilites::LOG_MSG;
using namespace gons::gradientsearch;

// equality constrained quadratic program
// quadratic cost function
// linear constraint
// DataType :  float or double
// Var_Size :  number of variables
// description:
//             f(x) = 0.5 * x' * P * x + q' * x
//             s.t.
//             Ax = b
template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class function_ecqp {
private:
  GONS_FLOAT rho_ = 10.0;
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> b_;

public:
  function_ecqp(Matrix<DataType, Var_Size, Var_Size> &P,
                Vector<DataType, Var_Size> &q,
                Matrix<DataType, Con_Size, Var_Size> &A,
                Vector<DataType, Con_Size> &b)
      : P_(P), q_(q), A_(A), b_(b) {} // Fixed order here
  DataType operator()(Vector<DataType, Var_Size> &x) const {
    return 0.5 * CostFunction(x) + 0.5 * rho_ * penaltyFunction(x);
  }
  Vector<DataType, Var_Size> gradient(Vector<DataType, Var_Size> &x) {
    return CostFunction_grad(x) + rho_ * grad_penaltyFunction(x);
  }
  void update_rho(GONS_FLOAT rho) { rho_ = rho; }

public:
  double CostFunction(Vector<DataType, Var_Size> &x) const {
    return x.dot(P_ * x.transpose()) + q_.dot(x);
  }
  Vector<DataType, Var_Size>
  CostFunction_grad(Vector<DataType, Var_Size> &x) const {
    return x * P_ + q_;
  }
  double penaltyFunction(Vector<DataType, Var_Size> &x) const {
    DataType sum = 0;
    Vector<DataType, Con_Size> Conflict = A_ * x.transpose() - b_;
    for (size_t i = 0; i < Con_Size; i++) {
      sum += Conflict(i) * Conflict(i);
    }
    return sum;
  }
  Vector<DataType, Var_Size> // To be Done  ... !!!
  grad_penaltyFunction(Vector<DataType, Var_Size> &x) {
    Vector<DataType, Con_Size> Ax = A_ * x.transpose();
    Vector<DataType, Con_Size> Conflict = Ax - b_;
    Vector<DataType, Var_Size> grad = Conflict * A_;
    return grad;
  }
};

// inequality constrained quadratic program
// quadratic cost function
// linear constraint
// DataType :  float or double
// Var_Size :  number of variables
// Con_Size :  number of constraints
// description:
//             f(x) = 0.5 * x' * P * x + q' * x
//             s.t.
//             l <= Ax <= b
template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class function_iecqp {
private:
  GONS_FLOAT rho_ = 10.0;
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> l_;
  Vector<DataType, Con_Size> u_;

public:
  function_iecqp(Matrix<DataType, Var_Size, Var_Size> &P,
                 Vector<DataType, Var_Size> &q,
                 Matrix<DataType, Con_Size, Var_Size> &A,
                 Vector<DataType, Con_Size> &l, Vector<DataType, Con_Size> &u)
      : P_(P), q_(q), A_(A), l_(l), u_(u) {}
  DataType operator()(Vector<DataType, Var_Size> &x) const {
    return 0.5 * CostFunction(x) + 0.5 * rho_ * penaltyFunction(x);
  }
  Vector<DataType, Var_Size> gradient(Vector<DataType, Var_Size> &x) {
    return CostFunction_grad(x) + rho_ * grad_penaltyFunction(x);
  }
  DataType CostFunction(Vector<DataType, Var_Size> &x) const {
    return x.dot(P_ * x.transpose()) + q_.dot(x);
  }
  Vector<DataType, Var_Size>
  CostFunction_grad(Vector<DataType, Var_Size> &x) const {
    return x * P_ + q_;
  }
  DataType penaltyFunction(Vector<DataType, Var_Size> &x) const {
    DataType sum = DataType(0);
    Vector<DataType, Con_Size> Ax = A_ * x.transpose();
    for (size_t i = 0; i < Con_Size; i++) {
      sum += std::exp2(std::max(Ax(i) - u_(i), DataType(0))) +
             std::exp2(std::max(l_(i) - Ax(i), DataType(0)));
    }
    return sum;
  }
  Vector<DataType, Var_Size>
  grad_penaltyFunction(Vector<DataType, Var_Size> &x) {
    Vector<DataType, Con_Size> Ax = A_ * x.transpose();
    Vector<DataType, Con_Size> Conflict_l = l_ - Ax;
    Vector<DataType, Con_Size> Conflict_u = Ax - u_;
    Vector<DataType, Var_Size> grad_l;
    Vector<DataType, Var_Size> grad_u;
    for (size_t i = 0; i < Con_Size; i++) {
      if (Conflict_l(i) > 0) {
        grad_l += -1.0 * Conflict_l(i) * A_.get_Row(i); // dont forget -1.0
      }
      if (Conflict_u(i) > 0) {
        grad_u += A_.get_Row(i) * Conflict_u(i);
      }
    }
    return grad_l + grad_u;
  }

  void update_rho(GONS_FLOAT rho) { rho_ = rho; }
};
// quadratic const function
// linear constraint
// DataType :  float or double
// Var_Size :  number of variables
// Con_Size :  number of constraints
// description:
//             f(x) = 0.5 * x' * P * x + q' * x
//             s.t.
//             Ax = b

template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class QuadPenaltyFunction_ECQP {
  using X = Vector<DataType, Var_Size>;
  using F = function_ecqp<DataType, Var_Size, Con_Size>;

public:
  struct quad_penalty_parmeters {
    DataType rho = 1;
    DataType alpha = 2; // penalty growth rate
    GONS_UINT max_iter = 1000;
    GONS_BOOL verbose = false;
    GONS_BOOL war_start = true;
    DataType tol = 1e-3;
  };
  enum class PenaltyOptStatus { Failure, MaxIterationReached, Success };

public:
  // for quadratic function:
  // f(x) = 0.5 * x' * P * x + q' * x
  // for linear inequality constraint:
  //  A * x = b

  QuadPenaltyFunction_ECQP(Matrix<DataType, Var_Size, Var_Size> &P,
                           Vector<DataType, Var_Size> &q,
                           Matrix<DataType, Con_Size, Var_Size> &A,
                           Vector<DataType, Con_Size> &b)
      : P_(P), q_(q), A_(A), b_(b), f_(P, q, A, b) {}
  Vector<DataType, Var_Size> get_result() const { return x_; }
  PenaltyOptStatus Optimize() {

    GONS_UINT iter = 0;
    do {
      quasi_newton::BFGSMethod<F, X> BFGS(f_, x_);
      BFGS.Optimize();
      x_ = BFGS.get_x();
      // 判断是否结束循环
      if (f_.grad_penaltyFunction(x_).Norm2() < parameters_.tol) {
        LOG("满足结束条件")
        LOG("迭代次数为：" << iter)
        LOG("当前解为：" << x_.transpose())
        return PenaltyOptStatus::Success;
      }
      // 不满足结束条件， 继续增大rho
      parameters_.rho *= parameters_.alpha;
      f_.update_rho(parameters_.rho);
      iter++;
    } while (iter < parameters_.max_iter);
    LOG_ERROR("迭代次数超过最大限制")
    return PenaltyOptStatus::MaxIterationReached;
  }

private:
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> b_;
  // optimal variables
  Vector<DataType, Var_Size> x_;

  quad_penalty_parmeters parameters_;
  // Create a function object
  function_ecqp<DataType, Var_Size, Con_Size> f_;
  // unconstrained search method
};
template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class QuadPenaltyFunction_IECQP {
  using X = Vector<DataType, Var_Size>;
  using F = function_iecqp<DataType, Var_Size, Con_Size>;

public:
  struct quad_penalty_parmeters {
    DataType rho = 1;
    DataType alpha = 2.0; // penalty growth rate
    GONS_UINT max_iter = 1000;
    GONS_BOOL verbose = false;
    GONS_BOOL war_start = true;
    DataType tol = 1e-3;
  };
  enum class PenaltyOptStatus {
    Failure,
    InternalUnConSoverMaxStep,
    InternalUnCOnSolverFail,
    MaxIterationReached,
    Success
  };
  QuadPenaltyFunction_IECQP(Matrix<DataType, Var_Size, Var_Size> &P,
                            Vector<DataType, Var_Size> &q,
                            Matrix<DataType, Con_Size, Var_Size> &A,
                            Vector<DataType, Con_Size> &l,
                            Vector<DataType, Con_Size> &u)
      : P_(P), q_(q), A_(A), l_(l), u_(u), f_(P, q, A, l, u) {}
  Vector<DataType, Var_Size> get_result() const { return x_; }
  PenaltyOptStatus Optimize() {
    GONS_UINT iter = 0;
    do {

      gradientsearch::BarzilaiBorwein<F, X> BB(f_, x_);
      auto status = BB.Optimize();
      if (status ==
          gradientsearch::BarzilaiBorwein<F,
                                          X>::BarzilaiBorweinStatus::FAILURE) {
        LOG_ERROR("Internal gradient method Failure")
        return PenaltyOptStatus::InternalUnCOnSolverFail;
      }
      if (status == gradientsearch::BarzilaiBorwein<
                        F, X>::BarzilaiBorweinStatus::MAX_ITERATION_REACHED) {
        LOG_ERROR("Internal gradient method Max Iteration Reached")
        return PenaltyOptStatus::InternalUnConSoverMaxStep;
      }
      Vector<DataType, Var_Size> x_new = BB.get_x();

      // 满足结束条件
      if (f_.grad_penaltyFunction(x_new).Norm2() < parameters_.tol) {
        x_ = x_new;
        LOG("满足结束条件")
        LOG("迭代次数为：" << iter)
        LOG("当前解为：" << x_.transpose())
        break;
      }
      ++iter;
      x_ = x_new;
      LOG("当前rho为：" << parameters_.rho)
      LOG("当前解为：" << x_.transpose())
      parameters_.rho *= parameters_.alpha;
      f_.update_rho(parameters_.rho);
    } while (iter < parameters_.max_iter);
    if (iter >= parameters_.max_iter) {
      LOG_ERROR("迭代次数超过最大限制")
      return PenaltyOptStatus::MaxIterationReached;
    }
    return PenaltyOptStatus::Failure;
  }

private:
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> l_;
  Vector<DataType, Con_Size> u_;
  // optimal variables
  Vector<DataType, Var_Size> x_;

  quad_penalty_parmeters parameters_;
  // Create a function object
  function_iecqp<DataType, Var_Size, Con_Size> f_;
};

} // namespace penalty_function
} // namespace gons

#endif