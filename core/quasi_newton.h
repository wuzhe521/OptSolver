#ifndef GONS_QUASSINEWTON_H
#define GONS_QUASSINEWTON_H

#include "algebra.h"
#include "function.h"
#include "linearsearch.h"
#include "matrix.h"
#include "utilites.h"
#include "vector.h"

namespace gons {
namespace qussinewton {

using namespace gons::utilites::LOG_MSG;
template <typename Function, typename X> class SRMethod {

public:
  enum class SRMethodStatus { Success, Failure, MaxIterationReached };
  enum class SRMethodType { SR1, SR2 };
  struct SRMethodParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 1000;
    GONS_BOOL verbose = false;
    SRMethodType MethodType = SRMethodType::SR1;
  };

  SRMethod(Function &f, X &x) : m_f(f), m_x(x) {}
  ~SRMethod() = default;

  void set_parameters(const SRMethodParameters &params) { m_params = params; }

  X get_x() const { return m_x; }
  double get_function_value() const { return m_f(m_x); }

private:
  template <typename T, GONS_SIZE Size>
  Matrix<T, Size, Size> init_hessian_approximation(const Vector<T, Size> &x) {
    UNUSED(x); // in case warning
    return Identity<T, Size>();
  }

public:
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> SROneUpdate(const Matrix<T, Size, Size> &B,
                                    const Vector<T, Size> &s,
                                    const Vector<T, Size> &y) {
    Matrix<T, Size, Size> result;
    // Compute the outer product of s and y
    Vector<T, Size> yk_min_BkSk = y - B * s.transpose();
    Vector<T, Size> yk_min_BkSk_T = yk_min_BkSk.transpose();
    T denominator = yk_min_BkSk_T.dot(s);
    // CHECK(FLT_EQUAL(denominator, 0.0), "Denominator is zero");
    Matrix<T, Size, Size> nominator = yk_min_BkSk.outerProduct(yk_min_BkSk);

    result = B + nominator * (1.0 / denominator);
    // Update the Hessian approximation
    return result;
  }

public:
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> SRTwoMethod(const Matrix<T, Size, Size> &H,
                                    const Vector<T, Size> &s,
                                    const Vector<T, Size> &y) {
    Matrix<T, Size, Size> result;
    // Compute the outer product of s and y
    Vector<T, Size> sk_min_Bkyk = s - H * y.transpose();
    Vector<T, Size> sk_min_Bkyk_T = sk_min_Bkyk.transpose();
    T denominator = sk_min_Bkyk_T.dot(y);
    // CHECK(FLT_EQUAL(denominator, 0.0), "Denominator is zero");
    Matrix<T, Size, Size> nominator = sk_min_Bkyk.outerProduct(sk_min_Bkyk);

    result = H + nominator * (1.0 / denominator);
    // Update the Hessian approximation
    return result;
  }

  SRMethodStatus Optimize() {

    GONS_UINT iter = 0;

    auto QuassiH =
        init_hessian_approximation(m_x); // Initial Hessian approximation

    while (true) {
      // calc direnction

      X pk = m_params.MethodType == SRMethodType::SR1
                 ? -1 * QuassiH.Inverse() * m_f.gradient(m_x).transpose()
                 : -1 * QuassiH * m_f.gradient(m_x).transpose();
      // line search to find step length Todo: change to wolfe search
      // update x
      X x_new = m_x + pk;
      // update B
      X sk = x_new - m_x;

      X yk = m_f.gradient(x_new) - m_f.gradient(m_x);

      QuassiH = m_params.MethodType == SRMethodType::SR1
                    ? SROneUpdate(QuassiH, sk, yk)
                    : SRTwoMethod(QuassiH, sk, yk);

      m_x = x_new;
      ++iter;

      if (m_f.gradient(m_x).Norm2() < m_params.tolerance) {

        LOG("SRMethod: Optimization using "
            << (m_params.MethodType == SRMethodType::SR1 ? "SR1" : "SR2")
            << " finished.");
        LOG("SRMethod: After Iterations: " << iter);
        LOG("SRMethod: Final X : " << m_x);
        LOG("SRMethod: Function Function value: " << m_f(m_x) << '\n');

        return SRMethodStatus::Success;
      }
      if (iter++ >= m_params.max_iterations) {
        LOG_ERROR("SRMethod: Maximum number of iterations reached.");

        return SRMethodStatus::MaxIterationReached;
      }
      if (m_params.verbose) {
        LOG("Iteration: " << iter);
        LOG("x: " << m_x);
      }
    }
  }

private:
  Function &m_f;
  X &m_x;
  SRMethodParameters m_params;
  // Line search method Todo: implement different line search methods
  // linearsearch::WolfeSearch<Function, X> m_wolfe_search;
}; // class SRMethod

template <typename Function, typename X> class BFGSMethod {

public:
  enum class ApproxType { Direct, Inverse };
  enum class BFGSMethodStatus { Success, Failure, MaxIterationReached };
  struct BFGSMethodParameters {
    double tolerance = 1.0e-6;
    GONS_UINT max_iterations = 5;
    GONS_BOOL verbose = false;
    ApproxType HessianApproxType = ApproxType::Inverse;
  };

public:
  BFGSMethod(Function &f, X &x) : m_f(f), m_x(x) {}
  ~BFGSMethod() = default;
  void set_parameters(const BFGSMethodParameters &params) { m_params = params; }
  X get_x() const { return m_x; }
  double get_function_value() const { return m_f(m_x); }

private:
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> init_hessian_approximation(const Vector<T, Size> &x) {
    UNUSED(x); // in case warning
    return Identity<T, Size>();
  }

public:
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> DirectUpdateApprox(const Matrix<T, Size, Size> &B,
                                           const Vector<T, Size> &s,
                                           const Vector<T, Size> &y) {

    Matrix<T, Size, Size> ykykT = y.outerProduct(y.transpose());
    double inv_skTyk = 1.0 / (s.transpose().dot(y));
    Matrix<T, Size, Size> ykykT_skTyk = ykykT * inv_skTyk;
    Vector<T, Size> Bksk = B * s.transpose();
    Matrix<T, Size, Size> Bksk_BkskT = Bksk.outerProduct(Bksk.transpose());
    double skTBksk = s.transpose().dot(Bksk);
    Matrix<T, Size, Size> Bksk_BkskT_skTBksk = Bksk_BkskT * (1.0 / skTBksk);
    Matrix<T, Size, Size> result = B + ykykT_skTyk - Bksk_BkskT_skTBksk;
    return result;
  }
  template <typename T, GONS_UINT Size>
  Matrix<T, Size, Size> InverseUpdateApprox(const Matrix<T, Size, Size> &H,
                                            const Vector<T, Size> &s,
                                            const Vector<T, Size> &y) {

    Matrix<T, Size, Size> I = Identity<T, Size>();
    double rhok = (1.0 / s.transpose().dot(y));
    Matrix<T, Size, Size> ykskT = y.outerProduct(s.transpose());
    Matrix<T, Size, Size> I_rhoykskT = I - rhok * ykskT;

    Matrix<T, Size, Size> rhok_skskT = rhok * s.outerProduct(s.transpose());
    Matrix<T, Size, Size> result =
        I_rhoykskT.Transpose() * H * I_rhoykskT + rhok_skskT;
    return result;
  }

  BFGSMethodStatus Optimize() {
    GONS_UINT iter = 0;
    auto ApproxHessian = init_hessian_approximation(m_x);
    while (true) {
      X pk = m_params.HessianApproxType == ApproxType::Direct
                 ? -1 * ApproxHessian.Inverse() * m_f.gradient(m_x).transpose()
                 : -1 * ApproxHessian * m_f.gradient(m_x).transpose();
      X x_new = m_x + pk;
      X sk = x_new - m_x;
      X yk = m_f.gradient(x_new).transpose() - m_f.gradient(m_x).transpose();

      ApproxHessian = m_params.HessianApproxType == ApproxType::Direct
                          ? DirectUpdateApprox(ApproxHessian, sk, yk)
                          : InverseUpdateApprox(ApproxHessian, sk, yk);
      m_x = x_new;
      ++iter;
      if (m_params.tolerance) {
        LOG("Iteration: " << iter << " x: " << m_x << " f(x): " << m_f(m_x));
      }
      if (iter > m_params.max_iterations) {
        LOG_WARNING("Max iterations reached");
        return BFGSMethodStatus::MaxIterationReached;
      }
      if (m_f.gradient(m_x).Norm2() < m_params.tolerance) {
        LOG("BFGS, Converged");
        LOG("After " << iter << " iterations");
        LOG("x: " << m_x << " f(x): " << m_f(m_x));
        return BFGSMethodStatus::Success;
      }
    }

    return BFGSMethodStatus::Success;
  }

private:
  Function &m_f;
  X &m_x;
  BFGSMethodParameters m_params;
};
} // namespace qussinewton

} // namespace gons

#endif // GONS_QUASSINEWTON_H