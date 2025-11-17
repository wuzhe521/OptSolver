#include "../core/quasi_newton.h"

using namespace gons::utilites::LOG_MSG;

using X = gons::Vector<double, 2>;
// 定义测试函数类，用于计算目标函数值和梯度
class TestFunction {
public:
  inline double operator()(const X &x) const {
    return x(0) * x(0) + (x(1) - 1) * (x(1) - 1);
  }
  inline X gradient(const X &x) { return {2 * x(0), 2 * (x(1) - 1)}; }
  inline gons::Matrix<double, 2, 2> hessian(const X &x) {
    return gons::Matrix<double, 2, 2>({{2, 0}, {0, 2}});
  }
};
int main() {
  TestFunction f;
  X x0 = {-100.0, 100.0};
  gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodParameters params;
  params.MethodType =
      gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodType::SR1;
  gons::quasi_newton::SRMethod<TestFunction, X> sr_method(f, x0);
  sr_method.set_parameters(params);
  sr_method.Optimize();
  // test SR2
  X x1 = {-100.0, -100.0};
  params.MethodType =
      gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodType::SR2;
  gons::quasi_newton::SRMethod<TestFunction, X> sr1_method(f, x1);
  sr1_method.set_parameters(params);
  sr1_method.Optimize();

  X x2 = {-100.0, -100.0};
  gons::quasi_newton::BFGSMethod<TestFunction, X> bfgs_method(f, x2);
  bfgs_method.Optimize();
  return 0;
}