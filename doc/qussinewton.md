# Quasi-Newton Methods
牛顿法在理论上和实践中均取得很好的效果．然而对于大规模问题，函数的海瑟矩阵计算代价特别大或者难以得到，即便得到海瑟矩阵我们还需要求解一个大规模线性方程组．那么能否使用海瑟矩阵或其逆矩阵的近似来进行牛顿迭代呢？拟牛顿法便是这样的算法，它能够在每一步以较小的计算代价生成近似矩阵，并且使用近似矩阵代替海瑟矩阵而产生的迭代序列仍具有超线性收敛的性质．  
拟牛顿方法不计算海瑟矩阵 $\nabla^2 f(x)$，而是构造其近似矩阵 $B^k$ 或其逆的近似矩阵 $H^k$ ．我们希望 $B^k$或$H^k$仍然保留海瑟矩阵的部分性质，例如使得 $d^k$ 仍然为下降方向．

在牛顿法中已知，在第k 步迭代中， $x^k$ 处， 梯度函数 $\nabla f(x)$ 满足：
$$f(x^k + d^k) \approx f(x^k) + \bigtriangledown f(x^k)(d^k)+ \frac{1}{2}(d^k)^T \bigtriangledown^2 f(x^k)(d^k)$$
步长 $d^k$ 最佳的选取是：$d^k = -\bigtriangledown^2 f(x^k)^{-1} \bigtriangledown f(x^k)$

拟牛顿法的理论基础就是找到一个矩阵 $B^k$，使得 $B^k d^k = - \bigtriangledown f(x^k)$

# 割线法



连续二阶可微函数 $f(x)$ 的 梯度函数 $\nabla f(x)$， 在 $x^{k+1}$ 处， 泰勒展式为：
$$
\nabla f(x) = \nabla f(x^{k+1})^T  + \nabla^2 f(x^{k+1}) (x - x^{k+1}) + O(||x - x^{k+1}||^2)
$$
转换一下：
$$
\nabla^2 f(x^{k+1}) (x - x^{k+1}) + O(||x - x^{k+1}||^2) = \nabla f(x) - \nabla f(x^{k+1})^T (x - x^k)
$$
令 $x = x^k$ ， $s^k = x^{k+1} - x^k$，$ y^k = \nabla f(x^{k+1}) - \nabla f(x^k)^T$，则可以得到：
$$
\nabla^2 f(x^{k+1}) s^k + O(||s^k||^2) = y^kx
$$

忽略高阶项：
$$
\nabla^2 f(x^{k+1}) s^k \approx y^kx
$$

我们构造的$f(x)$的 Hessian矩阵的近似矩阵$B$, 和近似逆矩阵$H$，仍满足：  
$$
B^k s^k = y^kx \\

s^k = H^k y^k  \tag{1}
$$

方程（1） 也被成为割线方程。


另外一个角度看割线方程：  
牛顿法的本质是对 $f(x)$在 点$x^k$处二阶泰勒展开， 然后根据展开后的公式，固定 $f(x)$ 和 $\nabla f(x)$， 求解 $d^{k}$ ， 求使得 $f(x+ d^{k})$ 最小化的 $d^{k}$  
当我们用割线法，用 $B^k$近似代替 $ \nabla^2 f(x^{k+1})$ 后， 将整个步长函数$M(d)$ 在点$k+1$ 的二阶泰勒近似：
$$
M_{k+1}(d) = f(x^{k+1}) + \nabla f(x^{k+1})^T d^k + \frac{1}{2}(d^k)^T B^k d^k
$$
$$
\nabla M_{k+1}(d) = \nabla f(x^{k+1}) + B^k d^k
$$
我们要求 在 $M_{k+1}(d)$ 在 $d = -s^k$和 $d = 0$ 处的梯度，和 $f(x)$ 在 $x=x^k$和$x = x^{k+1}$处的梯度保持一致。
$$
\nabla M_{k+1}(-s^k) = M_{k+1}(-s^k) =  \nabla f(x^{k+1}) + B^k s^k  = \nabla f(x^k) \\
\nabla M_{k+1}(0) = \nabla f(x^k) 
$$
可以看出，$\nabla M_{k+1}(0)$ 必然满足，主要 $\nabla f(x^{k+1}) + B^k s^k  = \nabla f(x^k) $成立，即可满足要求。
可得：
$$
B^k s^k  = \nabla f(x^k) - \nabla f(x^{k+1})  = y^k
$$

## 算法框架
1. 给定初始点 $x_0$， $B^0$ 或者 $H^0$
2. 循环迭代：
3. 计算方向 $ d^k = -(B^k)^{-1} \nabla f(x^k)$ 或者 $d^k = -H^k \nabla f(x^k)$
4. 线搜索确定步长 $\alpha^k$
5. 更新 $x^{k+1} = x^k + \alpha^k d^k$
6. 更新 $B^{k+1}$ 或者 $H^{k+1}$
7. 判断是否收敛，如果 $||d^k|| < \epsilon$, 则停止迭代，否则返回 2


在实际应用中基于$H^{k}$的拟牛顿法因其无需求解线性方程组，更加实用。基于 $B^{k}$ 的方法稳定性更好，但计算量更大。

那么，剩下的问题就成了如何选取 $B^{k+1}$ 或者 $H^{k+1}$ ？

## 拟牛顿矩阵更新形式

### [秩一更新](https://zhuanlan.zhihu.com/p/144736223)

这个方法继承的思路是：我希望在满足割线方程的条件下，每一步对于$B^k$的更新能够尽可能的小。这个“尽可能小”得含义，其实就是只修改矩阵的一个特征值：
$$
\begin{aligned}
&B^{k+1} = \min ||B- B^k||^2 \\
&s.t.\\
&B^{k+1} s^k = y^k
\end{aligned}
$$
有一种矩阵的变换方式叫做秩一更新，$ A + auu^T $, 在正定矩阵$A的基础上，加上一个秩为1的矩阵。
$$
B^{k+1} = B^k + auu^T \\
$$
其中, $auu^T$ 是一个正定，秩为1的矩阵， $a$ 是一个正数， $u$ 是一个单位向量。$a$和$u$是待定矩阵。
这种方法的好处是，可以套用 [Sherman-Morrison公式](https://en.wikipedia.org/wiki/Sherman%E2%80%93Morrison_formula)来求出 $B^{k+1}$ 的逆矩阵。将求逆 O(n^3) 的复杂度，降低为 O(n^2)。但前提是，海瑟矩阵 $B^k$ 是正定的 ，且满足割线条件。

根据割线方程：
$$
B^{k+1}s^k = B^ks^k + \alpha u  u^T s^k = -y^k \tag{1}
$$
$$
\alpha u^T  u  s^k = y^k - B^ks^k \tag{2}
$$
令 $u = y^k - B^ks^k$ ，则 (2) 可以写成：
$$
\alpha (y^k - B^ks^k)^Ts^k(y^k - B^ks^k) = y^k - B^ks^k \tag{3}
$$
得到， $\alpha$的更新格式:
$$
\alpha = \frac{(y^k - B^ks^k)(y^k - B^ks^k)^T}{(y^k - B^ks^k)^Ts^k}
$$
同理，可以得到基于$H^k$的更新格式：
$$
H^{k+1} = H^k + \frac{(s^k - H^ky^k)(s^k - H^ky^k)^T}{(s^k - H^ky^k)^T y^k}  \tag{4}
$$

SR1 公式虽然结构简单，但是有一个重大缺陷：它不能保证矩阵在迭代过程中保持正定, 因此不能保证矩阵的逆矩阵存在．因此在实际中较少使用 SR1 公式。

SR 方法的代码实现：
```cpp
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

    Matrix QuassiH =
        m_f.hessian(m_x).Identity(); // Initial Hessian approximation

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
};
```

测试函数 using $f(x) = x^2 + (y -1)^2 $,  代码：
```cpp
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
```
主函数：
```cpp
TestFunction f;
X x0 = {-100.0, 100.0};
gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodParameters params;
params.MethodType = gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodType::SR1;
gons::quasi_newton::SRMethod<TestFunction, X> sr_method(f, x0);
sr_method.set_parameters(params);
sr_method.Optimize();
// test SR2
X x1 = {-100.0, -100.0};
params.MethodType = gons::quasi_newton::SRMethod<TestFunction, X>::SRMethodType::SR2;
gons::quasi_newton::SRMethod<TestFunction, X> sr1_method(f, x1);
sr1_method.set_parameters(params);
sr1_method.Optimize();
```
输出结果：
```bash
SRMethod: Optimization using SR1 finished.
SRMethod: After Iterations: 3
SRMethod: Final X : 0 1 

SRMethod: Function Function value: 0

SRMethod: Optimization using SR2 finished.
SRMethod: After Iterations: 3
SRMethod: Final X : 0 1 

SRMethod: Function Function value: 8.07794e-28
```

## BFGS Method 

SR1 公式虽然结构简单，但是有一个重大缺陷：它不能保证矩阵在迭代过程中保持正定．为了客服这个缺陷，我们可以使用 BFGS 方法，BFGS 方法的矩阵是正定的，并且可以保证矩阵在迭代过程中保持正定．
BFGS 采用秩二更新的方法， Rank-2 update
$$
B^{k+1} = B^k + \alpha_k uu^T + \beta_k vv^T
$$
其中： $u^k$， $v^k$， $\alpha_k$， $\beta_k$ 是迭代过程中的待定系数。
根据割线方程：
$$
\begin{aligned}
B^{k+1} * s^k &= y^k \\
[B^k + \alpha_k uu^T + \beta_k vv^T]s^k &= y^k \\
B^k * s^k + \alpha_k uu^T * s^k + \beta_k vv^T * s^k &= y^k \\ 
u^T * s^k\quad and \quad v^T * s^k \quad is  \quad constant:\\
(\alpha_k u^T * s^k) * u + (\beta_k v^T * s^k) * v &= y^k - B^k * s^k 
\end{aligned}
$$
$u^k$， $v^k$， $\alpha_k$， $\beta_k$ 的取法非常多，暂令：
$$
u = y^k, \qquad \alpha_k u^T * s^k = 1\\
v = B^ks^k, \qquad \beta_k v^T * s^k = -1
$$
得到更新方程：
$$
\begin{aligned}
 B^{k+1} &= B^k + \alpha_k uu^T + \beta_k  vv^T \\
  &= B^k + \frac{y^k(y^k)^T}{(y^k)^Ts^k}  - \frac{B^ks^k(B^k*s^k)^T}{(B^ks^k)^Ts^k}
 \end{aligned}
$$
如果采用一种形式的割线方程，我们则可得到 更简单的更新方程：

$$
\begin{aligned}
H^{k+1} = (I - \rho^k y^k (s^k)^T) ^T H^k (I - \rho^k s^k (y^k)^T) + \rho^k s^k (s^k)^T
 \end{aligned}
 $$
其中 $\rho^k = \frac{1}{(s^k)^Ty^k}$ 

代码实现:
```cpp
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
      ApproxHessian.Print("ApproxHessian: \n");
      X pk = m_params.HessianApproxType == ApproxType::Direct
                 ? -1 * ApproxHessian.Inverse() * m_f.gradient(m_x).transpose()
                 : -1 * ApproxHessian * m_f.gradient(m_x).transpose();
      X x_new = m_x + pk;
      LOG("pk: \n" << pk);
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
```
测试函数仍采用，SR-秩1更新的测试函数。 $ f(x, y) = x^2 + (y - 1)^2 $.
测试代码如下：
```c++
  X x2 = {-100.0, -100.0};
  gons::quasi_newton::BFGSMethod<TestFunction, X> bfgs_method(f, x2);
  bfgs_method.Optimize();
```
测试结果如下：
```bash
BFGS, Converged
After 2 iterations
x: -1.42109e-14 1 
f(x): 1.00974e-27
```

## L-BFGS 
拟牛顿法虽然克服了计算海瑟矩阵的困难，但是它仍然无法应用在大规模优化问题上．一般来说，拟牛顿矩阵 $B^k$ 或$ H^k $是稠密矩阵，而存储稠密矩阵要消耗 $O(n^2)$ 的内存，这对于大规模问题显然是不可能实现的．

<font color="red"> 暂不实现 </font>