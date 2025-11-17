# 罚函数法
约束优化问题
一般约束优化问题的数学表达
$$
min f(x)\\s.t.\\
g_{l}(x)\le 0\\
h_{i}(x)=0\\
i=1\ 2\ 3\ 4...n\\
l=1\ 2\ 3\ 4...m$$
我们求解此类问题的一般思想是，将约束优化问题转化为无约束优化问题。
罚函数法
1. 等式约束问题
$$min f(x)\\s.t.\\
h_{i}(x)=0\\i=1\ 2\ 3\ 4...n$$
构造惩罚项$$p(x)$$
惩罚项要满足以下两点：
  - 当x在可行域S中
$$
x \in S=>h_{i}(x) = 0\\p(x) = 0\\
$$
  - 当x不在可行域S中
$$
x \not\in S=>h_{i}(x) \neq 0\\p(x) > 0\\$$
设计惩罚项
$$p(x) = \sum_{i=1}^l{h_{i}^2(x)}$$
根据原约束优化问题，构造新的罚函数 
$$P(x,\sigma) = f(x) + \sigma p(x)\\
\sigma > 0$$
同时产生一个无约束优化问题：
$$\min\limits_{x\in R^N} P(x,\sigma) = f(x) + \sigma p(x)
\\
\sigma > 0$$
原问题的另一种表达，在可行域S内，求f(x)关于x的最小值：
$$\min\limits_{x \in S} f(x)\\
s:\{x|h_{i}(x) = 0\ \ \ i=1\ 2\ 3...l\}$$
因为在可行域内，惩罚项为0，因此也就等价于
$$\min\limits_{x \in S} f(x) + \sigma p(x)\\
s:\{x|h_{i}(x) = 0\ \ \ i=1\ 2\ 3...l\}\\
\sigma > 0$$
又因为：
$$\min\limits_{x \in S} f(x) + \sigma p(x)\geq  \min\limits_{x \in R^N} f(x) + \sigma p(x) = P(x, \sigma)\\
\sigma > 0$$
即产生的无约束优化问题，是原约束优化问题的一个下届
显然，最小值的取值和 $$\sigma$$的取值有关
记 $$\theta(\sigma) = \min\limits_{x\in R^N} f(x) + \sigma p(x) \ 
\sigma > 0$$, $$\theta(\sigma) $$为原问题最优解的一个下届。在对偶分析中，已经了解，下届越大越好
求$$\max\limits_{\sigma} \theta(\sigma) \ 
\sigma > 0$$，因为 $$p(x) \geq 0$$ $$f(x) + \sigma p(x)$$是关于 $\sigma$单调递增的, min/max不影响单调性。
$\max\limits_{\sigma} \theta(\sigma) $ 等价于 $$ \lim_\{\sigma \to \inf \} \theta (\sigma) $$

例子：
$$min\ x_{1} + x_{2} \\
x_{2} - x_{1}^2 = 0$$
按照原来的方式，通过图形法或KKT点来求解，在此，运用罚函数法：
解：
构造罚函数：
$$P(x, \sigma) = x_{1} + x_{2} + \sigma(x_{2} - x_{1}^2)^2$$
原问题转化为求 $$\min\limits_{x}P(x, \sigma)$$，利用一阶条件- $$\frac{ \partial P }{ \partial x } = 0$$
$$\frac{ \partial P }{\partial x_{1} } = 1 + 2\sigma(x_{2} - x_{1}^2)*(-2x_{1}) = 0\\
\frac{ \partial P }{\partial x_{2} } =1 + 2\sigma(x_{2} - x_{1}^2) = 0$$
解得：
$$x_{1}(\sigma) = -\frac{1}{2}
\\
x_{2}(\sigma) =  \frac{1}{4} - \frac{1}{2\sigma}$$ $$Hessian = \begin{matrix} \frac{ \partial P }{ \partial x_{1} \partial x_{1}} & \frac{ \partial P }{ \partial x_{1} \partial x_{2}}\\ \frac{ \partial P }{ \partial x_{1} \partial x_{2}}&\frac{ \partial P }{ \partial x_{2} \partial x_{2}}\end{matrix} > 0=> minimal\ value$$
令 $$\sigma \rightarrow  inf$$
$$x^* = x(\sigma \rightarrow inf)= \begin{bmatrix} x_{1}\\ x_{2} \end{bmatrix} = \begin{bmatrix} -\frac{1}{2}\\ \frac{1}{4}\end{bmatrix}$$
通用求解过程：
1. 取初始点 $x_{0}$ 初始惩罚系数 $\sigma_{1} > 0$， 初始容差 $\varepsilon >0$， 初始迭代步 $k = 1$
2. 以 $x_{k-1}$为起始点，求解
  无约束优化问题 $\min f(x) + \sigma_{k} p(x)$
  得到 $x_{k}$
3. 判断 $\sigma_{k}p(x_{k}) < \varepsilon$，若成立，终止计算
4. 更新 $\sigma_{k+1}$, 满足条件： $\sigma_{k+1} > \sigma_{k}$
      - 线性递增 $\sigma_{k +1} = \beta\sigma_{k} \    \ \ \   \beta>1$
  -  根据 $\min P(x, \sigma)$求解的难易程度，1. 容易 $\beta = 10$, 2 难 $\beta =2$
5. 跳转至step 2

性质1：
假设 $$x_{k} $$是 $$\min f(x) + \sigma_{k} p(x)$$的全局最优解，在求解过程中，得到三个点列：
1. $$P_{k}(x_{k} + \sigma_{k}p(x_{k}))$$ 
    单调递增
    $$P_{k}(x_{k} + \sigma_{k}p(x_{k})) = \min P(x, \sigma_{k}) = \theta(\sigma_{k})$$
    因为 $\theta(\sigma)$单调递增，且 $$\sigma_{k+1} > \sigma_{k}$$,
    所以， $$\theta(\sigma_{k}) = P_{k}(x_{k} + \sigma_{k}p(x_{k}))$$单调递增
2. $$\{p(x_{k})\}$$单调递减
       $$x_{k} $$是 $$\min f(x) + \sigma_{k} p(x)$$的全局最优解
       $$x_{k+1} $$是 $$\min f(x) + \sigma_{k+1} p(x)$$的全局最优解
得到：
$$ f(x_{k}) + \sigma_{k} p(x_{k})\leq f(x_{k+1}) + \sigma_{k} p(x_{k+1})$$
$$f(x_{k+1}) + \sigma_{k+1} p(x_{k+1})\leq f(x_{k}) + \sigma_{k+1} p(x_{k})$$
化简一下：
$$\sigma_{k} p(x_{k}) + \sigma_{k+1} p(x_{k+1}) \leq \sigma_{k} p(x_{k+1}) + \sigma_{k+1} p(x_{k})$$
$$\sigma_{k} p(x_{k}) -\sigma_{k} p(x_{k+1})  \leq  \sigma_{k+1} p(x_{k}) - \sigma_{k+1} p(x_{k+1})$$
$$\sigma_{k}[ p(x_{k}) - p(x_{k+1})]  \leq  \sigma_{k+1}[ p(x_{k}) - p(x_{k+1})]$$
已知， $$\sigma_{k+1} > \sigma_{k}$$, 要使得等式成立：
$$ p(x_{k}) - p(x_{k+1}) > 0$$
即 $$\{p(x_{k})\}$$单调递减
3. $$\{f(x_{k})\}$$单调递增
根据2中的不等式 1：
$$ f(x_{k}) + \sigma_{k} p(x_{k})\leq f(x_{k+1}) + \sigma_{k} p(x_{k+1})$$
转化为
$$ f(x_{k}) - f(x_{k+1}) \leq   \sigma_{k}[  p(x_{k+1}) - p(x_{k})]$$
因为：
  $$\{p(x_{k})\}$$单调递减， $$ \ \ \ \ p(x_{k+1}) - p(x_{k}) < 0$$
  得到
  $$ f(x_{k}) - f(x_{k+1}) \leq   \sigma_{k}[  p(x_{k+1}) - p(x_{k})] < 0$$
  $$\{f(x_{k})\}$$单调递增

性质2：
假设 $$x_{k} $$是 $$\min f(x) + \sigma_{k} p(x)$$的全局最优解， $$\sigma_{k} \rightarrow inf$$ , $$\{x_{k}\}$$的任一聚点是 $$\min\limits_{x \in S} f(x) + \sigma p(x)$$的全局最优解
可行且相等，证明：
假设$$\min\limits_{x \in S} f(x) + \sigma p(x)$$的全局最优解为 $$x^*$$
因为，$$x_{k} $$是 $$\min f(x_{k}) + \sigma_{k} p(x_{k})$$的全局最优解, 且$$x^*$$一定在可行域内 ，$$p(x^*) = 0$$ ：
$$P_{k}(x_{k} ， \sigma_{k}) = f(x_{k}) + \sigma_{k}p(x_{k})\leq f(x^*) + \sigma p(x^*)  = f(x^*)$$
$$f(x_{k}) \leq f(x_{k}) + \sigma_{k}p(x_{k})\leq f(x^*) + \sigma p(x^*)  = f(x^*)$$
得：
$$P_{k}(x_{k} ， \sigma_{k})\leq f(x^*)$$ 点列 $$P_{k}$$单调上升有上界
$$f(x_{k}) \leq f(x^*)$$点列 $$f(x_k)$$单调上升有上界
有 $$P_{k}(x_{k} ， \sigma_{k}) = f(x_{k}) + \sigma_{k}p(x_{k})$$
$$p(x_{k}) = \frac{1}{\sigma_{k}}[P_{k}(x_{k} ， \sigma_{k})- f(x_{k})] $$ 
当 $$\sigma_{k} \rightarrow inf$$ $$p(x_{k}) \rightarrow 0$$
即点列 $$\{x_{k}\}$$ 的聚点 $$\bar{x}$$，满足 约束条件 ，是可行点。
因为 上面推导 ：$$f(x_{k}) \leq f(x^*)$$点列 $$f(x_k)$$单调上升有上界
$$f(\bar{x}) \leq f(x^*)$$
$$\min\limits_{x \in S} f(x) + \sigma p(x)$$的全局最优解为 $$x^*$$
$$f(x^*)  \leq f(\bar{x})$$
得到：
$$f(x^*)  =  f(\bar{x})$$

性质3 设 $$x_{k}$$满足 $$\bigtriangledown_{x} P(x_{k}, \sigma_{k})=0$$，点列 $$\{ x_{k}\} $$的聚点 $$\bar{x}$$，假设$$\bar{x}$$处， $$\bigtriangledown_{x}h_{i}(x)$$线性无关，则$$\bar{x}$$是原约束优化问题的KKT点
证明：
$$\min\limits_{x\in R^N} P(x_{k},\sigma_{k}) = f(x) + \sigma_{k} p(x) = f(x) + \sigma_{k} \sum_{i=1}^l{h_{i}^2(x)}\\
\bigtriangledown_{x} P(x_{k}, \sigma_{k}) = \bigtriangledown f(x_{k}) + 2\sigma_{k} \sum_{i=1}^l{h_{i}(x_{k})}\bigtriangledown h_{i}(x_{k})\\

\sum_{i=1}^l{h_{i}(x_{k})}\bigtriangledown h_{i}(x_{k}) = -\frac{\bigtriangledown f(x_{k}) }{2\sigma_{k}}$$
当 $$k \rightarrow \infty $$上式 右侧 趋向于 0 ， $$\sum_{i=1}^l{h_{i}(x_{k})}\bigtriangledown h_{i}(x_{k}) = 0$$ 得到 $$\sum_{i=1}^l{h_{i}(x_{k})}=0$$，即  $$\bar{x}$$是可行点
记 $$\lambda^*_{i} = 2\sigma_{k}h_{i}(x_{k})$$，得到门面方程组：
$$\bigtriangledown f(x_{k}) +\sum_{i=1}^l{\lambda^*_{i}\bigtriangledown h_{i}(x_{k}) }= 0$$
又因为， $$\bigtriangledown_{x}h_{i}(x)$$线性无关， $$\lambda^*_{i}$$的解唯一
当 $$k \rightarrow \infty $$上式，可以改写为 
$$\bigtriangledown f(\bar{x}) +\sum_{i=1}^l{\lambda_{i}\bigtriangledown h_{i}(\bar{x}) }= 0$$ 且 $$\sum_{i=1}^l{h_{i}(\bar{x})}=0$$  
KKT点

2. 不等式约束 - 一般形势的约束优化问题
$$min f(x)\\s.t.\\
g_{l}(x)\le 0\\
h_{i}(x)=0\\
i=1\ 2\ 3\ 4...n\\
l=1\ 2\ 3\ 4...m$$
构造一般形式的惩罚函数：
$$P(x) = f(x) + \sum_{i=1}^m{\max|0,g_{i}(x)|^\alpha} + \sum_{l=1}^n{|h_{l}(x)|^\beta}\ \ \ \ \ \alpha \ \beta\ge1$$
罚函数的缺点
罚函数法的思想是，通过不断地增加 $$\sigma_{k}$$，解$$\min f(x) + \sigma_{k} p(x)$$的全局最优解，使得点列 $$x_{k}$$逐步逼近最优解，从原理可知， $$\sigma_{k }\rightarrow \infty $$，点列的聚点 $$\bar{x}$$就是约束问题的全局最优解。
随着$$\sigma_{k }\rightarrow \infty $$ 惩罚函数法越病态，这使得问题随着惩罚参数的增加而变得更加难以求解，从而导致数值误差和收敛速度缓慢，并且解通常不可行，需要进一步的步骤来修正。此外，惩罚权重的微调也很困难，并且该方法会增加搜索空间的粗糙度，可能导致算法寻求局部最小值而不是全局最优值。
  - 随着惩罚参数的增加以更严格地执行约束，惩罚函数变得更加病态，这意味着它具有较大的梯度和突然的变化。
  - 这种病态会导致数值误差，并导致用于解决惩罚问题的无约束优化方法收敛缓慢。
  - 惩罚问题的解一般不是原约束问题的精确解，并且随着惩罚参数趋于无穷大，可能无法得到真正的解。
  - 由于该方法对可行域的偏差进行惩罚，因此被惩罚问题的解往往略微不可行。
  - 罚参数太低，会导致约束不充分，罚参数太大，会导致问题病态。
  - 罚函数法是对原始问题的近似，会导致“解析误差”
  - 罚函数法的解可能存在不可行域内，导致有些问题无解


  # 等式约束的二次罚函数法
  目标函数是二次函数，约束是线性约束的数学形式:
  $$ \min \frac{1}{2}x^TQx + c^Tx  $$
  线性不等式约束的数学形式:
  $$ Ax = b \qquad i = 1,2,\dots,m $$
    
  线性约束的罚函数形式:
  $$ \frac{1}{2}x^TQx + c^Tx + \frac{\rho}{2} (Ax - b)^2 $$
  
  罚函数的导数：
  $$ \frac{1}{2}Qx + c + \rho (Ax - b)A^T  $$

  二次规划，等式约束的罚函数法，代码实现：
  ```cpp
  namespace penalty_function {
using namespace gons::utilites::LOG_MSG;
using namespace gons::gradientsearch;

template <typename DataType, GONS_UINT Var_Size, GONS_UINT Con_Size>
class function {
private:
  GONS_FLOAT rho_ = 10.0;
  Matrix<DataType, Var_Size, Var_Size> P_;
  Vector<DataType, Var_Size> q_;
  Matrix<DataType, Con_Size, Var_Size> A_;
  Vector<DataType, Con_Size> b_;

public:
  function(Matrix<DataType, Var_Size, Var_Size> &P,
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
    return x.dot(P_ *x.transpose()) + q_.dot(x);
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
class QuadPenaltyFunction {
  using X = Vector<DataType, Var_Size>;
  using F = function<DataType, Var_Size, Con_Size>;

public:
  struct quad_penalty_parmeters {
    DataType rho = 1;
    DataType alpha = 2; // penalty growth rate
    GONS_UINT max_iter = 1000;
    GONS_BOOL verbose = false;
    GONS_BOOL war_start = true;
    DataType tol = 1e-5;
  };
  enum class PenaltyOptStatus { Failure, MaxIterationReached, Success };

public:
  // for quadratic function:
  // f(x) = 0.5 * x' * P * x + q' * x
  // for linear inequality constraint:
  // l <= g(x) = A * x + b <= u

  QuadPenaltyFunction(Matrix<DataType, Var_Size, Var_Size> &P,
                      Vector<DataType, Var_Size> &q,
                      Matrix<DataType, Con_Size, Var_Size> &A,
                      Vector<DataType, Con_Size> &b)
      : P_(P), q_(q), A_(A), b_(b), f_(P, q, A, b) {}
  Vector<DataType, Var_Size> get_result() const {
    return x_;
  }
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
  function<DataType, Var_Size, Con_Size> f_;
  // unconstrained search method
};

} // namespace penalty_function
```

测试代码：
```cpp
int main() {

  Matrix<double, 2, 2> P = {{4.0, 1.0}, {1.0, 2.0}};
  Vector<double, 2> q = {1.0, 1.0};
  Matrix<double, 1, 2> A = {{1.0, 1.0}};
  Vector<double, 1> b = {1.0};
  penalty_function::QuadPenaltyFunction<double, 2, 1> PenaltyFunc(P, q, A, b);
  PenaltyFunc.Optimize();
  X rlt = PenaltyFunc.get_result();
  std::cout << "result: " << rlt << std::endl;
  
  return 0;
}
```  
输出结果：
```bash
result: 0.249999 0.749996 
```


# 不等式约束的二次罚函数法
目标函数是二次函数，约束是线性约束的数学形式:
$$ \min \frac{1}{2}x^TQx + c^Tx  $$
线性不等式约束的数学形式:
$$ l \le Ax \le u \qquad i = 1,2,\dots,m $$
其中，等式约束是以不等式约束的形式表示的：  
$$ Ax = b \rightarrow  b \le Ax \le b $$

算法的设计：  
1. 设置初始罚系数 $\delta_0 > 0$ , 求解精度 $tol$ , 初始解$ x_0$ 罚系数增长率 $\alpha > 1$
2. 迭代：  
2.1 够构造二次罚函数： $  \min \frac{1}{2}x^TQx + c^Tx + \delta_k \sum_{i=1}^m \max(0, l_i - Ax_k)_i^2 + \delta_k \sum_{i=1}^m \max(0, Ax_k - u_i)_i^2 $  
2.2 用B优化二次罚函数，得到解$x_k$
3. 满足停止条件 - 罚函数的梯度2范数 < $tol$，或者迭代次数超过最大限制，则结束循环
4. 继续增大罚系数 $\delta_{k+1} = \alpha \delta_k$
5. 迭代步骤 2

代码如下：
```cpp
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
```
测试函数采用OSQP官网给出的例子：
$$ \begin{aligned}
&x\frac{1}{2} '\begin{bmatrix} 4& 1 \\ 1& 2\end{bmatrix}x + \begin{bmatrix}1 \\ 1\end{bmatrix}x \\
&s.t. \\
&\begin{bmatrix}1\\0\\0\end{bmatrix} \le \begin{bmatrix}1 & 1\\1 & 0\\0& 1\end{bmatrix}x \le \begin{bmatrix}1 \\ 0.7 \\ 0.7 \end{bmatrix}
\end{aligned}$$

测试代码：
```cpp
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
  ```
  输出：
  > result: 0.299622 0.700024  
  
  osqp的测试结果为：
  > 0.29877108 0.701228


  <strong>结果与OSQP解相差大概在 1e-3 这个数量级上。尝试继续增大罚系数来提高精度，但是结果没有得到改善，反而引起了内部梯度方法失败， 尝试了缩小BB步长增长速度，限制BB步长大小，加入BB步长的飞单调性约束，但是结果没有得到改善。</strong>   
  请教了知乎的用户 [@凡人一品](https://www.zhihu.com/question/1972328690827388648) 他给出了很多建议，同时也建议对于二次规划问题，最好使用内点法或者ADMM。尽量避开罚函数因罚系数增大而导致的“病态”问题。  
  
  
