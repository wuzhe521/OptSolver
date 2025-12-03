import numpy as np
import scipy.linalg
from typing import Tuple, Optional

def qp_primal_dual_interior_point(
    Q: np.ndarray, 
    c: np.ndarray, 
    A: Optional[np.ndarray] = None,
    b: Optional[np.ndarray] = None,
    A_eq: Optional[np.ndarray] = None,
    b_eq: Optional[np.ndarray] = None,
    x0: Optional[np.ndarray] = None,
    tol: float = 1e-8,
    max_iter: int = 100
) -> Tuple[np.ndarray, float, dict]:    # sourcery skip: low-code-quality
    """
    二次规划问题的原始-对偶内点法求解器
    
    求解问题:
        min  (1/2) x^T Q x + c^T x
        s.t. A x <= b      (不等式约束)
             A_eq x = b_eq  (等式约束)
             x >= 0
    
    参数:
        Q: 二次项矩阵 (n x n), 应该是对称半正定的
        c: 一次项系数 (n,)
        A: 不等式约束矩阵 (m_ineq x n)
        b: 不等式约束右端项 (m_ineq,)
        A_eq: 等式约束矩阵 (m_eq x n)  
        b_eq: 等式约束右端项 (m_eq,)
        x0: 初始点 (可选)
        tol: 收敛容忍度
        max_iter: 最大迭代次数
    
    返回:
        x_opt: 最优解
        f_opt: 最优目标函数值
        info: 求解信息字典
    """
    
    n = Q.shape[0]

    # 处理约束
    if A is None:
        A = np.zeros((0, n))
        b = np.zeros(0)
    if A_eq is None:
        A_eq = np.zeros((0, n))
        b_eq = np.zeros(0)

    m_ineq = A.shape[0]
    m_eq = A_eq.shape[0]

    # 初始化变量
    x = np.ones(n) if x0 is None else x0.copy()
    # 对偶变量和松弛变量
    lam = np.ones(m_ineq)  # 不等式约束的对偶变量
    s = np.ones(m_ineq)    # 不等式约束的松弛变量 (s = b - A x)
    nu = np.zeros(m_eq)    # 等式约束的对偶变量

    print(f"{'Iter':<4} {'Objective':<12} {'Duality Gap':<12} {'Mu':<12} {'Primal Feas':<12} {'Dual Feas':<12}")
    print("-" * 80)

    history = {
        'primal_obj': [],
        'duality_gap': [],
        'mu': [],
        'primal_feas': [],
        'dual_feas': []
    }

    for iter in range(max_iter):
        # 计算松弛变量
        s = b - A.dot(x)

        # 计算对偶间隙
        mu = np.dot(lam, s) / m_ineq if m_ineq > 0 else 0.0

        # 计算目标函数值
        f_val = 0.5 * x.dot(Q.dot(x)) + c.dot(x)

        # 计算残差
        # 对偶残差: Qx + c + A^T lam + A_eq^T nu
        r_dual = Q.dot(x) + c + A.T.dot(lam) + A_eq.T.dot(nu)

        # 原始残差
        r_pri_ineq = s - (b - A.dot(x)) if m_ineq > 0 else np.zeros(0)
        r_pri_eq = A_eq.dot(x) - b_eq if m_eq > 0 else np.zeros(0)

        # 互补松弛残差
        r_comp = lam * s if m_ineq > 0 else np.zeros(0)

        # 计算收敛性度量
        duality_gap = np.abs(mu * m_ineq) if m_ineq > 0 else 0.0
        primal_feas = (np.linalg.norm(r_pri_eq) + 
                      (np.linalg.norm(np.minimum(s, 0)) if m_ineq > 0 else 0.0))
        dual_feas = np.linalg.norm(r_dual)

        # 记录历史
        history['primal_obj'].append(f_val)
        history['duality_gap'].append(duality_gap)
        history['mu'].append(mu)
        history['primal_feas'].append(primal_feas)
        history['dual_feas'].append(dual_feas)

        print(f"{iter:<4} {f_val:<12.6f} {duality_gap:<12.2e} {mu:<12.2e} "
              f"{primal_feas:<12.2e} {dual_feas:<12.2e}")

        # 收敛检查
        if (primal_feas < tol and dual_feas < tol and 
            (m_ineq == 0 or duality_gap < tol)):
            print("\n收敛!")
            break

        # 构建牛顿系统
        sigma = 0.1  # 中心参数

        # 构造对角矩阵
        S = np.diag(s) if m_ineq > 0 else np.zeros((0, 0))
        L = np.diag(lam) if m_ineq > 0 else np.zeros((0, 0))
        S_inv = np.diag(1.0 / s) if m_ineq > 0 else np.zeros((0, 0))
        L_inv = np.diag(1.0 / lam) if m_ineq > 0 else np.zeros((0, 0))

        # 构造KKT系统矩阵
        # [ Q + A^T S^{-1} L A   A_eq^T ] [ dx ]   [ -r_dual + A^T S^{-1} (r_comp - sigma*mu*e) ]
        # [ A_eq                  0     ] [ dnu ] = [ -r_pri_eq                                ]

        M = Q + A.T.dot(S_inv).dot(L).dot(A) if m_ineq > 0 else Q.copy()
        # 构造KKT矩阵
        if m_eq > 0:
            KKT_top = np.hstack([M, A_eq.T])
            KKT_bottom = np.hstack([A_eq, np.zeros((m_eq, m_eq))])
            KKT_matrix = np.vstack([KKT_top, KKT_bottom])
        else:
            KKT_matrix = M

        # 构造右端项
        if m_ineq > 0:
            rhs_dual = -r_dual + A.T.dot(S_inv).dot(r_comp - sigma * mu * np.ones(m_ineq))
        else:
            rhs_dual = -r_dual

        if m_eq > 0:
            rhs_eq = -r_pri_eq
            rhs = np.concatenate([rhs_dual, rhs_eq])
        else:
            rhs = rhs_dual

        # 求解KKT系统
        try:
            if m_eq > 0:
                d_sol = scipy.linalg.solve(KKT_matrix, rhs, assume_a='sym')
                dx = d_sol[:n]
                dnu = d_sol[n:n+m_eq]
            else:
                dx = scipy.linalg.solve(KKT_matrix, rhs, assume_a='sym')
                dnu = np.zeros(m_eq)
        except Exception:
            # 如果直接求解失败，使用最小二乘
            if m_eq > 0:
                d_sol = np.linalg.lstsq(KKT_matrix, rhs, rcond=None)[0]
                dx = d_sol[:n]
                dnu = d_sol[n:n+m_eq]
            else:
                dx = np.linalg.lstsq(KKT_matrix, rhs, rcond=None)[0]
                dnu = np.zeros(m_eq)

        # 计算其他变量的步长
        if m_ineq > 0:
            ds = -A.dot(dx) - r_pri_ineq
            dlam = -S_inv.dot(L.dot(ds) + r_comp - sigma * mu * np.ones(m_ineq))
        else:
            ds = np.zeros(0)
            dlam = np.zeros(0)

        # 计算步长
        alpha_pri = 1.0
        alpha_dual = 1.0

        # 原始步长 (保持 x >= 0, s >= 0)
        for i in range(n):
            if dx[i] < 0:
                alpha_pri = min(alpha_pri, -0.99 * x[i] / dx[i])

        if m_ineq > 0:
            for i in range(m_ineq):
                if ds[i] < 0:
                    alpha_pri = min(alpha_pri, -0.99 * s[i] / ds[i])

            for i in range(m_ineq):
                if dlam[i] < 0:
                    alpha_dual = min(alpha_dual, -0.99 * lam[i] / dlam[i])

        # 应用安全因子
        eta = 0.95
        alpha_pri = min(1.0, eta * alpha_pri)
        alpha_dual = min(1.0, eta * alpha_dual)

        # 更新变量
        x += alpha_pri * dx
        if m_ineq > 0:
            s += alpha_pri * ds
            lam += alpha_dual * dlam
        if m_eq > 0:
            nu += alpha_dual * dnu

    else:
        print(f"\n达到最大迭代次数 {max_iter}")

    # 计算最终目标函数值
    f_opt = 0.5 * x.dot(Q.dot(x)) + c.dot(x)

    info = {
        'success': (primal_feas < tol and dual_feas < tol and 
                   (m_ineq == 0 or duality_gap < tol)),
        'message': '收敛' if iter < max_iter - 1 else '达到最大迭代次数',
        'iterations': iter + 1,
        'primal_feas': primal_feas,
        'dual_feas': dual_feas,
        'duality_gap': duality_gap,
        'history': history,
        'lam': lam if m_ineq > 0 else np.zeros(0),
        'nu': nu if m_eq > 0 else np.zeros(0)
    }

    return x, f_opt, info                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

# 测试例子
if __name__ == "__main__":
    print("=" * 60)
    print("测试1: 简单的二次规划问题")
    print("=" * 60)
    
    # 问题: min 0.5*(x1^2 + x2^2) + x1 + 2x2
    #       s.t. x1 + x2 >= 1
    #            x1, x2 >= 0
    
    # 转换为标准形式: min 0.5*x^T Q x + c^T x, s.t. A x <= b, x >= 0
    # 注意: x1 + x2 >= 1 转换为 -x1 - x2 <= -1
    
    Q = np.array([[1, 0],
                  [0, 1]])
    c = np.array([1, 2])
    A = np.array([[-1, -1]])  # 注意符号转换
    b = np.array([-1])
    
    print("问题: min 0.5*(x1^2 + x2^2) + x1 + 2x2")
    print("约束: x1 + x2 >= 1, x1,x2 >= 0")
    print("预期最优解在边界上\n")
    
    x_opt, f_opt, info = qp_primal_dual_interior_point(Q, c, A, b)
    
    print(f"\n最优解: x = [{x_opt[0]:.6f}, {x_opt[1]:.6f}]")
    print(f"最优目标值: {f_opt:.6f}")
    print(f"迭代次数: {info['iterations']}")
    print(f"对偶变量: λ = {info['lam']}")
    
    print("\n" + "=" * 60)
    print("测试2: 带等式约束的二次规划")
    print("=" * 60)
    
    # 问题: min 0.5*(x1^2 + x2^2) 
    #       s.t. x1 + x2 = 1
    #            x1, x2 >= 0
    
    Q = np.array([[1, 0],
                  [0, 1]])
    c = np.array([0, 0])
    A_eq = np.array([[1, 1]])
    b_eq = np.array([1])
    
    print("问题: min 0.5*(x1^2 + x2^2)")
    print("约束: x1 + x2 = 1, x1,x2 >= 0")
    print("预期最优解: x1=0.5, x2=0.5\n")
    
    x_opt, f_opt, info = qp_primal_dual_interior_point(Q, c, A_eq=A_eq, b_eq=b_eq)
    
    print(f"\n最优解: x = [{x_opt[0]:.6f}, {x_opt[1]:.6f}]")
    print(f"最优目标值: {f_opt:.6f}")
    print(f"等式约束对偶变量: ν = {info['nu']}")
    
    print("\n" + "=" * 60)
    print("测试3: 投资组合优化问题")
    print("=" * 60)
    
    # 简单的投资组合优化
    # min 0.5 * w^T Sigma w - mu^T w
    # s.t. sum(w) = 1, w >= 0
    
    np.random.seed(42)
    n_assets = 3
    Sigma = np.array([[0.1, 0.02, 0.01],
                      [0.02, 0.15, 0.03],
                      [0.01, 0.03, 0.08]])  # 协方差矩阵
    mu = np.array([0.08, 0.12, 0.10])      # 预期收益率
    
    A_eq = np.ones((1, n_assets))
    b_eq = np.array([1.0])
    
    print("投资组合优化问题")
    print("目标: 在风险最小化和收益最大化之间权衡")
    print("约束: 权重和为1, 不允许卖空\n")
    
    x_opt, f_opt, info = qp_primal_dual_interior_point(
        Q=Sigma, c=-mu, A_eq=A_eq, b_eq=b_eq  # 注意c=-mu因为我们要最大化收益
    )
    
    print(f"\n最优权重: {x_opt}")
    print(f"预期收益率: {mu.dot(x_opt):.4f}")
    print(f"风险(方差): {x_opt.dot(Sigma.dot(x_opt)):.4f}")
    print(f"夏普比率(假设无风险利率=0): {mu.dot(x_opt)/np.sqrt(x_opt.dot(Sigma.dot(x_opt))):.4f}")