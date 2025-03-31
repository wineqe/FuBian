import sympy as sp

# 定义符号变量
r, theta, U, a = sp.symbols('r theta U a', real=True, positive=True)
z = r * sp.exp(sp.I * theta)  # 极坐标下的复数表示

# 定义复势函数
Phi = U * (z + a**2 / z)

# 分离复势的实部和虚部
varphi, psi = sp.re(Phi), sp.im(Phi)

# 计算偏导数
dvarphi_dr = sp.diff(varphi, r)
dpsi_dtheta = sp.diff(psi, theta)
dpsi_dr = sp.diff(psi, r)
dvarphi_dtheta = sp.diff(varphi, theta)

# 根据柯西-黎曼方程计算残差
residual1 = dvarphi_dr - (1 / r) * dpsi_dtheta  # 柯西-黎曼方程的第一式
residual2 = dpsi_dr + (1 / r) * dvarphi_dtheta  # 柯西-黎曼方程的第二式

# 计算最大残差
max_residual = sp.Max(sp.Abs(residual1), sp.Abs(residual2))

# 数值评估
# 假设 U = 1, a = 1，并对 r 和 theta 赋予具体值
values = {U: 1, a: 1}

# 定义一个函数来计算最大残差的数值
def evaluate_residual(r_val, theta_val):
    # 输入：r_val(径向坐标), theta_val(角度坐标)
    # 输出：两个残差在给定点上的最大绝对值
    eval_residual1 = residual1.subs({**values, r: r_val, theta: theta_val}).evalf()
    eval_residual2 = residual2.subs({**values, r: r_val, theta: theta_val}).evalf()
    return max(abs(eval_residual1), abs(eval_residual2))

# 测试一些具体的 (r, theta) 值
test_values = [1, 2, 5]  # r 的测试值
theta_test = sp.pi / 4  # theta 的测试值（45度）

for r_val in test_values:
    max_residual_value = evaluate_residual(r_val, theta_test)
    # 将 max_residual_value 转换为浮点数后再格式化
    max_residual_value_float = float(max_residual_value)
    print(f"r = {r_val}, theta = {theta_test:.3f}, 最大残差 = {max_residual_value_float:.2e}")
    if max_residual_value_float < 1e-6:
        print("满足要求：最大残差小于 1e-6")
    else:
        print("不满足要求：最大残差大于或等于 1e-6")
