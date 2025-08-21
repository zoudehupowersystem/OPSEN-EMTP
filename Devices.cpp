#include "Devices.h"
#include "Grid.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

constexpr double PI = 3.14159265358979323846;

// --- Resistor ---
Resistor::Resistor(int n1, int n2, double r)
    : n1(n1)
    , n2(n2)
{
    if (std::abs(r) < 1e-9)
        throw std::runtime_error("Resistor value cannot be zero.");
    G = 1.0 / r;
}

void Resistor::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 电阻的印章法：
    // 对角线元素 G_ii 和 G_jj 增加 G
    // 非对角线元素 G_ij 和 G_ji 减去 G
    if (n1 > 0)
        triplets.emplace_back(n1 - 1, n1 - 1, G);
    if (n2 > 0)
        triplets.emplace_back(n2 - 1, n2 - 1, G);
    if (n1 > 0 && n2 > 0) {
        triplets.emplace_back(n1 - 1, n2 - 1, -G);
        triplets.emplace_back(n2 - 1, n1 - 1, -G);
    }
}

void Resistor::updateState(const Eigen::VectorXd& V, double dt)
{
    // 根据节点电压计算当前流过电阻的电流 I = (V1 - V2) * G
    double v1 = (n1 > 0) ? V(n1 - 1) : 0.0;
    double v2 = (n2 > 0) ? V(n2 - 1) : 0.0;
    double v_new = v1 - v2;
    i_now = G * v_new;
}

// --- series_RL ---
series_RL::series_RL(int n1, int n2, double r, double l)
    : n1(n1)
    , n2(n2)
    , R(r)
    , L(l)
    , i_hist(0.0)
    , v_hist(0.0)
{
    if (std::abs(r) < 1e-12 || std::abs(l) < 1e-12)
        throw std::runtime_error("R or L value in series_RL cannot be zero.");
}

void series_RL::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 梯形法则离散化: i(t) = G_eq * v(t) + I_hist
    // 等效电导 G_eq = 1 / (R + 2L/dt)
    G_eq = 1.0 / (R + (2.0 * L) / dt);
    // 历史电流更新系数
    value1_RL = (1 - dt * R / 2.0 / L) / (1 + dt * R / 2.0 / L);

    if (n1 > 0)
        triplets.emplace_back(n1 - 1, n1 - 1, G_eq);
    if (n2 > 0)
        triplets.emplace_back(n2 - 1, n2 - 1, G_eq);
    if (n1 > 0 && n2 > 0) {
        triplets.emplace_back(n1 - 1, n2 - 1, -G_eq);
        triplets.emplace_back(n2 - 1, n1 - 1, -G_eq);
    }
}

void series_RL::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // 更新历史电流源 I_hist(t) = f(i(t-dt), v(t-dt))
    // 此处的递推公式是梯形法则的一种简化形式
    i_hist = (i_hist + v_hist * G_eq) * value1_RL + v_hist * G_eq;
    // 将历史电流注入到节点
    if (n1 > 0)
        I(n1 - 1) -= i_hist;
    if (n2 > 0)
        I(n2 - 1) += i_hist;
}

void series_RL::updateState(const Eigen::VectorXd& V, double dt)
{
    // 获取当前支路两端电压
    double v1 = (n1 > 0) ? V(n1 - 1) : 0.0;
    double v2 = (n2 > 0) ? V(n2 - 1) : 0.0;
    double v_new = v1 - v2;
    // 更新当前电流 i(t) = G_eq * v(t) + I_hist
    i_now = i_hist + G_eq * v_new;
    // 保存当前电压，用于下一步历史项计算
    v_hist = v_new;
}

// --- Inductor ---
Inductor::Inductor(int n1, int n2, double l)
    : n1(n1)
    , n2(n2)
    , L(l)
    , i_hist(0.0)
    , v_hist(0.0)
{
    if (std::abs(l) < 1e-12)
        throw std::runtime_error("Inductor value cannot be zero.");
}

void Inductor::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 梯形法则离散化: i(t) = (dt/2L) * v(t) + I_hist
    // 等效电导 G_eq = dt / (2L)
    G_eq = dt / (2.0 * L);
    if (n1 > 0)
        triplets.emplace_back(n1 - 1, n1 - 1, G_eq);
    if (n2 > 0)
        triplets.emplace_back(n2 - 1, n2 - 1, G_eq);
    if (n1 > 0 && n2 > 0) {
        triplets.emplace_back(n1 - 1, n2 - 1, -G_eq);
        triplets.emplace_back(n2 - 1, n1 - 1, -G_eq);
    }
}

void Inductor::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // 更新历史电流源: I_hist(t) = i(t-dt) + (dt/2L) * v(t-dt)
    // 这里的 v_hist * G_eq = v(t-dt) * (dt/2L)
    i_hist = i_hist + 2 * v_hist * G_eq;
    if (n1 > 0)
        I(n1 - 1) -= i_hist;
    if (n2 > 0)
        I(n2 - 1) += i_hist;
}

void Inductor::updateState(const Eigen::VectorXd& V, double dt)
{
    double v1 = (n1 > 0) ? V(n1 - 1) : 0.0;
    double v2 = (n2 > 0) ? V(n2 - 1) : 0.0;
    double v_new = v1 - v2;
    // 计算当前电流 i(t) = G_eq * v(t) + I_hist
    i_now = i_hist + G_eq * v_new;
    // 保存当前电压，用于下一步历史项计算
    v_hist = v_new;
}

// --- Capacitor ---
Capacitor::Capacitor(int n1, int n2, double c)
    : n1(n1)
    , n2(n2)
    , C(c)
    , i_hist(0.0)
    , v_hist(0.0)
{
    if (std::abs(c) < 1e-15)
        throw std::runtime_error("Capacitor value cannot be zero.");
}

void Capacitor::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 梯形法则离散化: i(t) = (2C/dt) * v(t) + I_hist
    // 等效电导 G_eq = 2C / dt
    G_eq = 2.0 * C / dt;
    if (n1 > 0)
        triplets.emplace_back(n1 - 1, n1 - 1, G_eq);
    if (n2 > 0)
        triplets.emplace_back(n2 - 1, n2 - 1, G_eq);
    if (n1 > 0 && n2 > 0) {
        triplets.emplace_back(n1 - 1, n2 - 1, -G_eq);
        triplets.emplace_back(n2 - 1, n1 - 1, -G_eq);
    }
}

void Capacitor::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // 更新历史电流源: I_hist(t) = -i(t-dt) - (2C/dt) * v(t-dt)
    i_hist = -i_hist - 2 * v_hist * G_eq;
    if (n1 > 0)
        I(n1 - 1) -= i_hist;
    if (n2 > 0)
        I(n2 - 1) += i_hist;
}

void Capacitor::updateState(const Eigen::VectorXd& V, double dt)
{
    double v1 = (n1 > 0) ? V(n1 - 1) : 0.0;
    double v2 = (n2 > 0) ? V(n2 - 1) : 0.0;
    double v_new = v1 - v2;
    // 计算当前电流 i(t) = G_eq * v(t) + I_hist
    i_now = G_eq * v_new + i_hist;
    // 保存当前电压，用于下一步历史项计算
    v_hist = v_new;
}

// --- TwoValueResistor ---
TwoValueResistor::TwoValueResistor(int n1, int n2, double r_closed, double r_open)
    : n1(n1)
    , n2(n2)
    , is_closed(false) // 默认初始为断开状态
{
    G_closed = 1.0 / r_closed;
    G_open = 1.0 / r_open;
}

void TwoValueResistor::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 根据当前开断状态选择对应的电导值进行印章
    double G = is_closed ? G_closed : G_open;
    if (n1 > 0) {
        triplets.emplace_back(n1 - 1, n1 - 1, G);
    }
    if (n2 > 0) {
        triplets.emplace_back(n2 - 1, n2 - 1, G);
    }
    if (n1 > 0 && n2 > 0) {
        triplets.emplace_back(n1 - 1, n2 - 1, -G);
        triplets.emplace_back(n2 - 1, n1 - 1, -G);
    }
}

void TwoValueResistor::updateState(const Eigen::VectorXd& V, double dt)
{
    double v1 = (n1 > 0) ? V(n1 - 1) : 0.0;
    double v2 = (n2 > 0) ? V(n2 - 1) : 0.0;
    double v_new = v1 - v2;
    double G = is_closed ? G_closed : G_open;
    i_now = G * v_new;
}

void TwoValueResistor::setState(bool closed) { is_closed = closed; }
double TwoValueResistor::getCurrentResistance() const { return 1.0 / (is_closed ? G_closed : G_open); }

// --- CircuitBreakerPhase ---
CircuitBreakerPhase::CircuitBreakerPhase(int n1, int n2,
    double r_closed, double r_open,
    bool initiallyClosed, double freq_Hz)
    : n1_(n1)
    , n2_(n2)
    , freq_(freq_Hz)
    , sw_(n1, n2, r_closed, r_open)
    , closed_(initiallyClosed)
{
    sw_.setState(closed_);
}

void CircuitBreakerPhase::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 断路器的印章操作委托给内部的双值电阻
    sw_.stamp(triplets, dt);
}

void CircuitBreakerPhase::updateState(const Eigen::VectorXd& V, double dt)
{
    // 在更新状态前，记录上一步的电流，用于过零检测
    i_prev_ = sw_.get_I();
    sw_.updateState(V, dt);
}

void CircuitBreakerPhase::schedule(double t_action, Change change)
{
    // 将一个操作指令添加到待执行列表中
    ops_.push_back(Op { t_action, change, false });
}

bool CircuitBreakerPhase::applyScheduledAt(double t_now)
{
    bool changed = false;
    for (auto& op : ops_) {
        // 检查是否有未执行且已到时的操作
        if (!op.done && t_now + 1e-12 >= op.t_effective) {
            op.done = true;
            if (op.change == Change::OpenToClose) {
                // 合闸操作：如果当前是断开的，则立即闭合
                if (!closed_) {
                    setClosed(true);
                    changed = true; // 导纳矩阵改变，需要重建
                }
            } else { // CloseToOpen
                // 分闸操作：如果当前是闭合的
                if (closed_) {
                    double i_now = sw_.get_I();
                    if (std::abs(i_now) <= I_EPS_) {
                        // 如果电流已经很小（已过零），立即分闸
                        setClosed(false);
                        changed = true;
                    } else {
                        // 如果电流不为零，则进入“等待过零”状态
                        armed_open_ = true;
                    }
                }
            }
        }
    }
    return changed;
}

bool CircuitBreakerPhase::checkZeroCrossAndOpen()
{
    // 如果不处于“等待过零”状态，直接返回
    if (!armed_open_ || !closed_)
        return false;

    double i_now = sw_.get_I();
    // 判断是否过零：前后两步电流异号，或当前电流足够小
    bool crossed = (i_prev_ * i_now <= 0.0) || (std::abs(i_now) <= I_EPS_);
    if (crossed) {
        setClosed(false); // 执行分闸
        armed_open_ = false; // 解除“等待过零”状态
        return true; // 导纳矩阵改变，需要重建
    }
    return false;
}

// --- VoltageSource ---
VoltageSource::VoltageSource(int n_pos, int n_neu,
    double v_phase_peak, double freq, double phase_deg,
    double rs, double isc_peak)
    : n_pos(n_pos)
    , n_neu(n_neu)
    , V_peak(v_phase_peak)
    , omega(2.0 * PI * freq)
    , phase_rad(phase_deg * PI / 180.0)
    , Rs(rs)
{
    if (Rs <= 0.0)
        throw std::runtime_error("VoltageSource: Rs must be > 0.");
    if (isc_peak <= 0.0)
        throw std::runtime_error("VoltageSource: isc_peak must be > 0.");
    // 根据短路电流计算内电感: Ls = V_peak / (ω * I_sc_peak)
    Ls = V_peak / (omega * isc_peak);
}

void VoltageSource::allocateNodes(Grid& grid)
{
    if (n_src > 0)
        return; // 防止重复分配
    // 申请一个内部节点 n_src，用于连接理想源和内阻抗
    n_src = grid.addNode();

    // 创建内部元件
    // 串联电感：n_pos --- Ls --- n_src
    series_L = std::make_unique<Inductor>(n_pos, n_src, Ls);
    // 并联内阻：n_src --- Rs --- n_neu (接地)
    rs_to_gnd = std::make_unique<Resistor>(n_src, n_neu, Rs);
}

void VoltageSource::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 电压源的印章操作委托给其内部的 R 和 L 元件
    if (rs_to_gnd)
        rs_to_gnd->stamp(triplets, dt);
    if (series_L)
        series_L->stamp(triplets, dt);
}

void VoltageSource::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // 软启动，在1秒内电压从0线性增加到额定值，避免初始冲击
    double k { 1.0 };
    if (t < 1.0) {
        k = t;
    }

    // 计算理想电压源的诺顿等效电流源: I_src = V_src / Rs
    const double v_src = k * V_peak * std::sin(omega * t + phase_rad);
    const double I_src = v_src / Rs;
    // 将该电流注入到内部节点 n_src 和 n_neu
    if (n_src > 0)
        I(n_src - 1) += I_src;
    if (n_neu > 0)
        I(n_neu - 1) -= I_src;

    // 更新内部 R 和 L 元件的历史项
    if (rs_to_gnd)
        rs_to_gnd->updateHistory(I, t, dt); // Resistor无历史项
    if (series_L)
        series_L->updateHistory(I, t, dt);
}

void VoltageSource::updateState(const Eigen::VectorXd& V, double dt)
{
    // 更新内部元件的状态
    if (rs_to_gnd)
        rs_to_gnd->updateState(V, dt);
    if (series_L)
        series_L->updateState(V, dt);
}

// --- Load ---
Load::Load(std::vector<int> nodes_abc,
    double P_total_MW,
    double Q_total_MVAr,
    double V_LL_kV,
    double freq_Hz,
    int n_neu)
    : nodes(std::move(nodes_abc))
    , n_gnd(n_neu)
    , P_total(P_total_MW)
    , Q_total(std::max(0.0, Q_total_MVAr)) // 仅考虑感性无功
    , V_LL(V_LL_kV)
    , freq(freq_Hz)
{
    if (nodes.size() != 3)
        throw std::runtime_error("Load: nodes_abc must contain 3 node numbers");

    // 计算每相等效 R 和 L
    // R_phase = V_ph^2 / P_ph = (V_LL/√3)^2 / (P_total/3) = V_LL^2 / P_total
    if (P_total > 1e-12) {
        R_phase = (V_LL * V_LL) / P_total; // 单位: [kV^2/MW] -> [Ω]
    }

    if (Q_total > 1e-12) {
        // X_L = V_LL^2 / Q_total
        double X_L = (V_LL * V_LL) / Q_total; // 单位: [kV^2/MVAr] -> [Ω]
        double omega = 2.0 * PI * freq;
        L_phase = X_L / omega; // H
    }

    // 根据计算出的 R, L 值创建内部元件
    if (R_phase > 0) {
        rA = std::make_unique<Resistor>(nodes[0], n_gnd, R_phase);
        rB = std::make_unique<Resistor>(nodes[1], n_gnd, R_phase);
        rC = std::make_unique<Resistor>(nodes[2], n_gnd, R_phase);
    }
    if (L_phase > 0) {
        lA = std::make_unique<Inductor>(nodes[0], n_gnd, L_phase);
        lB = std::make_unique<Inductor>(nodes[1], n_gnd, L_phase);
        lC = std::make_unique<Inductor>(nodes[2], n_gnd, L_phase);
    }
}

// 负荷的 stamp, updateHistory, updateState 均委托给其内部的 R, L 元件
void Load::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    if (rA)
        rA->stamp(triplets, dt);
    if (rB)
        rB->stamp(triplets, dt);
    if (rC)
        rC->stamp(triplets, dt);
    if (lA)
        lA->stamp(triplets, dt);
    if (lB)
        lB->stamp(triplets, dt);
    if (lC)
        lC->stamp(triplets, dt);
}

void Load::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    if (rA)
        rA->updateHistory(I, t, dt);
    if (rB)
        rB->updateHistory(I, t, dt);
    if (rC)
        rC->updateHistory(I, t, dt);
    if (lA)
        lA->updateHistory(I, t, dt);
    if (lB)
        lB->updateHistory(I, t, dt);
    if (lC)
        lC->updateHistory(I, t, dt);
}

void Load::updateState(const Eigen::VectorXd& V, double dt)
{
    if (rA)
        rA->updateState(V, dt);
    if (rB)
        rB->updateState(V, dt);
    if (rC)
        rC->updateState(V, dt);
    if (lA)
        lA->updateState(V, dt);
    if (lB)
        lB->updateState(V, dt);
    if (lC)
        lC->updateState(V, dt);
}

// 获取各相总电流（R支路电流 + L支路电流）
double Load::get_Ia() const
{
    double i = 0.0;
    if (rA)
        i += rA->get_I();
    if (lA)
        i += lA->get_I();
    return i;
}

double Load::get_Ib() const
{
    double i = 0.0;
    if (rB)
        i += rB->get_I();
    if (lB)
        i += lB->get_I();
    return i;
}

double Load::get_Ic() const
{
    double i = 0.0;
    if (rC)
        i += rC->get_I();
    if (lC)
        i += lC->get_I();
    return i;
}

// --- SimpleLine ---
SimpleLine::SimpleLine(
    std::vector<int> nodes_i, std::vector<int> nodes_j,
    double R1, double L1, double C1,
    double R0, double L0, double C0,
    double freq, double dt)
    : nodes_i(std::move(nodes_i))
    , nodes_j(std::move(nodes_j))
{
    // 简单模型只使用正序参数，忽略耦合
    R = R1;
    L = L1;
    C = 0.5 * C1; // PI模型，电容分两半
}

void SimpleLine::allocateNodes(Grid& grid)
{
    // 为三相分别创建串联RL支路
    internal_devices.push_back(std::make_unique<series_RL>(nodes_i[0], nodes_j[0], R, L));
    internal_devices.push_back(std::make_unique<series_RL>(nodes_i[1], nodes_j[1], R, L));
    internal_devices.push_back(std::make_unique<series_RL>(nodes_i[2], nodes_j[2], R, L));

    // 在i端创建三相的对地电容
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_i[0], 0, C));
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_i[1], 0, C));
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_i[2], 0, C));

    // 在j端创建三相的对地电容
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_j[0], 0, C));
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_j[1], 0, C));
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_j[2], 0, C));
}

// SimpleLine 的 stamp, updateHistory, updateState 均委托给其内部元件
void SimpleLine::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    for (auto& dev : internal_devices)
        dev->stamp(triplets, dt);
}

void SimpleLine::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    for (auto& dev : internal_devices)
        dev->updateHistory(I, t, dt);
}

void SimpleLine::updateState(const Eigen::VectorXd& V, double dt)
{
    for (auto& dev : internal_devices)
        dev->updateState(V, dt);
}

// --- PI_line ---

// 辅助函数：根据序参数计算相参数矩阵 (适用于完全换位的线路)
inline Eigen::Matrix3d make_seq_to_phase_real(double M1, double M0)
{
    // M_abc = M1 * I + ((M0 - M1)/3) * J
    // I 是单位矩阵, J 是全1矩阵
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d J;
    J.setOnes();
    return M1 * I + ((M0 - M1) / 3.0) * J;
}

PI_line::PI_line(std::vector<int> nodes_i, std::vector<int> nodes_j,
    double R1, double L1, double C1,
    double R0, double L0, double C0,
    double freq, double dt)
    : nodes_i(std::move(nodes_i))
    , nodes_j(std::move(nodes_j))
{
    // 计算相域参数矩阵
    Rabc = make_seq_to_phase_real(R1, R0);
    Labc = make_seq_to_phase_real(L1, L0);
    Cabc = make_seq_to_phase_real(C1, C0);
    // 计算初始的等效矩阵
    computeMatrices(dt);
}

void PI_line::computeMatrices(double dt)
{
    // 串联RL支路的诺顿等效导纳矩阵: G_series = (R_abc + (2/dt) * L_abc)^-1
    Eigen::Matrix3d A = Rabc + (2.0 / dt) * Labc;
    G_series = A.inverse();

    // 并联电容支路的诺顿等效导纳矩阵 (两端各一半): G_sh_half = (C_abc/2) / (dt/2) = C_abc/dt
    G_sh_half = Cabc / dt;

    // 历史电流衰减矩阵 Alpha = (I - dt/2 * L^-1*R) * (I + dt/2 * L^-1*R)^-1
    Eigen::Matrix3d Linv = Labc.inverse();
    Eigen::Matrix3d K = (dt * 0.5) * (Linv * Rabc);
    Alpha = (Eigen::Matrix3d::Identity() - K) * (Eigen::Matrix3d::Identity() + K).inverse();
}

void PI_line::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 如果时间步长改变，需要重新计算矩阵
    computeMatrices(dt);

    auto stamp_block = [&](int r, int c, double val) {
        if (r > 0 && c > 0)
            triplets.emplace_back(r - 1, c - 1, val);
    };

    // 串联支路 3x3 导纳矩阵的印章
    // G_ii, G_jj 块增加 G_series
    // G_ij, G_ji 块减去 G_series
    for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) {
            double Ypq = G_series(p, q);
            stamp_block(nodes_i[p], nodes_i[q], +Ypq);
            stamp_block(nodes_j[p], nodes_j[q], +Ypq);
            stamp_block(nodes_i[p], nodes_j[q], -Ypq);
            stamp_block(nodes_j[p], nodes_i[q], -Ypq);
        }
    }

    // 并联电容 3x3 导纳矩阵的印章 (两端对地)
    // G_ii, G_jj 块分别增加 G_sh_half
    for (int p = 0; p < 3; ++p) {
        for (int q = 0; q < 3; ++q) {
            double Ypq = G_sh_half(p, q);
            stamp_block(nodes_i[p], nodes_i[q], +Ypq);
            stamp_block(nodes_j[p], nodes_j[q], +Ypq);
        }
    }
}

void PI_line::updateHistory(Eigen::VectorXd& I, double /*t*/, double dt)
{
    // 更新串联支路的历史电流注入
    // i_hist(t) = Alpha * i_hist(t-dt) + (I+Alpha) * G_series * v(t-dt)
    i_hist_series = Alpha * i_hist_series + (Eigen::Matrix3d::Identity() + Alpha) * (G_series * v_prev_series);

    for (int p = 0; p < 3; ++p) {
        int ni = nodes_i[p], nj = nodes_j[p];
        if (ni > 0)
            I(ni - 1) -= i_hist_series(p);
        if (nj > 0)
            I(nj - 1) += i_hist_series(p);
    }

    // 更新并联电容的历史电流注入
    // i_hist(t) = -i_hist(t-dt) - 2 * G_sh * v(t-dt)
    i_hist_i = -i_hist_i - 2.0 * (G_sh_half * v_prev_i);
    i_hist_j = -i_hist_j - 2.0 * (G_sh_half * v_prev_j);
    for (int p = 0; p < 3; ++p) {
        int ni = nodes_i[p], nj = nodes_j[p];
        if (ni > 0)
            I(ni - 1) -= i_hist_i(p);
        if (nj > 0)
            I(nj - 1) -= i_hist_j(p);
    }
}

void PI_line::updateState(const Eigen::VectorXd& V, double /*dt*/)
{
    // 从总电压向量 V 中提取本线路两端的节点电压
    Eigen::Vector3d Vi, Vj;
    for (int p = 0; p < 3; ++p) {
        int ni = nodes_i[p], nj = nodes_j[p];
        Vi(p) = (ni > 0) ? V(ni - 1) : 0.0;
        Vj(p) = (nj > 0) ? V(nj - 1) : 0.0;
    }

    // 更新串联支路状态
    Eigen::Vector3d v_now_series = Vi - Vj;
    i_now_series = G_series * v_now_series + i_hist_series;
    v_prev_series = v_now_series; // 保存当前电压差，用于下一步历史项计算

    // 更新并联支路状态
    v_prev_i = Vi;
    v_prev_j = Vj;
}

// --- The Full HB-π-ZIM Line Model ---本混合模型是邹德虎提出并于2025年8月撰写论文投稿---jin
// ---------- HBPiZimLine 构造（只算参数，不建器件/节点） ----------
HBPiZimLine::HBPiZimLine(
    std::vector<int> nodes_i, std::vector<int> nodes_j,
    double R1, double L1, double C1,
    double R0, double L0, double C0,
    double freq, double dt)
    : nodes_i(std::move(nodes_i))
    , nodes_j(std::move(nodes_j))
{
    double omega = 2.0 * PI * freq;
    double k = 0.99; // Bergeron 部分补偿系数

    // 1. Bergeron 段参数
    z_c = k * L1 / dt;

    // 2. Shunt compensation (并联补偿) 段参数
    double Cb_total = dt * dt / (k * L1);
    double delta_Cs = Cb_total - C0;
    if (delta_Cs <= 0)
        throw std::runtime_error("HBPiZimLine: delta_Cs <= 0, model parameters invalid.");
    delta_L = 1.0 / (omega * omega * delta_Cs); // 补偿电感
    C_m = (C1 - C0) / 3.0; // 相间互容

    // 3. Series synthesis (串联综合) 段参数
    double RS1 = R1, XS1 = (1.0 - k) * L1 * omega;
    double RS0 = R0, XS0 = (L0 - k * L1) * omega;

    // 计算综合网络的RL值 (A,B,C,D)
    double denS1 = RS1 * RS1 + XS1 * XS1;
    double denS0 = RS0 * RS0 + XS0 * XS0;
    A = 3.0 / ((2.0 * RS1) / denS1 + RS0 / denS0);
    B = 3.0 / (omega * ((2.0 * XS1) / denS1 + XS0 / denS0));
    C = 3.0 / (RS0 / denS0 - RS1 / denS1);
    D = 3.0 / (omega * (XS0 / denS0 - XS1 / denS1));

    // 初始化历史状态为零
    v_i_hist.setZero();
    i_i_hist.setZero();
    v_k_hist.setZero();
    i_k_hist.setZero();
}

void HBPiZimLine::allocateNodes(Grid& grid)
{
    if (!nodes_k.empty())
        return; // 防止重复分配

    // 申请三个内部节点 kA, kB, kC
    nodes_k = { grid.addNode(), grid.addNode(), grid.addNode() };

    // 创建并联补偿部分的元件 (连接在 k 节点上)
    for (int n : nodes_k) {
        internal_devices.push_back(std::make_unique<Inductor>(n, 0, delta_L)); // delta_L 到地
    }
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_k[0], nodes_k[1], C_m)); // C_m 跨接
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_k[1], nodes_k[2], C_m));
    internal_devices.push_back(std::make_unique<Capacitor>(nodes_k[2], nodes_k[0], C_m));

    // 创建串联综合部分的元件 (连接在 k 和 j 节点之间)
    auto add_rl_parallel = [&](int n1, int n2, double r, double l) {
        internal_devices.push_back(std::make_unique<Resistor>(n1, n2, r));
        internal_devices.push_back(std::make_unique<Inductor>(n1, n2, l));
    };
    // kX - jX (自支路)
    add_rl_parallel(nodes_k[0], nodes_j[0], A, B);
    add_rl_parallel(nodes_k[1], nodes_j[1], A, B);
    add_rl_parallel(nodes_k[2], nodes_j[2], A, B);
    // kX - jY (互支路)
    add_rl_parallel(nodes_k[0], nodes_j[1], C, D);
    add_rl_parallel(nodes_k[0], nodes_j[2], C, D);
    add_rl_parallel(nodes_k[1], nodes_j[0], C, D);
    add_rl_parallel(nodes_k[1], nodes_j[2], C, D);
    add_rl_parallel(nodes_k[2], nodes_j[0], C, D);
    add_rl_parallel(nodes_k[2], nodes_j[1], C, D);
    // 端内耦合支路
    add_rl_parallel(nodes_k[0], nodes_k[1], -C, -D);
    add_rl_parallel(nodes_k[0], nodes_k[2], -C, -D);
    add_rl_parallel(nodes_k[1], nodes_k[2], -C, -D);
    add_rl_parallel(nodes_j[0], nodes_j[1], -C, -D);
    add_rl_parallel(nodes_j[0], nodes_j[2], -C, -D);
    add_rl_parallel(nodes_j[1], nodes_j[2], -C, -D);
}

void HBPiZimLine::stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt)
{
    // 1. Bergeron 段的印章：在 i 端和 k 端对角线上增加等效电导 1/zc
    double G_b = 1.0 / z_c;
    for (int idx = 0; idx < 3; ++idx) {
        if (nodes_i[idx] > 0)
            triplets.emplace_back(nodes_i[idx] - 1, nodes_i[idx] - 1, G_b);
        if (!nodes_k.empty() && nodes_k[idx] > 0)
            triplets.emplace_back(nodes_k[idx] - 1, nodes_k[idx] - 1, G_b);
    }
    // 2. PI-ZIM 网络部分的印章委托给其内部元件
    for (auto& dev : internal_devices)
        dev->stamp(triplets, dt);
}

void HBPiZimLine::updateHistory(Eigen::VectorXd& I, double t, double dt)
{
    // 1. Bergeron 段的历史电流注入
    // I_hist_i(t) = (1/zc)*V_k(t-Δt) - I_k(t-Δt)
    // I_hist_k(t) = (1/zc)*V_i(t-Δt) - I_i(t-Δt)
    // 注意：这里的电流方向定义可能与传统 Bergeron 模型有符号差异，需与 updateState 保持一致
    Eigen::Vector3d I_hist_i = (1.0 / z_c) * v_k_hist + i_k_hist;
    Eigen::Vector3d I_hist_k = (1.0 / z_c) * v_i_hist + i_i_hist;

    for (int m = 0; m < 3; ++m) {
        if (nodes_i[m] > 0)
            I(nodes_i[m] - 1) += I_hist_i(m);
        if (!nodes_k.empty() && nodes_k[m] > 0)
            I(nodes_k[m] - 1) += I_hist_k(m);
    }

    // 2. PI-ZIM 网络部分的历史注入委托给其内部元件
    for (auto& dev : internal_devices)
        dev->updateHistory(I, t, dt);
}

void HBPiZimLine::updateState(const Eigen::VectorXd& V, double dt)
{
    // 提取 i 端和 k 端的当前电压
    Eigen::Vector3d v_i_new, v_k_new;
    for (int m = 0; m < 3; ++m) {
        v_i_new(m) = (nodes_i[m] > 0) ? V(nodes_i[m] - 1) : 0.0;
        v_k_new(m) = (!nodes_k.empty() && nodes_k[m] > 0) ? V(nodes_k[m] - 1) : 0.0;
    }

    // 计算 i 端和 k 端的当前电流
    // i(t) = (1/zc)*v(t) - I_hist(t)
    Eigen::Vector3d I_hist_i_t = (1.0 / z_c) * v_k_hist + i_k_hist;
    Eigen::Vector3d I_hist_k_t = (1.0 / z_c) * v_i_hist + i_i_hist;

    Eigen::Vector3d i_i_new = (1.0 / z_c) * v_i_new - I_hist_i_t;
    Eigen::Vector3d i_k_new = (1.0 / z_c) * v_k_new - I_hist_k_t;

    // 更新历史状态，用于下一个时间步
    v_i_hist = v_i_new;
    i_i_hist = i_i_new;
    v_k_hist = v_k_new;
    i_k_hist = i_k_new;

    // 更新 PI-ZIM 网络内部元件的状态
    for (auto& dev : internal_devices)
        dev->updateState(V, dt);
}