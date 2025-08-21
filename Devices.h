#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <string>
#include <vector>

// 前向声明 Grid 类，避免循环包含
class Grid;

/**
 * @class Device
 * @brief 所有电路元件的抽象基类。
 *
 * 定义了所有元件必须实现的通用接口，包括：
 * - stamp: 将元件对系统导纳矩阵的贡献（“印章”）添加到三元组列表中。
 * - updateHistory: 更新与元件相关的历史等效电流源注入。
 * - updateState: 在求解节点电压后，更新元件的内部状态（如电流、电压）。
 * - allocateNodes: 为复合元件（如电压源、线路模型）分配内部节点。
 */
class Device {
public:
    /// @brief 虚析构函数，确保派生类对象能被正确销毁。
    virtual ~Device() = default;

    /**
     * @brief 将元件的导纳矩阵贡献添加到三元组列表。
     * @param triplets 用于构建稀疏矩阵的三元组向量。
     * @param dt 仿真时间步长 (s)。
     */
    virtual void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) = 0;

    /**
     * @brief 更新历史等效电流源注入向量 I。
     * @param I 节点电流注入向量，将被此方法修改。
     * @param t 当前仿真时间 (s)。
     * @param dt 仿真时间步长 (s)。
     */
    virtual void updateHistory(Eigen::VectorXd& I, double t, double dt) = 0;

    /**
     * @brief 根据新计算出的节点电压 V 更新元件内部状态。
     * @param V 节点电压向量。
     * @param dt 仿真时间步长 (s)。
     */
    virtual void updateState(const Eigen::VectorXd& V, double dt) = 0;

    /**
     * @brief 为需要内部节点的复合元件分配节点。
     * @param grid Grid 对象的引用，用于添加新节点。
     */
    virtual void allocateNodes(Grid& grid) { }

    /**
     * @brief 获取流过元件的当前电流值。
     * @return 当前电流 (A)。
     */
    virtual double get_I() { return 0.0; }
};

// --- 基本元件 ---

/**
 * @class Resistor
 * @brief 线性电阻元件。
 */
class Resistor : public Device {
public:
    /**
     * @brief 构造一个电阻。
     * @param n1 节点1的编号。
     * @param n2 节点2的编号。
     * @param r 电阻值 (Ω)。
     */
    Resistor(int n1, int n2, double r);
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override { } // 纯电阻无历史项
    void updateState(const Eigen::VectorXd& V, double dt);
    double get_I() { return i_now; }

private:
    int n1, n2; ///< 连接节点
    double G; ///< 电导 (S)
    double i_now { 0.0 }; ///< 当前流过电阻的电流 (A)
};

/**
 * @class series_RL
 * @brief 串联 RL 支路。
 * 采用梯形积分法进行离散化。
 */
class series_RL : public Device {
public:
    /**
     * @brief 构造串联RL支路。
     * @param n1 节点1。
     * @param n2 节点2。
     * @param r 电阻 (Ω)。
     * @param l 电感 (H)。
     */
    series_RL(int n1, int n2, double r, double l);
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
    double get_I() { return i_now; }

private:
    int n1, n2; ///< 连接节点
    double R; ///< 电阻值 (Ω)
    double L; ///< 电感值 (H)
    double G_eq; ///< 诺顿等效电导 (S)
    double value1_RL; ///< 历史电流计算系数
    double i_hist { 0.0 }; ///< 历史等效电流源 (A)
    double v_hist { 0.0 }; ///< 上一步的支路电压 (V)
    double i_now { 0.0 }; ///< 当前支路电流 (A)
};

/**
 * @class Inductor
 * @brief 理想电感元件。
 * 采用梯形积分法离散化，等效为诺顿电路（等效电导 + 历史电流源）。
 */
class Inductor : public Device {
public:
    /**
     * @brief 构造一个电感。
     * @param n1 节点1。
     * @param n2 节点2。
     * @param l 电感值 (H)。
     */
    Inductor(int n1, int n2, double l);
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
    double get_I() { return i_now; }

private:
    int n1, n2; ///< 连接节点
    double L; ///< 电感值 (H)
    double G_eq; ///< 等效电导 (S)
    double i_hist { 0.0 }; ///< 历史等效电流源 (A)
    double v_hist { 0.0 }; ///< 上一步的元件电压 (V)
    double i_now { 0.0 }; ///< 当前流过电感的电流 (A)
};

/**
 * @class Capacitor
 * @brief 理想电容元件。
 * 采用梯形积分法离散化，等效为诺顿电路（等效电导 + 历史电流源）。
 */
class Capacitor : public Device {
public:
    /**
     * @brief 构造一个电容。
     * @param n1 节点1。
     * @param n2 节点2。
     * @param c 电容值 (F)。
     */
    Capacitor(int n1, int n2, double c);
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
    double get_I() { return i_now; }

private:
    int n1, n2; ///< 连接节点
    double C; ///< 电容值 (F)
    double G_eq; ///< 等效电导 (S)
    double i_hist { 0.0 }; ///< 历史等效电流源 (A)
    double v_hist { 0.0 }; ///< 上一步的元件电压 (V)
    double i_now { 0.0 }; ///< 当前流过电容的电流 (A)
};

// --- 复合与开关元件 ---

/**
 * @class TwoValueResistor
 * @brief 具有两种阻值（例如开/关）的电阻。
 * 用于模拟开关、故障等。
 */
class TwoValueResistor : public Device {
public:
    /**
     * @brief 构造一个双值电阻。
     * @param n1 节点1。
     * @param n2 节点2。
     * @param r_closed 闭合状态的电阻 (Ω)。
     * @param r_open 打开状态的电阻 (Ω)。
     */
    TwoValueResistor(int n1, int n2, double r_closed, double r_open);
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override { }
    void updateState(const Eigen::VectorXd& V, double dt);
    void setState(bool closed);

    /// @brief 获取当前等效电阻值。
    double getCurrentResistance() const;
    double get_I() { return i_now; }

private:
    double i_now { 0.0 }; ///< 当前电流 (A)
    int n1, n2; ///< 连接节点
    double G_closed, G_open; ///< 闭合和打开状态的电导 (S)
    bool is_closed; ///< 当前是否处于闭合状态
};

/**
 * @class CircuitBreakerPhase
 * @brief 单相断路器模型。
 * 内部使用 TwoValueResistor 实现，并增加了时序控制和电流过零开断逻辑。
 */
class CircuitBreakerPhase : public Device {
public:
    /// @brief 定义断路器动作类型
    enum class Change { CloseToOpen,
        OpenToClose };

    /**
     * @brief 构造一个单相断路器。
     * @param n1 节点1。
     * @param n2 节点2。
     * @param r_closed 合闸等效电阻 (Ω)，应非常小。
     * @param r_open 分闸等效电阻 (Ω)，应非常大。
     * @param initiallyClosed 初始状态是否为闭合。
     * @param freq_Hz 系统频率 (Hz)，用于内部逻辑（当前未使用，但为未来扩展保留）。
     */
    CircuitBreakerPhase(int n1, int n2,
        double r_closed, double r_open,
        bool initiallyClosed, double freq_Hz);

    // --- Device 接口实现 ---
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override
    { /* 无历史项 */
    }
    void updateState(const Eigen::VectorXd& V, double dt) override;

    /**
     * @brief 安排一次未来的动作。
     * @param t_action 动作发生的目标时间 (s)。
     * @param change 动作类型 (合闸或分闸)。
     */
    void schedule(double t_action, Change change);

    /**
     * @brief 在当前时间步检查并应用已到期的预定动作。
     * @param t_now 当前仿真时间 (s)。
     * @return 如果状态改变导致导纳矩阵需要重建，则返回 true。
     */
    bool applyScheduledAt(double t_now);

    /**
     * @brief 在电压求解后调用，检查是否满足过零开断条件。
     * 如果断路器处于“等待过零”状态且电流已过零，则执行开断。
     * @return 如果状态改变导致导纳矩阵需要重建，则返回 true。
     */
    bool checkZeroCrossAndOpen();

    // --- 状态查询 ---
    bool isClosed() const { return closed_; }
    double get_I() { return sw_.get_I(); }

private:
    /// @brief 预定操作的结构体
    struct Op {
        double t_effective; ///< 操作生效时间
        Change change; ///< 操作类型
        bool done { false }; ///< 是否已完成
    };

    int n1_, n2_; ///< 连接节点
    double freq_; ///< 系统频率
    TwoValueResistor sw_; ///< 内部开关元件
    bool closed_; ///< 当前是否闭合
    double i_prev_ { 0.0 }; ///< 上一步的电流，用于过零检测
    bool armed_open_ { false }; ///< 是否已收到分闸指令并等待电流过零
    const double I_EPS_ = 1e-6; ///< 电流过零判断阈值 (A)
    std::vector<Op> ops_; ///< 预定操作列表

    /// @brief 内部辅助函数，用于改变开关状态
    void setClosed(bool c)
    {
        closed_ = c;
        sw_.setState(c);
    }
};

/**
 * @class VoltageSource
 * @brief 理想正弦电压源，带戴维南等效内阻抗 (Rs + Ls)。
 * 内部通过诺顿等效电路实现：一个理想电流源并联一个电阻 Rs，再串联一个电感 Ls。
 */
class VoltageSource : public Device {
public:
    /**
     * @brief 构造一个电压源。
     * @param n_pos 正极节点。
     * @param n_neu 负极/中性点节点。
     * @param v_phase_peak 相电压峰值 (V)。
     * @param freq 频率 (Hz)。
     * @param phase_deg 相位角 (度)。
     * @param rs 串联内阻 (Ω)。
     * @param isc_peak 短路电流峰值 (A)，用于计算内电感 Ls。
     */
    VoltageSource(int n_pos, int n_neu,
        double v_phase_peak, double freq, double phase_deg,
        double rs, double isc_peak);

    void allocateNodes(Grid& grid) override;
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

    double get_I()
    {
        return series_L->get_I();
    }

private:
    int n_pos, n_neu; ///< 外部连接节点
    int n_src = -1; ///< 内部理想源节点

    double V_peak, omega, phase_rad; ///< 电压源参数
    double Rs, Ls; ///< 内阻抗参数

    // 内部实现所用的基本元件
    std::unique_ptr<class Resistor> rs_to_gnd; ///< 并联内阻
    std::unique_ptr<class Inductor> series_L; ///< 串联内感
};

/**
 * @class Load
 * @brief 三相等平衡RL负荷模型。
 * 根据总有功、无功和额定线电压计算每相的 R 和 L，并创建相应的并联支路。
 */
class Load : public Device {
public:
    /**
     * @brief 构造一个三相负荷。
     * @param nodes_abc 三相连接节点 [A, B, C]。
     * @param P_total_MW 三相总有功功率 (MW)。
     * @param Q_total_MVAr 三相总无功功率 (MVAr)。
     * @param V_LL_kV 额定线电压 (kV)。
     * @param freq_Hz 额定频率 (Hz)。
     * @param n_neu 中性点/接地节点 (默认为0)。
     */
    Load(std::vector<int> nodes_abc,
        double P_total_MW,
        double Q_total_MVAr,
        double V_LL_kV,
        double freq_Hz,
        int n_neu = 0);

    void allocateNodes(Grid& grid) override { } // 无内部新增节点
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

    // --- 电流查询 ---
    double get_Ia() const;
    double get_Ib() const;
    double get_Ic() const;

private:
    std::vector<int> nodes; ///< [na, nb, nc]
    int n_gnd { 0 }; ///< 接地节点
    double P_total, Q_total, V_LL, freq;
    double R_phase { 0.0 }; ///< 每相等效电阻 (Ω)
    double L_phase { 0.0 }; ///< 每相等效电感 (H)

    // 每相的内部元件（如果 P 或 Q 为0，则对应指针为空）
    std::unique_ptr<class Resistor> rA, rB, rC;
    std::unique_ptr<class Inductor> lA, lB, lC;
};

/**
 * @class SimpleLine
 * @brief 简化的三相PI型线路模型（忽略相间耦合）。
 * 每相都等效为一个独立的PI型电路（串联RL，两端各接一半对地电容）。
 * @deprecated 该模型已被耦合的 PI_line 替代，仅用于教学或对比。
 */
class SimpleLine : public Device {
public:
    SimpleLine(
        std::vector<int> nodes_i, std::vector<int> nodes_j,
        double R1, double L1, double C1,
        double R0, double L0, double C0,
        double freq, double dt);

    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
    void allocateNodes(Grid& grid) override;

    // --- 电流查询 ---
    double get_Ia() { return internal_devices[0]->get_I(); }
    double get_Ib() { return internal_devices[1]->get_I(); }
    double get_Ic() { return internal_devices[2]->get_I(); }

private:
    std::vector<int> nodes_i, nodes_j;
    std::vector<std::unique_ptr<Device>> internal_devices;
    double R, L, C; // 正序参数
};

/**
 * @class PI_line
 * @brief 考虑相间耦合的集中参数三相PI型线路模型。
 * 使用3x3的相域参数矩阵 (Rabc, Labc, Cabc) 来表示串联阻抗和并联导纳。
 */
class PI_line : public Device {
public:
    /**
     * @brief 构造一个耦合PI线路模型。
     * @param nodes_i i端三相节点 [iA, iB, iC]。
     * @param nodes_j j端三相节点 [jA, jB, jC]。
     * @param R1, L1, C1 正序参数。
     * @param R0, L0, C0 零序参数。
     * @param freq 频率 (Hz)。
     * @param dt 时间步长 (s)。
     */
    PI_line(std::vector<int> nodes_i, std::vector<int> nodes_j,
        double R1, double L1, double C1,
        double R0, double L0, double C0,
        double freq, double dt);

    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;
    void allocateNodes(Grid& grid) override { } // 无内部新节点

    // --- 电流查询 (串联支路电流) ---
    double get_Ia() { return i_now_series(0); }
    double get_Ib() { return i_now_series(1); }
    double get_Ic() { return i_now_series(2); }

private:
    std::vector<int> nodes_i, nodes_j; // i端(A,B,C), j端(A,B,C)

    // 相域参数矩阵
    Eigen::Matrix3d Rabc, Labc, Cabc;
    // 离散化后的等效导纳和系数矩阵
    Eigen::Matrix3d G_series; // 串联支路诺顿等效电导 (R + 2/dt L)^-1
    Eigen::Matrix3d G_sh_half; // 半个并联支路诺顿等效电导 C/dt
    Eigen::Matrix3d Alpha; // 串联支路历史电流衰减矩阵

    // 状态变量
    Eigen::Vector3d v_prev_series = Eigen::Vector3d::Zero();
    Eigen::Vector3d i_hist_series = Eigen::Vector3d::Zero();
    Eigen::Vector3d i_now_series = Eigen::Vector3d::Zero();

    Eigen::Vector3d v_prev_i = Eigen::Vector3d::Zero(), i_hist_i = Eigen::Vector3d::Zero();
    Eigen::Vector3d v_prev_j = Eigen::Vector3d::Zero(), i_hist_j = Eigen::Vector3d::Zero();

    /// @brief 根据时间步长 dt 计算内部等效矩阵。
    void computeMatrices(double dt);
};

/**
 * @class HBPiZimLine
 * @brief 混合 Bergeron-PI-ZIM 线路模型。
 * 该模型由邹德虎提出，旨在结合分布参数模型和集总参数模型的优点，已经于2025年8月撰写论文并投稿，不得侵犯知识产权!
 * 结构：i端 <--- Bergeron ---> k端(内部) <--- PI-ZIM 网络 ---> j端
 */
class HBPiZimLine : public Device {
public:
    HBPiZimLine(
        std::vector<int> nodes_i, std::vector<int> nodes_j,
        double R1, double L1, double C1,
        double R0, double L0, double C0,
        double freq, double dt);

    void allocateNodes(Grid& grid) override;
    void stamp(std::vector<Eigen::Triplet<double>>& triplets, double dt) override;
    void updateHistory(Eigen::VectorXd& I, double t, double dt) override;
    void updateState(const Eigen::VectorXd& V, double dt) override;

private:
    std::vector<int> nodes_i, nodes_j;
    std::vector<int> nodes_k; // 内部三相节点 kA, kB, kC

    // Bergeron 段 (i <-> k) 的参数和历史状态
    double z_c; // 特性阻抗
    Eigen::Vector3d v_i_hist, i_i_hist; // i 端的历史电压和电流
    Eigen::Vector3d v_k_hist, i_k_hist; // k 端的历史电压和电流

    // PI-ZIM 网络的参数 (在构造时计算)
    double A, B, C, D;
    double delta_L, C_m;

    // 存储 PI-ZIM 网络所有元件的容器
    std::vector<std::unique_ptr<Device>> internal_devices;
};