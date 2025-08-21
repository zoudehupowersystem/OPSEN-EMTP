#pragma once
#include "Control.h"
#include "Curve.h"
#include "Devices.h"
#include "Grid.h"
#include <Eigen/Sparse>
#include <vector>

/**
 * @class Simulation
 * @brief 仿真器主类，负责驱动整个电磁暂态仿真过程。
 *
 * 管理主时间步循环，处理事件，调用 Grid 和求解器，并记录数据。
 */
class Simulation {
public:
    /**
     * @brief 构造仿真器。
     * @param ctrl 仿真控制参数对象。
     * @param grid 电路网络对象。
     * @param curve 曲线记录对象。
     * @param fault_switch 指向用于施加故障的开关元件的指针。
     * @param fault_node 故障施加的节点编号。
     */
    Simulation(Control& ctrl, Grid& grid, Curve& curve, TwoValueResistor* fault_switch, int fault_node);

    /// @brief 启动并运行整个仿真。
    void run();

    /**
     * @brief 注册一个断路器到仿真器中，以便在循环中处理其时序逻辑。
     * @param brk 指向 CircuitBreakerPhase 对象的指针。
     */
    void addBreaker(CircuitBreakerPhase* brk) { breakers.push_back(brk); }

private:
    Control& ctrl;
    Grid& grid;
    Curve& curve;
    TwoValueResistor* fault_switch; ///< 指向故障开关
    int fault_node; ///< 故障节点

    /// @brief 稀疏矩阵线性方程组求解器 (LU分解)
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    /// @brief 存储所有需要时序控制的断路器
    std::vector<CircuitBreakerPhase*> breakers;
};
