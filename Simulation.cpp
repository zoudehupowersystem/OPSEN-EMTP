#include "Simulation.h"
#include <chrono>
#include <iostream>
#include <sstream> 
Simulation::Simulation(Control& ctrl, Grid& grid, Curve& curve, TwoValueResistor* fault_switch, int fault_node)
    : ctrl(ctrl)
    , grid(grid)
    , curve(curve)
    , fault_switch(fault_switch)
    , fault_node(fault_node)
{
}

void Simulation::run()
{
    std::cout << "Starting simulation..." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now();

    // 仿真初始化
    grid.allocateInternalNodes(); // 确定最终节点数
    grid.buildMatrix(ctrl.dt); // 构建初始导纳矩阵
    // grid.dumpMatrix("G at t=0 (pre-fault)"); // 调试用，打印初始矩阵

    // 对矩阵进行一次符号分解 (analyzePattern) 和数值分解 (factorize)
    // 如果矩阵结构不变，后续只需调用 factorize
    solver.analyzePattern(grid.getG());
    solver.factorize(grid.getG());

    int num_steps = static_cast<int>(ctrl.t_end / ctrl.dt);
    Eigen::VectorXd V(grid.getNumNodes()); // 节点电压向量
    V.setZero();
    Eigen::VectorXd I(grid.getNumNodes()); // 节点注入电流向量

    // --- 主时间步循环 ---
    for (int i = 0; i <= num_steps; ++i) {
        I.setZero();
        V.setZero();
        double t = i * ctrl.dt;
        bool matrix_changed = false; // 标记本步导纳矩阵是否发生变化

        // 1. 处理预定事件（例如故障）
        for (const auto& event_variant : ctrl.events) {
            const auto& event = std::get<FaultEvent>(event_variant);
            // 检查是否到达故障发生时间
            if (event.node == this->fault_node && std::abs(t - event.time) < ctrl.dt / 2.0) {
                if (fault_switch) {
                    std::cout << "Time: " << t << "s - Fault event on node " << event.node
                              << ": " << (event.apply ? "ON" : "OFF") << std::endl;
                    fault_switch->setState(event.apply); // 改变故障开关状态
                    matrix_changed = true; // 矩阵需要重建
                }
            }
        }

        // 2. 处理断路器的预定动作
        for (auto* brk : breakers) {
            if (brk && brk->applyScheduledAt(t)) {
                matrix_changed = true; // 如果断路器状态改变，矩阵需要重建
            }
        }

        // 3. 如果导纳矩阵改变，则重建并重新分解
        if (matrix_changed) {
            grid.buildMatrix(ctrl.dt);
            // 矩阵结构可能不变，但数值已变，需重新数值分解
            solver.factorize(grid.getG());
        }

        // 4. 更新历史电流注入向量 I
        grid.updateHistoryVector(I, t, ctrl.dt);

        // 5. 求解线性方程组 G * V = I
        V = solver.solve(I);

        // 6. 更新所有设备的状态
        grid.updateDeviceStates(V, ctrl.dt);

        // 7. 处理断路器过零开断逻辑
        bool changed_post = false;
        for (auto* brk : breakers) {
            if (brk && brk->checkZeroCrossAndOpen()) {
                changed_post = true; // 如果发生过零开断，矩阵需要重建
            }
        }
        if (changed_post) {
            grid.buildMatrix(ctrl.dt);
            solver.factorize(grid.getG());
        }

        // 8. 采样并记录数据
        curve.sample(t, V);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    std::cout << "Simulation finished in " << elapsed.count() << " seconds." << std::endl;
}