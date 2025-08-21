#pragma once
#include "Devices.h"
#include <Eigen/Sparse>
#include <memory>
#include <vector>

/**
 * @class Grid
 * @brief 管理整个电路网络的核心类。
 *
 * 负责存储所有设备，构建和维护系统导纳矩阵 G，
 * 并提供接口来更新历史电流向量和设备状态。
 */
class Grid {
public:
    /**
     * @brief 构造一个电网。
     * @param num_nodes 初始的外部节点数量（不包括元件内部可能增加的节点）。
     */
    Grid(int num_nodes);

    /**
     * @brief 添加一个设备到电网中。
     * @param device 指向设备对象的 unique_ptr，所有权将转移给 Grid。
     */
    void addDevice(std::unique_ptr<Device> device);

    /**
     * @brief 遍历所有设备，执行它们的 allocateNodes 方法。
     * 这一步必须在所有设备添加后、构建矩阵前调用，以确定最终的节点总数。
     */
    void allocateInternalNodes();

    /**
     * @brief 增加一个新节点到系统中。
     * @return 新节点的编号（1-based）。
     */
    int addNode();

    /**
     * @brief 根据所有设备的状态构建（或重建）系统导纳矩阵 G。
     * @param dt 当前仿真时间步长。
     */
    void buildMatrix(double dt);

    /**
     * @brief 计算并更新整个系统的历史电流注入向量 I。
     * @param I 将被填充的电流向量。
     * @param t 当前仿真时间。
     * @param dt 仿真时间步长。
     */
    void updateHistoryVector(Eigen::VectorXd& I, double t, double dt);

    /**
     * @brief 在求解出节点电压后，更新所有设备各自的内部状态。
     * @param V 节点电压解向量。
     * @param dt 仿真时间步长。
     */
    void updateDeviceStates(const Eigen::VectorXd& V, double dt);

    /// @brief 获取对系统导纳矩阵 G 的 const 引用。
    const Eigen::SparseMatrix<double>& getG() const { return G; }

    /// @brief 获取当前系统中的总节点数。
    int getNumNodes() const { return num_nodes; }

    /**
     * @brief 打印导纳矩阵的内容到控制台，用于调试。
     * @param title 打印时显示的标题。
     */
    void dumpMatrix(const std::string& title) const;

private:
    int num_nodes; ///< 系统中的总节点数
    Eigen::SparseMatrix<double> G; ///< 系统导纳矩阵
    std::vector<std::unique_ptr<Device>> devices; ///< 存储所有设备的容器
};