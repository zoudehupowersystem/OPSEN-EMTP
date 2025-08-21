#include "Grid.h"
#include <iomanip>
#include <iostream>

Grid::Grid(int num_nodes)
    : num_nodes(num_nodes)
{
    // 初始化导纳矩阵的大小
    G.resize(num_nodes, num_nodes);
}

void Grid::allocateInternalNodes()
{
    // 遍历所有设备，让需要内部节点的复合设备（如电压源、复杂线路模型）
    // 向Grid申请并注册它们需要的额外节点。
    for (auto& dev : devices)
        dev->allocateNodes(*this);
}

int Grid::addNode()
{
    // 增加节点计数，并返回新的节点编号（1-based）
    ++num_nodes;
    return num_nodes;
}

void Grid::addDevice(std::unique_ptr<Device> device)
{
    // 将设备的所有权转移到Grid的设备列表中
    devices.push_back(std::move(device));
}

void Grid::buildMatrix(double dt)
{
    // 使用三元组（row, col, value）列表来高效构建稀疏矩阵
    std::vector<Eigen::Triplet<double>> triplets;
    // 遍历所有设备，调用它们的stamp方法，将各自的导纳贡献填充到三元组列表中
    for (const auto& device : devices)
        device->stamp(triplets, dt);

    // 根据最新的节点总数（可能在allocateInternalNodes后改变）重置矩阵大小
    G.resize(num_nodes, num_nodes);
    // 从三元组列表构建稀疏导纳矩阵G。相同位置的元素会自动累加。
    G.setFromTriplets(triplets.begin(), triplets.end());
}

void Grid::updateHistoryVector(Eigen::VectorXd& I, double t, double dt)
{
    // 初始化节点电流注入向量为零
    I.setZero(num_nodes);
    // 遍历所有设备，调用它们的updateHistory方法，将历史电流源注入到向量I中
    for (auto& device : devices) {
        device->updateHistory(I, t, dt);
    }
}

void Grid::updateDeviceStates(const Eigen::VectorXd& V, double dt)
{
    // 在节点电压V求解完毕后，遍历所有设备，更新它们各自的内部状态
    // （例如，计算当前流过元件的电流）
    for (auto& device : devices) {
        device->updateState(V, dt);
    }
}

void Grid::dumpMatrix(const std::string& title) const
{
    // 打印矩阵信息，用于调试
    std::cout << "\n==== " << title << " | size = "
              << G.rows() << "x" << G.cols() << " ====\n";

    // 打印稠密矩阵形式（仅适用于小系统）
    Eigen::MatrixXd D = Eigen::MatrixXd(G);
    std::cout << std::fixed << std::setprecision(6);
    for (int i = 0; i < D.rows(); ++i) {
        for (int j = 0; j < D.cols(); ++j) {
            std::cout << std::setw(14) << D(i, j);
        }
        std::cout << '\n';
    }

    // 打印稀疏矩阵的非零元素列表（更通用）
    std::cout << "\n-- nonzeros (row col value) --\n";
    std::cout << std::setprecision(12);
    for (int k = 0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G, k); it; ++it) {
            std::cout << std::setw(3) << (it.row() + 1) << " "
                      << std::setw(3) << (it.col() + 1) << "  "
                      << it.value() << "\n";
        }
    }
    std::cout.flush();
}