#pragma once
#include "Control.h"
#include "Devices.h"
#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

/**
 * @class Curve
 * @brief 负责在仿真过程中采样数据并将其写入文件。
 */
class Curve {
public:
    /**
     * @brief 构造曲线记录器。
     * @param ctrl 仿真控制参数，用于获取绘图时间范围和电压曲线设置。
     * @param fault_switch (未使用，可移除) 指向故障开关的指针。
     */
    Curve(const Control& ctrl, const TwoValueResistor* fault_switch);

    /// @brief 析构函数，确保文件流被正确关闭。
    ~Curve();

    /**
     * @brief 在每个时间步调用，对需要的数据进行采样。
     * @param t 当前仿真时间。
     * @param V 当前的节点电压解向量。
     */
    void sample(double t, const Eigen::VectorXd& V);

    /**
     * @brief 注册一个需要记录的电流波形。
     * @param name 电流曲线的名称（将作为文件头）。
     * @param reader 一个函数对象（如此处的lambda），调用它时返回当前的电流值。
     */
    void addCurrentTrace(const std::string& name, std::function<double()> reader)
    {
        curr_traces.push_back({ name, std::move(reader) });
    }

private:
    /// @brief 用于存储电流曲线信息的结构体
    struct CurrTrace {
        std::string name;
        std::function<double()> reader; // 用于获取实时电流值的回调函数
    };

    const Control& ctrl;
    const TwoValueResistor* fault_switch;

    std::ofstream fileV, fileI; // 电压和电流数据文件的输出流
    bool header_written_V { false }; // 标记电压文件头是否已写入
    bool header_written_I { false }; // 标记电流文件头是否已写入

    std::vector<CurrTrace> curr_traces; // 存储所有已注册的电流曲线
};