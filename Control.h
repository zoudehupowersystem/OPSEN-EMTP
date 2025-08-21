#pragma once
#include <string>
#include <variant>
#include <vector>

/**
 * @struct FaultEvent
 * @brief 定义一个故障事件。
 */
struct FaultEvent {
    double time; ///< 事件发生的时间 (s)
    int node; ///< 事件发生的节点编号
    bool apply; ///< true 表示施加故障，false 表示移除故障
};

// 使用 variant 可以方便地在未来扩展更多类型的事件
using Event = std::variant<FaultEvent>;

/**
 * @class Control
 * @brief 存储和管理仿真的控制参数和事件序列。
 *
 * 相当于一个“剧本”，定义了仿真的时长、步长、绘图范围以及期间发生的各种事件。
 */
class Control {
public:
    // --- 仿真基本参数 ---
    double dt = 50e-6; ///< 仿真时间步长 (s)
    double t_end = 3.3; ///< 仿真结束时间 (s)
    double plot_t_start = 1.9; ///< 曲线记录开始时间 (s)
    double plot_t_end = 3.3; ///< 曲线记录结束时间 (s)

    // --- 事件和绘图设置 ---
    /// @brief 存储所有预定事件的列表
    std::vector<Event> events;

    /// @brief 存储需要记录电压波形的节点信息 (曲线名称, 节点编号)
    std::vector<std::pair<std::string, int>> traces;

    /**
     * @brief 设置一个简单的测试用例。
     * 在此函数中定义具体的故障事件和要观察的波形。
     */
    void setupSimpleCase();
};