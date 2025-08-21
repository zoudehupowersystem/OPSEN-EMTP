#include "Curve.h"
#include <iomanip>
#include <iostream>

Curve::Curve(const Control& ctrl, const TwoValueResistor* fault_switch)
    : ctrl(ctrl)
    , fault_switch(fault_switch)
{
    // 打开用于写入电压和电流数据的文件
    fileV.open("curve_V.dat");
    if (!fileV.is_open()) {
        std::cerr << "Error: Could not open curve_V.dat for writing." << std::endl;
    }
    fileI.open("curve_I.dat");
    if (!fileI.is_open()) {
        std::cerr << "Error: Could not open curve_I.dat for writing." << std::endl;
    }
}

Curve::~Curve()
{
    // 在对象销毁时关闭文件
    if (fileV.is_open())
        fileV.close();
    if (fileI.is_open())
        fileI.close();
}

void Curve::sample(double t, const Eigen::VectorXd& V)
{
    // 调试代码块，可以打印指定步数内的所有节点电压
    static int dbg_count = 0;
    constexpr int DBG_MAX = 0; // 设置为0以禁用
    if (dbg_count < DBG_MAX) {
        // ... 调试打印逻辑 ...
        ++dbg_count;
    }

    // 检查当前时间是否在指定的绘图范围内
    if (t >= ctrl.plot_t_start && t <= ctrl.plot_t_end) {
        // --- 写入电压文件 ---
        if (!ctrl.traces.empty() && fileV.is_open()) {
            // 如果文件头还没写，则先写入
            if (!header_written_V) {
                fileV << "Time";
                for (const auto& tr : ctrl.traces)
                    fileV << "\t" << tr.first; // 写入曲线名称
                fileV << std::endl;
                header_written_V = true;
            }
            // 写入当前时间和对应的节点电压值
            fileV << t;
            for (const auto& tr : ctrl.traces) {
                double v = 0.0;
                int node_idx = tr.second;
                if (node_idx > 0 && node_idx <= V.size())
                    v = V(node_idx - 1); // 节点编号是1-based，向量索引是0-based
                fileV << "\t" << v;
            }
            fileV << std::endl;
        }

        // --- 写入电流文件 ---
        if (!curr_traces.empty() && fileI.is_open()) {
            // 写入文件头
            if (!header_written_I) {
                fileI << "Time";
                for (const auto& ct : curr_traces)
                    fileI << "\t" << ct.name;
                fileI << std::endl;
                header_written_I = true;
            }
            // 写入当前时间和对应的电流值
            fileI << t;
            for (const auto& ct : curr_traces) {
                double i = 0.0;
                if (ct.reader)
                    i = ct.reader(); // 调用回调函数获取电流值
                fileI << "\t" << i;
            }
            fileI << std::endl;
        }
    }
}