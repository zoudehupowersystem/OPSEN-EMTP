#include "Control.h"

void Control::setupSimpleCase()
{
    // --- 定义事件 ---
    // 在 2.0 秒时，在节点 4 施加一个永久性故障
    events.push_back(FaultEvent { 2.0, 4, true });

    // --- 定义需要记录的电压波形 ---
    traces.clear();
    traces.emplace_back("voltage_A", 4); // 记录节点4的电压，命名为 "voltage_A"
    traces.emplace_back("voltage_B", 5); // 记录节点5的电压，命名为 "voltage_B"
    traces.emplace_back("voltage_C", 6); // 记录节点6的电压，命名为 "voltage_C"
}