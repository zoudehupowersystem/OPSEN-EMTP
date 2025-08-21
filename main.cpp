// main.cpp
#include "Control.h"
#include "Curve.h"
#include "Devices.h"
#include "Grid.h"
#include "Simulation.h"
#include <iostream>
#include <memory>
#include <vector>

int main()
{
    std::cout << "--- EMTP Simulation for  Simple Case ---" << std::endl;

    Control control;
    control.setupSimpleCase();

    // --- System and Line Parameters from Paper Case 4.1 ---
    const double freq = 50.0;
    const double v_line = 110.0; // 线电压(kV)
    const double line_length_km = 5.0;
    const double R1_per_km = 0.2, R0_per_km = 0.6;
    const double X1_per_km = 0.26, X0_per_km = 1.13;
    const double C1_per_km = 0.012e-6, C0_per_km = 0.006e-6;

    // Total line parameters
    const double R1 = R1_per_km * line_length_km;
    const double R0 = R0_per_km * line_length_km;
    const double L1 = (X1_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double L0 = (X0_per_km * line_length_km) / (2 * 3.14159 * freq);
    const double C1 = C1_per_km * line_length_km;
    const double C0 = C0_per_km * line_length_km;

    // 0:GND, 1-3:Bus1(A,B,C), 4-6:Bus2(A,B,C), 7-9:中间节点(A,B,C)
    Grid grid(9);

    // 1) 源仍接在 1/2/3 对地
    double v_peak = v_line * sqrt(2.0 / 3.0);
    double rs = 0.001, isc_peak = 48;
    auto srcA = std::make_unique<VoltageSource>(1, 0, v_peak, freq, 0.0, rs, isc_peak);
    auto srcB = std::make_unique<VoltageSource>(2, 0, v_peak, freq, -120.0, rs, isc_peak);
    auto srcC = std::make_unique<VoltageSource>(3, 0, v_peak, freq, +120.0, rs, isc_peak);
    auto srcA_ptr = srcA.get(), srcB_ptr = srcB.get(), srcC_ptr = srcC.get();
    grid.addDevice(std::move(srcA));
    grid.addDevice(std::move(srcB));
    grid.addDevice(std::move(srcC));

    // 2) 断路器：直接用固定的 7/8/9 做中间节点
    const double R_closed = 1e-5, R_open = 1e9;
    auto brkA = std::make_unique<CircuitBreakerPhase>(1, 7, R_closed, R_open, /*initiallyClosed*/ true, freq);
    auto brkB = std::make_unique<CircuitBreakerPhase>(2, 8, R_closed, R_open, /*initiallyClosed*/ true, freq);
    auto brkC = std::make_unique<CircuitBreakerPhase>(3, 9, R_closed, R_open, /*initiallyClosed*/ true, freq);
    auto brkA_ptr = brkA.get(), brkB_ptr = brkB.get(), brkC_ptr = brkC.get();
    grid.addDevice(std::move(brkA));
    grid.addDevice(std::move(brkB));
    grid.addDevice(std::move(brkC));

    // 3) 线路：从 7/8/9 → 4/5/6
    auto line = std::make_unique<PI_line>(
        std::vector<int> { 7, 8, 9 },
        std::vector<int> { 4, 5, 6 },
        R1, L1, C1, R0, L0, C0, freq, control.dt);
    grid.addDevice(std::move(line));

    // 4) 负荷
    auto load = std::make_unique<Load>(std::vector<int> { 4, 5, 6 }, 60, 0, v_line, freq, 0);
    grid.addDevice(std::move(load));

    // 5)故障
    const int fault_node = 4; // Bus2 A 相
    auto fault_switch = std::make_unique<TwoValueResistor>(fault_node, 0, 1e-3, 1e9);
    TwoValueResistor* fault_switch_ptr = fault_switch.get();
    grid.addDevice(std::move(fault_switch));

    // 6) 曲线与仿真器
    Curve curve(control, fault_switch_ptr);
    curve.addCurrentTrace("I_src_A", [srcA_ptr]() { return srcA_ptr->get_I(); });
    curve.addCurrentTrace("I_src_B", [srcB_ptr]() { return srcB_ptr->get_I(); });
    curve.addCurrentTrace("I_src_C", [srcC_ptr]() { return srcC_ptr->get_I(); });

    Simulation simulation(control, grid, curve, fault_switch_ptr, fault_node);
    simulation.addBreaker(brkA_ptr);
    simulation.addBreaker(brkB_ptr);
    simulation.addBreaker(brkC_ptr);

    // 7) 动作时序（示例：故障→A相跳→1s后重合→0.1s后三相跳）
    const double t_trip_A = 2.10;
    const double t_recloseA = 3.10;
    const double t_trip_3ph = 3.20;

    brkA_ptr->schedule(t_trip_A, CircuitBreakerPhase::Change::CloseToOpen);
    brkA_ptr->schedule(t_recloseA, CircuitBreakerPhase::Change::OpenToClose);
    brkA_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);
    brkB_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);
    brkC_ptr->schedule(t_trip_3ph, CircuitBreakerPhase::Change::CloseToOpen);

    simulation.run();

    return 0;
}