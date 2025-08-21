#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
curve.py
默认：在脚本同目录下读取 'curve_V.dat' 和 'curve_I.dat'，先各自绘图与保存，再一次性 show() 打开两个窗口。
兼容：-i/--input 指定单个文件时仅绘该文件并弹出一个窗口。

首行是列名；若最后两个 token 为数字，则作为 Y 轴 (ymin, ymax)。
"""

from pathlib import Path
import argparse
import sys
from io import StringIO

# 依赖检查
try:
    import pandas as pd
except ImportError:
    print("缺少依赖：pandas。请先安装：pip install pandas", file=sys.stderr)
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("缺少依赖：matplotlib。请先安装：pip install matplotlib", file=sys.stderr)
    sys.exit(1)


def _is_number(token: str) -> bool:
    try:
        float(token); return True
    except Exception:
        return False


def parse_header_and_limits(file_path: Path):
    with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    header_idx = None
    for i, line in enumerate(lines):
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        header_idx = i
        header_line = s
        break

    if header_idx is None:
        raise ValueError(f"[{file_path.name}] 未找到有效的表头行。")

    tokens = header_line.split()
    if len(tokens) < 2:
        raise ValueError(f"[{file_path.name}] 表头列数不足（至少2列）。")

    y_limits = None
    if len(tokens) >= 3 and _is_number(tokens[-1]) and _is_number(tokens[-2]):
        y_limits = (float(tokens[-2]), float(tokens[-1]))
        colnames = tokens[:-2]
    else:
        colnames = tokens

    if len(colnames) < 2:
        raise ValueError(f"[{file_path.name}] 有效列不足（至少 X + 1 条 Y）。")

    data_buf = StringIO("".join(lines[header_idx + 1:]))
    return colnames, y_limits, data_buf


def read_table(file_path: Path):
    colnames, y_limits, data_buf = parse_header_and_limits(file_path)
    df = pd.read_csv(
        data_buf,
        sep=r"\s+",
        engine="python",
        comment="#",
        header=None,
        names=colnames,
    )
    df = df.dropna(axis=1, how="all")
    df.columns = [str(c).strip() for c in df.columns]
    if df.shape[1] < 2:
        raise ValueError(f"[{file_path.name}] 数据列数不足。")

    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(how="any")
    if df.empty:
        raise ValueError(f"[{file_path.name}] 清洗后为空。")

    return df, y_limits


def plot_curves(df, xcol, out_png: Path, dpi: int, title: str | None,
                y_limits=None, window_title: str | None = None):
    if xcol not in df.columns:
        raise KeyError(f"找不到X轴列 '{xcol}'。可用：{list(df.columns)}")

    x = df[xcol]
    y_cols = [c for c in df.columns if c != xcol]
    if not y_cols:
        raise ValueError("没有可用的Y轴列。")

    fig = plt.figure(figsize=(9.5, 5.5))
    if window_title:
        try:
            fig.canvas.manager.set_window_title(window_title)
        except Exception:
            pass

    for col in y_cols:
        plt.plot(x, df[col], label=col, linewidth=1.8)

    plt.xlabel(xcol)
    plt.ylabel("Value")
    if title:
        plt.title(title)
    if y_limits is not None:
        ymin, ymax = y_limits
        if ymin >= ymax:
            raise ValueError(f"非法Y轴范围：({ymin}, {ymax})")
        plt.ylim(ymin, ymax)

    plt.grid(True, which="both", linestyle="--", alpha=0.4)
    plt.legend(loc="best", frameon=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=dpi)
    return fig  # 不在此处 show/close


def plot_file(file_path: Path, out_png: Path, dpi: int, title: str | None,
              xcol_hint: str | None, window_title: str):
    df, y_limits = read_table(file_path)
    xcol = xcol_hint if xcol_hint else df.columns[0]
    fig = plot_curves(
        df, xcol,
        out_png=out_png,
        dpi=dpi,
        title=title,
        y_limits=y_limits,
        window_title=window_title,
    )
    print(f"✅ 绘图完成：{out_png.name}  ←  {file_path.name}")
    print(f"   列：{list(df.columns)}")
    if y_limits is not None:
        print(f"   使用Y轴范围：{y_limits}")
    else:
        print("   未在表头检测到Y轴范围，采用自动缩放。")
    return fig


def main():
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description="同时绘制 curve_V.dat 和 curve_I.dat（两个窗口同时打开）。兼容 -i 单文件模式。"
    )
    parser.add_argument("-i", "--input", type=Path, default=None,
                        help="单文件模式：指定一个数据文件（curve_V.dat / curve_I.dat / 旧版 curve.dat）")
    parser.add_argument("-x", "--xcol", type=str, default=None,
                        help="X轴列名；默认各文件第1列。")
    parser.add_argument("--dpi", type=int, default=150,
                        help="输出PNG的DPI（默认150）")
    parser.add_argument("--titleV", type=str, default="Voltages",
                        help="电压图标题")
    parser.add_argument("--titleI", type=str, default="Currents",
                        help="电流图标题")
    parser.add_argument("--no-show", action="store_true",
                        help="仅保存图片，不弹窗。")

    args = parser.parse_args()
    show_flag = (not args.no_show)

    figs = []

    # 单文件模式
    if args.input is not None:
        if not args.input.exists():
            print(f"找不到输入文件：{args.input}", file=sys.stderr)
            sys.exit(1)
        out_name = args.input.with_suffix(".png").name
        figs.append(
            plot_file(
                args.input,
                out_png=script_dir / out_name,
                dpi=args.dpi,
                title=None,
                xcol_hint=args.xcol,
                window_title=args.input.name,
            )
        )
    else:
        # 默认：双文件模式
        v_path = script_dir / "curve_V.dat"
        i_path = script_dir / "curve_I.dat"
        any_plotted = False

        if v_path.exists():
            figs.append(
                plot_file(
                    v_path,
                    out_png=script_dir / "curve_V.png",
                    dpi=args.dpi,
                    title=args.titleV,
                    xcol_hint=args.xcol,
                    window_title="Voltages (curve_V.dat)",
                )
            )
            any_plotted = True
        else:
            print("ℹ️ 未找到 curve_V.dat，跳过电压图。")

        if i_path.exists():
            figs.append(
                plot_file(
                    i_path,
                    out_png=script_dir / "curve_I.png",
                    dpi=args.dpi,
                    title=args.titleI,
                    xcol_hint=args.xcol,
                    window_title="Currents (curve_I.dat)",
                )
            )
            any_plotted = True
        else:
            print("ℹ️ 未找到 curve_I.dat，跳过电流图。")

        # 兜底：都没有时回退老文件
        if not any_plotted:
            legacy = script_dir / "curve.dat"
            if legacy.exists():
                print("⚠️ 未找到新文件，回退到历史 curve.dat（单窗口）")
                figs.append(
                    plot_file(
                        legacy,
                        out_png=script_dir / "curve.png",
                        dpi=args.dpi,
                        title=None,
                        xcol_hint=args.xcol,
                        window_title="curve.dat",
                    )
                )
            else:
                print("❌ 未发现可绘制的数据：curve_V.dat / curve_I.dat / curve.dat", file=sys.stderr)
                sys.exit(1)

    # 关键：一次性 show()，两个窗口会同时打开
    if show_flag and figs:
        try:
            plt.show()  # 阻塞直到窗口关闭
        except Exception as e:
            print(f"⚠️ 无法弹出窗口（可能是无显示环境）：{e}", file=sys.stderr)

    # 清理
    for fig in figs:
        try:
            plt.close(fig)
        except Exception:
            pass


if __name__ == "__main__":
    main()
