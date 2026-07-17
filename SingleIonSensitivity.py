#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
逐帧投影 + 峰面积与本底分析 + 逐帧绘图。

对 TH2F 的每一帧 (Y bin) 做 X 轴投影 → TH1，
在 309.732~309.736 MHz 范围内找峰，
以峰中心为中心、df 宽度计算峰面积，
以峰中心右边 delta_f 处、df 宽度计算本底面积，
在时间区间 [t1, t2] 内统计平均值和误差。
对每一帧的投影谱画图，标注峰/本底积分区间和面积。

用法（WSL conda）：
    conda activate base
    python3 SingleIonSensitivity.py --t1 0 --t2 2.9
    python3 SingleIonSensitivity.py --t1 0 --t2 1.35
"""

import argparse
import sys
from array import array
from pathlib import Path

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kError
gStyle = ROOT.gStyle
gStyle.SetPalette(51)
gStyle.SetNumberContours(255)
gStyle.SetOptStat(0)


def find_peak_in_range(h1, x_low, x_high):
    """在 [x_low, x_high] MHz 范围内寻找最高峰, 返回峰中心 (MHz)"""
    xaxis = h1.GetXaxis()
    bin_low = xaxis.FindBin(x_low)
    bin_high = xaxis.FindBin(x_high)
    if bin_low < 1:
        bin_low = 1
    if bin_high > h1.GetNbinsX():
        bin_high = h1.GetNbinsX()

    max_bin = bin_low
    max_val = h1.GetBinContent(bin_low)
    for b in range(bin_low, bin_high + 1):
        val = h1.GetBinContent(b)
        if val > max_val:
            max_val = val
            max_bin = b
    return xaxis.GetBinCenter(max_bin)


def calc_area(h1, center_MHz, df_MHz):
    """
    以 center_MHz 为中心、半宽 df_MHz/2，计算峰面积（bin 含量之和）。
    误差用积分区间内各 bin 含量的标准差。
    返回 (area, area_err)
    """
    xaxis = h1.GetXaxis()
    bin_low = xaxis.FindBin(center_MHz - df_MHz / 2.0)
    bin_high = xaxis.FindBin(center_MHz + df_MHz / 2.0)
    if bin_low < 1:
        bin_low = 1
    if bin_high > h1.GetNbinsX():
        bin_high = h1.GetNbinsX()

    vals = []
    total = 0.0
    for b in range(bin_low, bin_high + 1):
        val = h1.GetBinContent(b)
        total += val
        vals.append(val)

    if len(vals) <= 1:
        return total, 0.0

    # 用积分区间内 bin 含量的标准差作为面积误差
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
    std = var ** 0.5
    return total, std


def add_boundary_lines(h1, center, half_w, color, line_style=7):
    """
    在 center±half_w 处画两条虚线竖线，标记积分区间的频率边界。
    """
    x_low = center - half_w
    x_high = center + half_w
    y_max_frame = h1.GetMaximum() * 1.0

    line_low = ROOT.TLine(x_low, 0.0, x_low, y_max_frame)
    line_low.SetLineColor(color)
    line_low.SetLineWidth(2)
    line_low.SetLineStyle(line_style)
    line_low.DrawClone()

    line_high = ROOT.TLine(x_high, 0.0, x_high, y_max_frame)
    line_high.SetLineColor(color)
    line_high.SetLineWidth(2)
    line_high.SetLineStyle(line_style)
    line_high.DrawClone()


def get_region_edges(center, half_w):
    """返回积分区间 [center-half_w, center+half_w] 的频率边界 (MHz)"""
    return center - half_w, center + half_w


def draw_frame(h1, iy, time_center, peak_center, bg_center, df_MHz,
               p_area, p_err, b_area, b_err, net, net_err, output_dir, root_stem):
    """
    画出单帧投影谱，标注峰/本底积分区间和面积值，保存 PNG。
    """
    # X 轴显示范围：包含了峰和本底区域，留一些边距
    half_span = (bg_center - peak_center) + df_MHz + 0.001
    x_low = peak_center - half_span
    x_high = bg_center + half_span

    c = ROOT.TCanvas(f"c_frame_{iy}", f"Frame {iy} t={time_center:.4f}s", 1000, 600)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.13)
    c.SetTopMargin(0.10)

    # 设置直方图样式
    h1.SetLineColor(ROOT.kBlue)
    h1.SetLineWidth(2)
    h1.SetTitle(f"Frame Y Bin {iy}  t={time_center:.4f} s;Frequency [MHz];Counts")
    h1.GetXaxis().SetRangeUser(x_low, x_high)
    h1.GetXaxis().SetTitleSize(0.05)
    h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.SetStats(0)

    # 获取 y 范围供稍后画线和文字用
    y_min = h1.GetMinimum()
    y_max = h1.GetMaximum()
    if y_max <= 0:
        y_max = 1
    y_range = y_max - y_min
    h1.GetYaxis().SetRangeUser(y_min - 0.1 * y_range, y_max + 0.3 * y_range)

    h1.Draw("HIST")

    # ---- 峰积分区间的频率边界（绿色虚线）----
    add_boundary_lines(h1, peak_center, df_MHz / 2.0, ROOT.kGreen + 2)

    # ---- 本底积分区间的频率边界（红色虚线）----
    add_boundary_lines(h1, bg_center, df_MHz / 2.0, ROOT.kRed)

    # ---- 峰中心位置（粉色虚线）----
    line_center = ROOT.TLine(peak_center, 0, peak_center, y_max)
    line_center.SetLineColor(ROOT.kMagenta + 1)
    line_center.SetLineWidth(1)
    line_center.SetLineStyle(3)
    line_center.DrawClone()

    # ---- 标注文字（右下 TPaveText）----
    txt = ROOT.TPaveText(0.55, 0.50, 0.88, 0.88, "NDC")
    txt.SetFillColor(0)
    txt.SetTextSize(0.035)
    txt.SetBorderSize(1)
    txt.AddText(f"Peak: {p_area:.5f} #pm {p_err:.5f}")
    txt.AddText(f"BG:   {b_area:.5f} #pm {b_err:.5f}")
    txt.AddText(f"Net:  {net:.5f} #pm {net_err:.5f}")
    txt.AddText(f"t = {time_center:.4f} s")
    txt.DrawClone()

    c.Update()
    png_path = output_dir / f"frame_{iy:03d}_t{time_center:.4f}s_{root_stem}.png"
    c.SaveAs(str(png_path))
    c.Close()
    print(f"  -> {png_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description="逐帧投影 + 峰面积和本底分析 + 逐帧绘图")
    parser.add_argument("--t1", type=float, required=True,
                        help="时间区间下限 (s)")
    parser.add_argument("--t2", type=float, required=True,
                        help="时间区间上限 (s)")
    parser.add_argument("--root-file", type=str, default=None,
                        help="ROOT 文件路径（默认 data 目录下的指定文件）")
    parser.add_argument("--hist-name", type=str, default="h2d",
                        help="TH2F 名称（默认 h2d）")
    parser.add_argument("--search-low", type=float, default=309.732,
                        help="找峰范围下限 (MHz)，默认 309.732")
    parser.add_argument("--search-high", type=float, default=309.736,
                        help="找峰范围上限 (MHz)，默认 309.736")
    parser.add_argument("--df", type=float, default=0.001,
                        help="峰/本底积分宽度 (MHz)，默认 0.001 (1 kHz)")
    parser.add_argument("--delta-f", dest="delta_f", type=float, default=0.003,
                        help="峰中心到本底区域的偏移 (MHz)，默认 0.003 (3 kHz)")
    parser.add_argument("--plot-dir", type=str, default="frame_plots",
                        help="帧图片输出子目录名（默认 frame_plots）")
    args = parser.parse_args()

    # ====================== 文件路径 ======================
    if args.root_file:
        root_file = Path(args.root_file)
    else:
        script_dir = Path(__file__).parent
        root_file = (script_dir / "203Pt78_decay_BGSub" / "data"
                     / "8243_PY82ch1_0242_trigger_11_2026-04-08T01-41-32.root")

    if not root_file.exists():
        print(f"错误：找不到 ROOT 文件: {root_file}")
        sys.exit(1)

    print(f"读取 ROOT 文件: {root_file}")

    # ====================== 读取 TH2F ======================
    f = ROOT.TFile.Open(str(root_file))
    if not f or f.IsZombie():
        print("错误：无法打开 ROOT 文件")
        sys.exit(1)

    h2 = f.Get(args.hist_name)
    if not h2 or not isinstance(h2, ROOT.TH2):
        print(f"错误：找不到 TH2 '{args.hist_name}'")
        f.Close()
        sys.exit(1)

    n_y = h2.GetNbinsY()
    yaxis = h2.GetYaxis()
    print(f"\nTH2F: {h2.GetName()}")
    print(f"  X轴: {h2.GetNbinsX()} bins, [{h2.GetXaxis().GetXmin():.3f}, {h2.GetXaxis().GetXmax():.3f}] MHz")
    print(f"  Y轴: {n_y} bins, [{yaxis.GetXmin():.4f}, {yaxis.GetXmax():.4f}] s")

    # ====================== 确定时间 bin 区间 ======================
    bin_t1 = yaxis.FindBin(args.t1)
    bin_t2 = yaxis.FindBin(args.t2)
    if bin_t1 < 1:
        bin_t1 = 1
    if bin_t2 > n_y:
        bin_t2 = n_y
    if bin_t1 > bin_t2:
        print("错误：t1 应小于 t2")
        f.Close()
        sys.exit(1)

    print(f"\n时间区间: [{args.t1}, {args.t2}] s -> Y bin [{bin_t1}, {bin_t2}]")

    # ====================== 先用平均谱找峰 ======================
    h1_sum = h2.ProjectionX("h1_sum_for_peak", bin_t1, bin_t2)
    peak_center = find_peak_in_range(h1_sum, args.search_low, args.search_high)
    print(f"\n在 [{args.search_low}, {args.search_high}] MHz 找到峰中心: {peak_center:.6f} MHz")

    bg_center = peak_center + args.delta_f
    print(f"本底区域中心: {bg_center:.6f} MHz")

    df_MHz = args.df
    print(f"积分宽度 df: {df_MHz*1e3:.1f} kHz = {df_MHz:.6f} MHz")

    # ====================== 创建图片输出目录 ======================
    plot_dir = Path(__file__).parent / args.plot_dir
    plot_dir.mkdir(exist_ok=True)
    print(f"\n帧图片保存目录: {plot_dir}/")

    # ====================== 逐帧处理 + 绘图 ======================
    results = []

    for iy in range(bin_t1, bin_t2 + 1):
        time_center = yaxis.GetBinCenter(iy)

        h1 = h2.ProjectionX(f"h1_frame_{iy}", iy, iy)

        p_area, p_err = calc_area(h1, peak_center, df_MHz)
        b_area, b_err = calc_area(h1, bg_center, df_MHz)
        net = p_area - b_area
        net_err = (p_err ** 2 + b_err ** 2) ** 0.5

        # 峰/本底积分区间的频率边界
        p_lo, p_hi = get_region_edges(peak_center, df_MHz / 2.0)
        b_lo, b_hi = get_region_edges(bg_center, df_MHz / 2.0)

        # 相对误差
        p_rel_err = p_err / p_area if p_area != 0 else 0.0
        b_rel_err = b_err / b_area if b_area != 0 else 0.0
        n_rel_err = net_err / net if net != 0 else 0.0

        results.append((time_center, p_area, p_err, b_area, b_err, net, net_err,
                        p_lo, p_hi, b_lo, b_hi,
                        p_rel_err, b_rel_err, n_rel_err))

        # ---- 画图 ----
        draw_frame(h1.Clone(f"h1_frame_{iy}_plot"), iy, time_center,
                   peak_center, bg_center, df_MHz,
                   p_area, p_err, b_area, b_err, net, net_err,
                   plot_dir, root_file.stem)

    # ====================== 统计平均 ======================
    n_frames = len(results)
    if n_frames == 0:
        print("错误：没有符合条件的帧")
        f.Close()
        sys.exit(1)

    def avg_and_err(vals):
        n = len(vals)
        if n == 0:
            return 0, 0, 0
        mean = sum(vals) / n
        if n == 1:
            return mean, 0, 0
        var = sum((v - mean) ** 2 for v in vals) / (n - 1)
        std = var ** 0.5
        return mean, std / (n ** 0.5), std

    p_vals = [r[1] for r in results]
    b_vals = [r[3] for r in results]
    n_vals = [r[5] for r in results]
    p_rel_vals = [r[11] for r in results]
    b_rel_vals = [r[12] for r in results]
    n_rel_vals = [r[13] for r in results]

    p_mean, p_mean_err, p_std = avg_and_err(p_vals)
    b_mean, b_mean_err, b_std = avg_and_err(b_vals)
    n_mean, n_mean_err, n_std = avg_and_err(n_vals)
    p_rel_mean, p_rel_mean_err, p_rel_std = avg_and_err(p_rel_vals)
    b_rel_mean, b_rel_mean_err, b_rel_std = avg_and_err(b_rel_vals)
    n_rel_mean, n_rel_mean_err, n_rel_std = avg_and_err(n_rel_vals)

    # ====================== 打印结果 ======================
    print(f"\n{'='*160}")
    print(f"逐帧分析结果: t=[{args.t1}, {args.t2}] s, 共 {n_frames} 帧")
    print(f"峰中心: {peak_center:.6f} MHz,  本底中心: {bg_center:.6f} MHz,  df: {df_MHz*1e3:.1f} kHz")
    print(f"{'='*160}")
    print(f"{'帧中心(s)':>10s} {'峰面积':>8s} {'峰误差':>8s} {'本底面积':>8s} {'本底误差':>8s} {'净面积':>8s} {'净误差':>8s} "
          f"{'峰低边':>10s} {'峰高边':>10s} {'本底低边':>10s} {'本底高边':>10s} "
          f"{'峰相对误差':>10s} {'本底相对误差':>10s} {'净相对误差':>10s}")
    print(f"{'-'*160}")
    for r in results:
        print(f"{r[0]:>10.4f} {r[1]:>8.5f} {r[2]:>8.5f} {r[3]:>8.5f} {r[4]:>8.5f} {r[5]:>8.5f} {r[6]:>8.5f} "
              f"{r[7]:>10.6f} {r[8]:>10.6f} {r[9]:>10.6f} {r[10]:>10.6f} "
              f"{r[11]:>10.3f} {r[12]:>10.3f} {r[13]:>10.3f}")

    print(f"{'─'*160}")
    print(f"{'平均':>10s} {p_mean:>8.5f} {p_mean_err:>8.5f} {b_mean:>8.5f} {b_mean_err:>8.5f} {n_mean:>8.5f} {n_mean_err:>8.5f}")
    print(f"{'标准差':>10s} {p_std:>8.5f} {'':>8s} {b_std:>8.5f} {'':>8s} {n_std:>8.5f}")
    print(f"{'平均相对误差':>10s} {p_rel_mean:>8.3f} {p_rel_mean_err:>8.3f} {b_rel_mean:>8.3f} {b_rel_mean_err:>8.3f} {n_rel_mean:>8.3f} {n_rel_mean_err:>8.3f}")
    print(f"{'='*150}")

    print(f"\n# ONE-LINE RESULT: t1={args.t1} t2={args.t2} "
          f"peak_center={peak_center:.6f} bg_center={bg_center:.6f} "
          f"df_kHz={df_MHz*1e3:.1f} "
          f"peak_mean={p_mean:.5f} peak_mean_err={p_mean_err:.5f} "
          f"bg_mean={b_mean:.5f} bg_mean_err={b_mean_err:.5f} "
          f"net_mean={n_mean:.5f} net_mean_err={n_mean_err:.5f} "
          f"n_frames={n_frames}")

    # ====================== 保存 CSV ======================
    csv_path = script_dir / f"results_t{args.t1:.2f}_{args.t2:.2f}s_{root_file.stem}.csv"
    with open(csv_path, 'w', encoding='utf-8') as csv_f:
        csv_f.write("time_s,peak_area,peak_err,bg_area,bg_err,net_area,net_err,"
                    "peak_low_MHz,peak_high_MHz,bg_low_MHz,bg_high_MHz,"
                    "peak_rel_err,bg_rel_err,net_rel_err\n")
        for r in results:
            csv_f.write(f"{r[0]:.6f},{r[1]:.8f},{r[2]:.8f},{r[3]:.8f},{r[4]:.8f},{r[5]:.8f},{r[6]:.8f},"
                        f"{r[7]:.8f},{r[8]:.8f},{r[9]:.8f},{r[10]:.8f},"
                        f"{r[11]:.6f},{r[12]:.6f},{r[13]:.6f}\n")
        csv_f.write(f"\n# summary: t1={args.t1},t2={args.t2},peak_center={peak_center:.6f},bg_center={bg_center:.6f},df_kHz={df_MHz*1e3:.1f}\n")
        csv_f.write(f"# mean,{p_mean:.8f},{p_mean_err:.8f},{b_mean:.8f},{b_mean_err:.8f},{n_mean:.8f},{n_mean_err:.8f}\n")
        csv_f.write(f"# std,{p_std:.8f},,{b_std:.8f},,{n_std:.8f},\n")
    print(f"\nCSV 已保存: {csv_path}")

    # ====================== 绘制面积随时间变化图 ======================
    n = n_frames
    t_arr = array('d', [r[0] for r in results])
    p_y = array('d', [r[1] for r in results])
    p_ey = array('d', [r[2] for r in results])
    b_y = array('d', [r[3] for r in results])
    b_ey = array('d', [r[4] for r in results])
    n_y_vals = array('d', [r[5] for r in results])
    n_ey = array('d', [r[6] for r in results])

    ex_zero = array('d', [0.0] * n)

    gr_peak = ROOT.TGraphErrors(n, t_arr, p_y, ex_zero, p_ey)
    gr_peak.SetName("gr_peak_area")
    gr_peak.SetTitle("Peak Area")
    gr_peak.SetLineColor(ROOT.kGreen + 2)
    gr_peak.SetLineWidth(2)
    gr_peak.SetLineStyle(2)       # 宽虚线
    gr_peak.SetMarkerColor(ROOT.kGreen + 2)
    gr_peak.SetMarkerStyle(24)    # 空心大圆圈
    gr_peak.SetMarkerSize(1.5)

    gr_bg = ROOT.TGraphErrors(n, t_arr, b_y, ex_zero, b_ey)
    gr_bg.SetName("gr_bg_area")
    gr_bg.SetTitle("BG Area")
    gr_bg.SetLineColor(ROOT.kRed)
    gr_bg.SetLineWidth(2)
    gr_bg.SetMarkerColor(ROOT.kRed)
    gr_bg.SetMarkerStyle(21)
    gr_bg.SetMarkerSize(0.8)

    gr_net = ROOT.TGraphErrors(n, t_arr, n_y_vals, ex_zero, n_ey)
    gr_net.SetName("gr_net_area")
    gr_net.SetTitle("Net Area")
    gr_net.SetLineColor(ROOT.kBlue)
    gr_net.SetLineWidth(2)
    gr_net.SetMarkerColor(ROOT.kBlue)
    gr_net.SetMarkerStyle(22)
    gr_net.SetMarkerSize(0.8)

    # 找出统一的 Y 范围
    all_y = list(p_y) + list(b_y) + list(n_y_vals)
    all_ey = list(p_ey) + list(b_ey) + list(n_ey)
    y_min_all = min(all_y) - max(all_ey) * 1.2
    y_max_all = max(all_y) + max(all_ey) * 1.2
    if y_min_all >= 0:
        y_min_all = -max(all_y) * 0.05

    c_evo = ROOT.TCanvas("c_evolution", "Area vs Time", 1200, 700)
    c_evo.SetLeftMargin(0.12)
    c_evo.SetRightMargin(0.05)
    c_evo.SetBottomMargin(0.13)
    c_evo.SetTopMargin(0.12)

    frame = c_evo.DrawFrame(t_arr[0], y_min_all, t_arr[-1], y_max_all,
                            f"Area Evolution  t=[{args.t1}, {args.t2}] s;Time [s];Area")
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)

    gr_peak.Draw("LP SAME")
    gr_bg.Draw("P SAME")
    gr_net.Draw("LP SAME")

    legend = ROOT.TLegend(0.15, 0.70, 0.38, 0.88)
    legend.SetFillColor(0)
    legend.SetBorderSize(1)
    legend.AddEntry(gr_peak, "Peak", "pl")
    legend.AddEntry(gr_bg, "BG", "pl")
    legend.AddEntry(gr_net, "Net", "pl")
    legend.DrawClone()

    c_evo.Modified()
    c_evo.Update()

    # ---- 保存 PNG ----
    evo_png = script_dir / f"area_evolution_t{args.t1:.2f}_{args.t2:.2f}s_{root_file.stem}.png"
    c_evo.SaveAs(str(evo_png))
    print(f"\n面积演化图已保存: {evo_png}")

    # ---- 保存 ROOT 文件 ----
    output_root = script_dir / f"area_evolution_t{args.t1:.2f}_{args.t2:.2f}s_{root_file.stem}.root"
    out_root = ROOT.TFile(str(output_root), "RECREATE")
    gr_peak.Write()
    gr_bg.Write()
    gr_net.Write()
    c_evo.Write("c_evolution")
    out_root.Close()
    print(f"ROOT 文件已保存: {output_root}")

    f.Close()


if __name__ == "__main__":
    main()
