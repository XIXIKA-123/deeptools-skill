# 方法、使用与验证

> 工作流规范以 [SKILL.md](../SKILL.md) 为准：TSS 仅 center（锚定 TSS），其余类型默认 center + body；保留零行、可选联合全零过滤、共同 center 排序。公共接口使用 deepTools 3.5.1 的 mean 分箱语义。

## 与 deepTools 的对应

本实现以 deepTools 3.5.1 为对照，参考实现为该版本的 `deeptools/heatmapper.py`。
重点函数：`compute_sub_matrix_worker`、`coverage_from_big_wig`、`coverage_from_array`、`_matrix.removeempty`。
在线来源：[3.5.0 heatmapper 源码文档](https://deeptools.readthedocs.io/en/3.5.0/_modules/deeptools/heatmapper.html)、
[computeMatrix 官方说明](https://deeptools.readthedocs.io/en/3.5.1/content/tools/computeMatrix.html)。

| 步骤 | deepTools 3.5.1 | 本脚本 |
|---|---|---|
| 信号输入 | 已生成的 BigWig | 相同，不重新归一化 |
| 区域取值 | 每条 region 提取逐碱基数组 | 邻近 region 合并读取，减少重复解压与查询 |
| bin 均值 | 按 zone 用 linspace 取整边界，逐 bin masked mean | 相同边界，前缀和积分/有效碱基数 |
| 短区域 | body 小于 binSize 整行缺失/零；其余空 bin 至少取 1 bp | 同口径 |
| 负链 | 先按基因组方向分箱，再整行反转 | 同口径，避免先反向再取整产生偏移 |
| 缺失值 | 默认 NaN，可 `missingDataAsZero` | 用户偏好默认填零，可 `keep-missing` |
| skipZeros | 默认关闭；启用后对 masked matrix 的全部样本、全部 bin 求行均值，去掉无有效值或均值为零的行 | 默认关闭独立过滤，可选后处理仅删除跨模式逐值全零行 |
| 矩阵输出 | gzip 文本矩阵 | float64 NPY memmap + 独立保留行掩码 + 元信息 |
| profile | 对区域按列汇总 | mean，可选跨区域 SEM，并导出逐 bin 有效 n |
| heatmap | 可配置排序、聚类等 | 跨模式统一 center 排序，所有样本共用行序 |

**缺失不等于零。** 对区间 `[2, missing]`，保留缺失时均值为 2；填零时均值为 1。
`--missingDataAsZero --skipZeros` 并不剔除所有部分缺失区域：填零后只要整行均值非零就保留。
`--keep-missing --skipZeros` 忽略 NaN 求均值，仍保留含部分 NaN 的区域。

## 缓存

`compute -o <新目录>` 不覆盖已有目录。完整成功后写入 `metadata.npz`；中断目录缺少完整元信息，不能用作绘图缓存。

- `matrix.npy`：原始 BED 顺序，shape 为 `(全部区域数, 样本数, bin数)`，float64；包含随后被过滤的行。
- `keep.npy`：最终保留行的布尔掩码；消费矩阵时必须应用。
- `regions.tsv`：chrom/start/end/name/strand/group/kept；重复名称不作为唯一键，行号才是矩阵索引。
- `metadata.npz`：单位由调用者描述；记录参数、样本/分组、输入路径/大小/mtime、运行时长和版本。`metadata` 字段为 JSON 编码的 Unicode 标量，读取使用 `allow_pickle=False`；不生成独立 JSON 文件。

计算阶段不自动写 `profile.tsv` 等数值汇总表；绘图直接使用矩阵缓存，实际绘图数值仅在用户明确需要时由 `save_profile_data` 导出。

近似磁盘占用为 `8 × Nregion × Nsample × Nbin` bytes。默认 TE 侧翼各 3 kb、本体缩放 5 kb、bin 50 bp，共 220 bin；TSS ±2 kb 共 80 bin。
一次计算所有对比样本以固定区域集合；只改颜色、字体、排序样本或 y 轴显示不需重算。
计算参数、BigWig 或 BED 改变则生成新缓存，不自动复用过期缓存。

默认 `--chunk-bp 2000000 --batch-rows 2048 --merge-gap 10000` 控制读入跨度、坐标矩阵大小与允许跨过的无关间隙。
`-p/--processors` 使用多进程并行，按样本、染色体及最多 20,000 区域的任务块调度；仅并行计算，汇总和绘图在主进程。
每个 worker 独立打开 BigWig，写共享 memmap 的不同区域；矩阵顺序与筛选不依赖任务完成次序。
从 `-p 4` 起按任务规模测量，限制在可用 CPU 数内，并设置 `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1` 避免嵌套线程争抢。
小数据的进程启动/通信或共享存储 I/O 可能抵消收益，不能假定核数增加会线性加速。
稀疏 BED 与连续密集 BED 的加速程度不同；非常长的单一区域仍需读取整个窗口。
前缀和为 float64 累加；数值结果允许浮点舍入误差，不承诺 bitwise identical。

## 双模式联合筛选与排序规范

原始 center/body 矩阵记为 C/B，形状均为 `(Nregion, Nsample, Nbin)`，bin 数允许不同；缺失按本 skill 默认填零。须先检查原始 BED 行号、坐标、链及样本 ID/顺序逐项对应。

公共函数 paired_rows 实现以下后处理规则，并核对缓存行与输入指纹：

```python
center_nonzero = np.any(C != 0, axis=(1, 2))
body_nonzero = np.any(B != 0, axis=(1, 2))
keep = (center_nonzero | body_nonzero) if filter_all_zero else np.ones(len(C), dtype=bool)
rows = np.flatnonzero(keep)
score = C[rows, sort_sample_index, :].mean(axis=1)
order = rows[np.argsort(-score, kind="stable")]
# 两种模式用同一 keep 汇总 profile，同一 order 绘制 heatmap。
```

- `filter_all_zero` 默认 False。不能使用 `center_nonzero & body_nonzero`，也不能对两种模式分别 skipZeros；只有一种模式时可选全零过滤依据该模式。
- 输入不可用区域在计算前共同排除并记录原因；计数恒等式：输入数 = 不可用数 + 可计算数；可计算数 = center-only + body-only + 两种非零 + 两种全零；最终数 = 可计算数 − 实际过滤数。
- 共同排序按选定 center sample 的完整窗口均值计算；并列按原始 BED 行号稳定排列，不对 body 重排。若只运行 body，才使用 body 均值。超过 2,000 行时共享相同的聚合边界。
- 保留的某模式全零行必须参与该模式 profile 的均值/SEM，不能在绘图层再次删除；图中的 N 是实际参与汇总的共同位点数。
- 每张图独立做现有 shared_max 归一化和色阶计算，不将 center/body 合并取最大值。共同位点和排序不等于跨图颜色具有相同绝对尺度。
- 双模式验证覆盖：center-only 非零、body-only 非零、两种均非零、两种均全零、重复名称、排序并列、短 body、负值抵消；验证默认保留全部可计算行、可选过滤只删联合全零行、所有模式/样本行序和显示聚合一致，profile 包含保留零行。

## Python 公共接口

使用含 NumPy、pyBigWig、Matplotlib 的 Python 环境；仅绘图时无需加载 pyBigWig。
以下为分开的 notebook cell，图对象预览后再保存。

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path("skills/deeptools/scripts").resolve()))
# 以上相对路径以仓库根目录为当前目录；安装后改为实际 skill 的 scripts 目录。
from fast_chip_profiles import compute_profiles, make_plots, save_plots, save_profile_data

samples = [
    dict(path="/data/WT.bw", sample="WT", group="WT", label="WT"),
    dict(path="/data/KO.bw", sample="KO", group="KO", label="KO"),
]
```

```python
caches = compute_profiles(samples, {"ERVK": "/data/ERVK.bed"},
    "results/profiles/ERVK", kind="te", workers=4)
# 返回 {"center": Path(.../center), "body": Path(.../body)}。
# TSS: kind="tss"，仅返回 center，底层锚定 TSS；peak: kind="peak"。
# 单模式示例 modes="body"；窗口示例 window={"upstream": 2000}。
```

```python
plots = make_plots(caches, samples=samples,
    style={"colors": {"WT": "#3B6FB6", "KO": "#D98058"},
           "sort_sample": "WT", "signal_label": "ChIP (RPGC)"},
    filter_all_zero=False)
```

```python
for mode, views in plots.items():
    print(mode)
    display(views["ERVK"]["profile"])
    display(views["ERVK"]["heatmap"])
```

```python
save_plots(plots, "Figs/profiles/ERVK", dpi=600)
```

- samples 可为 records 列表或 DataFrame，sample 必须匹配全部缓存 ID。group 控制排列/颜色，label 修改显示名称；不自动合并重复。
- compute_profiles 的 regions 为 `{名称: BED路径}`；make_plots 的 regions 为 `{名称: Boolean行掩码}`，长度等于完整缓存行数。名称可以为 Broad/Narrow，子集允许重叠。
- paired_rows 校验两模式的全部输入行、样本顺序及 BED/BigWig 指纹一致；双模式只接受缺失已填零且 keep 全 True 的原始缓存。每个 view 单独计算共同保留掩码与行序；默认不做信号过滤。
- plots[mode][view]["info"] 记录归一化分母、分位数、色阶、原始 row_indices、display_edges、共同过滤数量、排序来源和物理尺寸。data 保存该图实际使用的 mean/sem/n_valid。
- save_plots 每种模式仅保存 profile/heatmap PDF、PNG，并将导出 dpi 记入图对象的 info；不生成 plot_settings.json，也不另写同类参数文件。
- save_profile_data 仅在用户明确需要实际绘图数值时调用，写共同 `<view>_keep.npy`、`<view>_rows.tsv`，以及每种模式的 `<view>_profile.tsv`；不附带过滤统计或参数文件。normalized=True 才额外写完整归一化矩阵，行序对应共同 rows.tsv。输入不可用与区域选择数量仅作内部核对，不另写统计文件。
- style 主要使用 cmap、colors、signal_label、profile_size、heatmap_size；可选 quantile=.99、sort_sample、max_rows=2000、sem=False、title。颜色优先按 sample 匹配，再按 group。
- 非负轨道默认 shared_max；负信号须显式设 normalization="none" 与发散色板，例如 cmap="RdBu_r"、clim=(-2,2)。

### 归一化与整数色标例子

`M = max(本图全部所选位点、样本和 bin)`，`H = X / M`；M=0 时保留零。NaN 保持为缺失、统计时忽略，默认计算阶段已经填零。
取完整 H 的第 99 百分位（包含零）作为显示上限候选，再向上取整到简洁整数刻度：例如 0.001604 → 0.002，色标标 0、1、2 和 ×10⁻³。
矩阵 H 不改；0.002 以上的颜色显示为最深色，色标端部仍是矩形。分位为零时回退到最大值，全零图显示 0–1。
`quantile=1` 使用完整范围候选；非负共享最大值归一化后最大值为 1。Profile 始终使用原始单位。

### CLI

CLI 默认保留零行；显式 --skipZeros 使用 deepTools 的单模式过滤语义。以下命令输出单模式缓存，不完成跨模式共同排序；新分析优先使用上面的 Python 公共接口。

```bash
SCRIPT=skills/deeptools/scripts/fast_chip_profiles.py
python "$SCRIPT" compute -S WT.bw KO.bw -R ERVK.bed \
  --sample-labels WT KO --group-labels ERVK --mode scale-regions \
  --missingDataAsZero --keep-zero-rows -p 4 -o results/profiles/ERVK/matrix
python "$SCRIPT" plot --matrix results/profiles/ERVK/matrix \
  -o Figs/profiles/ERVK --signal-label 'ChIP (RPGC)'
```

不指定侧翼时 CLI 的 scale-regions 默认 3 kb，reference-point 默认 2 kb。计算 TE/peak center 时必须显式设 `--mode reference-point --reference-point center -b 3000 -a 3000`；TSS 则用 `--reference-point TSS -b 2000 -a 2000`。
CLI 保留底层计算控制；自定义样本组别、显示名和样式建议用上面的 Python 接口。
`plot --quantile` 改分位数，`--normalization none` 保留原始值；`--vmin`/`--vmax` 必须一起提供。
Python 多模式输出按 center/body 分层；CLI 输出单模式结果。

## 对照测试

随包 `scripts/validate_fast_chip_profiles.py` 用真实的已安装 deepTools CLI 作 oracle，覆盖：
正负链、奇数长度、1–99 bp 短区域、边界截断、BigWig gap、显式零、两样本联合过滤，
scale-regions、TSS/TES/center，填零/保留缺失与 skipZeros 的组合，以及缓存、profile 汇总与实际 PDF/PNG 渲染。
验证运行只写调用者指定的临时工作目录，不改正式结果；对照检查结果返回终端，不另写统计报告文件。

```bash
RUN_TMP=$(mktemp -d)
PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
MPLCONFIGDIR="$RUN_TMP/mpl" \
python skills/deeptools/scripts/validate_fast_chip_profiles.py \
  --workdir "$RUN_TMP/validation"
```

可额外传入 `--benchmark-bed`、`--benchmark-bw` 比较真实数据。
报告区分本脚本 core 计算时间和 deepTools CLI 总时间；公平端到端测速需另外运行本脚本 CLI，
包括解释器启动、读 BED 和写缓存，并在相同机器、CPU 数、文件和参数下比较。
过滤后区域名称、NaN 位置必须一致；deepTools 导出文本为 6 位小数，比较容差为 `atol=5.1e-7, rtol=1e-7`。
底层未舍入数值另以 `1e-10` 容差对照。不能仅看曲线相似就宣称兼容。
