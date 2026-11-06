---
name: deeptools
license: MIT
description: 从已有 ChIP-seq、CUT&Tag 或 CUT&RUN BigWig 和 BED 快速计算 profile/heatmap，支持多进程、矩阵缓存、样本表及论文小 panel 输出。使用自定义快速实现，覆盖已验证的 deepTools mean 分箱子集；不负责 FASTQ/BAM、文库归一化或 peak calling。
---

# deeptools：快速 profile / heatmap

本 skill 名为 `deeptools`，实现位于 `scripts/fast_chip_profiles.py`，不是官方 Python `deeptools` 包。
计算对照 deepTools **3.5.1** 的 BED3–6、mean 分箱语义；不声称替代全部功能。
使用公共函数，不在项目中复制分箱、归一化或绘图实现。

## 执行优先级与输出约定

- deepTools 绘图严格遵循：**用户明确要求 > 本 skill > 项目原始实现**。当前目录的 notebook、helper、脚本或项目默认规则与本 skill 冲突时，按本 skill 修改和调用。
- 只有用户亲自明确指定的参数或样式可以覆盖对应规则；未指定的部分继续遵循本 skill。项目已有取值不等于用户要求，不因沿用项目设置而保留冲突行为。
- 默认保存图形和必要计算缓存。图形保存仅写 PDF/PNG；不额外输出绘图参数、过滤统计等说明文件，不改换成 TSV、JSON 或其他格式，也不向终端倾倒参数字典。内部核对所需信息保留在返回对象 `info` 中；缓存复用所需元信息保留为 NPZ。实际绘图数值仅在用户明确要求时另行导出。

## 输入与边界

- 样本表：`path`（BigWig）、`sample`（唯一缓存 ID）、`group`（生物学组别），可选 `label`（显示名）。每个 sample 分别绘图，不自动平均/合并重复；组别控制排列与 profile 颜色，同组重复用不同线型。默认顺序按组首次出现、组内按样本表。
- 位点输入为 `{区域集名称: BED路径}`。TE/GTF 提取、gene-level TSS 选择、consensus peak 求交留在项目 notebook；不根据文件名推断 signal normalization、Input 校正或 genome。
- BED 为 0-based 半开；负链反转、`.` 按正向。GTF 的 start 转 BED 时减 1。TSS 输入用带链完整 gene BED，锚定其 5′ 边界，不能误把已扩展 promoter BED 再当 gene。
- 不支持 BED12 exon 拼接、unscaled ends、nanAfterEnd、聚类、官方 `.mat.gz` 输出。需要时用官方 deepTools，不静默近似。

## 默认模式与横轴

区域类型与对齐模式分开定义；默认一次生成该类型所需的全部 profile 和 heatmap。

| 区域类型 | 默认输出 | 锚点与窗口 |
|---|---|---|
| TSS | 仅 center | 以链方向定义的 TSS 为中心，上下游各 2 kb；底层使用 reference-point TSS，不能使用完整 gene 的中点 |
| TE | center + body | center 为区间中点 ±3 kb；body 为本体缩放 5 kb、两侧各 3 kb |
| Peak / 其他普通区间 | center + body | center 为区间中点 ±3 kb；body 为本体缩放 5 kb、两侧各 3 kb |

- bin 默认 50 bp。区间中点为 `start + floor((end-start)/2)`，不是 peak summit。负链遵循分箱内核的反转规则，无链区间按基因组正向；TSS 必须有明确链方向。
- 横轴由计算元信息生成，profile 与 heatmap 共用：TSS 标 `−2 kb / TSS / +2 kb`；center 标 `−3 kb / TE center（或 Peak center、Region center）/ +3 kb`；body 标 `−3 kb / Start / End / +3 kb`，并注明 `Scaled body`。TE 的 Start/End 按链方向解释。
- body 的 5 kb 是缩放后的计算长度，不表示各位点真实长度；本体不能标成真实基因组距离。变更模式需计算相应矩阵，不能仅替换横轴标签。
- 可显式选择单一模式或覆盖窗口；不要按文件名猜测类型，也不要把 TSS 的展示模式 center 误传成底层 `--reference-point center`。

## 共同位点、过滤与排序

1. **同一份 BED、同一行索引。** 同一目标/区域子集的 center 与 body 从相同 BED、相同样本表开始；不得分别做 overlap 筛选、去重或信号过滤。行身份由输入文件与原始行号确定，不能仅凭可能重复的名称匹配。不同模式允许 bin 数不同，区域轴必须一致。
2. **默认保留零信号行。** 缺失值沿用本 skill 的填零偏好；两种模式均关闭独立 `skipZeros` 和 `skip_nan`。官方 deepTools 的 `--skipZeros` 默认也关闭；缺失填零则是本 skill 的选择，并非官方默认。参见 [computeMatrix 3.5.1 官方说明](https://deeptools.readthedocs.io/en/3.5.1/content/tools/computeMatrix.html)。
3. **可选联合过滤。** 仅在启用过滤时，删除 center 和 body 在全部样本、全部 bin 中均为零的位点。保留条件为 `any(center != 0) OR any(body != 0)`；不是两套非零集合的交集，也不是按 WT 或某一个样本决定。零值判断在缺失填零后的原始矩阵上进行，不用行均值为零代替逐值全零，以免正负抵消。单模式时可选过滤仅依据该模式。非法区间、缺失 contig 等输入不可用情况在计算前共同处理并单列统计，不冒充零信号过滤。
4. **共同排序。** 默认用样本表第一组第一个 sample 的 center 全窗口均值降序排序；若希望按 WT 排序，应在样本表中明确排列或指定 sort_sample。若指定 `sort_sample`，按指定 sample；不暗中平均重复。仅 body 模式才按 body 排序。并列值保持原始 BED 顺序。只计算一次行序，应用到两种模式和所有样本；不能让 body 再按自己的均值排序。
5. **共同汇总与显示。** 两种模式的 profile 与 heatmap 都应用同一最终保留掩码。Profile 的均值/SEM 包含保留的零行，不按单个 bin 再删零。超过显示行数时，两种模式使用相同行分块边界，相邻行求均值仅用于显示。
6. **内部核对。** 输入位点数、输入不可用数、可计算数、center-only 非零数、body-only 非零数、两种均非零数、两种均全零数、过滤开关、实际删除数和最终保留数用于核对，保留在内存对象中，不另写统计文件；默认不滤零时实际删除数为 0。共同保留掩码、原始行索引和排序结果供绘图复用；Broad/Narrow 等不同子集分别核对。

## 计算与缓存

- 原始 BigWig 单位与原始矩阵不改变。精确逐碱基积分，不使用 BigWig zoom 近似。body 短于 bin_size 的区间在当前内核中整行置零/缺失；其他取整空 bin 扩为 1 bp，与已验证的 deepTools 语义一致。短 body 造成的零行不应独立删掉；它不是原始轨道必然无信号的证据。
- 默认 4 个 worker，执行前确认不超过可用 CPU 数；可串行。多进程加速受 I/O 影响，不承诺固定倍数。
- 每个模式缓存完整原始矩阵；共同过滤/排序在内存中后处理，改过滤开关或排序样本不应重读 BigWig。复用核对参数、样本 ID、区域坐标与输入路径/大小/mtime；输入变化需新缓存。
- 双模式绘图需要保留完整输入行且未独立过滤的矩阵缓存；共同筛选从完整矩阵计算，不能对两个模式分别过滤后再取交集。

## 默认归一化与绘图

1. 每张热图分别取**所选完整位点 × 全部 sample × 全部 bin**的共同最大值 `M`，`H = signal / M`；全零图保留为零。所有样本使用同一分母。不是按行、按 sample 或按组单独缩放，不覆盖原始矩阵。
2. 归一化只用于 heatmap；profile 是原始单位的跨位点均值。可选 SEM 为跨位点 SEM，不是生物学重复误差。
3. 默认线性色阶下限 0，上限取 `H` 的 **第 99 百分位，向上取整为简洁刻度上限**。百分位包含零值，使用完整矩阵，在显示行合并之前计算；若为零则回退到最大值。完整归一化矩阵不裁剪，仅超出色阶的颜色饱和显示。
4. 色标为**矩形，禁止箭头/三角扩展端**；刻度从 0 开始，整数配统一倍率，例如 `0、1、2` 与 `×10⁻³`。不得把 0.001 直接标成整数 0，也不得把裁切后的上限重新标成 1。百分位与实际上限保留在返回对象 `info` 中供内部核对，不另写参数说明文件。
5. 默认 heatmap 色板固定为 `#F8FBFD → #D4E6F1 → #96C3DB → #5B9EC2 → #2C78A3 → #145474 → #10384D`。可传 Matplotlib 色板名或颜色列表；profile 默认蓝/暖橙起始，可按 group 或 sample 自定义颜色。
6. 所有 sample、center/body 使用上节规定的共同排序。每图最多显示 2,000 行，超出时按共同分块对相邻行求均值，轴标明 display bins。原始行序索引保留在返回信息中。
7. 默认 profile **60 × 45 mm**；两样本 heatmap **整体 90 × 80 mm**，多样本按列增加宽度。文字 6–7 pt、线宽约 0.5 pt，保留绘图库字体，不靠缩小字体塞入更多内容。
8. 默认 PDF + **600 dpi PNG**。PDF 的文字/线条保留矢量，heatmap 栅格层按最终物理尺寸以 600 dpi 输出；不使用 tight 裁切改变画布尺寸。具体期刊要求仍须按投稿目标核实。

`shared_max` 用于非负轨道；负信号报错而不静默丢弃，需明确用 `normalization="none"` 和合适的发散色板/范围。不同热图（包括 center 与 body）的颜色不能直接比较绝对强度，比较原始 profile；每张图分别记录归一化分母与色阶上限。

## 简洁接口

- `compute_profiles(samples, regions, cache, kind="te", modes=None, workers=4, window=None)`：返回 `{mode: cache_path}`。默认 te/peak/region 同时计算 center、body，tss 仅 center；modes 可指定单模式。每个模式缓存到 `cache/<mode>/`，保留全部输入行。
- `make_plots(caches, samples=None, regions=None, style=None, filter_all_zero=False)`：接受上述缓存映射，返回 `{mode: {区域集: {profile, heatmap, info, data}}}`。计算共同过滤与 center 排序，不写文件。regions 为原始行序的 Boolean 掩码字典，可同时生成 Broad/Narrow 等子集。
- `save_plots(plots, output, dpi=600)`：显式保存已预览的同一图对象，自动创建 `output/center/`、`body/`，每种模式仅写 profile/heatmap PDF、PNG；导出 dpi 记入该图的 `info`，不另写绘图参数文件。
- `save_profile_data(plots, caches, output, normalized=False)`：仅用户明确需要实际绘图数值时调用，保存各模式的 profile 均值/SEM/有效 n，以及共同保留掩码和行序，不附带过滤统计或参数文件。需要完整归一化矩阵时显式设 normalized=True。

样本表控制 sample/group/label；`style` 主要关注 cmap、colors、signal_label、profile_size、heatmap_size，排序 sample 可通过 sort_sample 修改。窗口仅覆盖 upstream/downstream/body_length/bin_size。filter_all_zero 默认 False；不要在项目 notebook 复制联合筛选、排序、分箱实现，也不要暴露底层 I/O 参数。

Notebook 保持计算 → 绘图 → 预览 → 保存独立，改颜色/名称不重读 BigWig。每个元件按模式分子文件夹，例如 `Figs/element_profiles/ERVK/center/`、`body/`；TSS 用 `Gene_TSS/center/`，但元信息记录 `reference_point=TSS`。缓存放对应 `results/<analysis>/<target>/<mode>/`；仅在用户要求导出实际绘图数值时，将共同保留掩码和行序置于该目标的共同父目录，不附带统计或参数说明文件；不同区域子集通过名称区分。

样本表显式传给绘图，组别/label 改动不重新计算。缓存内 sample_table 是计算时的记录，传入表可覆盖显示元数据。

API、CLI 与数值验证见 [使用与方法](references/methods-and-usage.md)。交付说明图形路径和验证范围，不主动附加参数清单或过滤统计报告。
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               