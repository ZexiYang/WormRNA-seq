# patches — 上游源码关键补丁

重上环境（重新安装 zUMIs / umite）后必须重新应用，否则流程会在特定数据上崩溃。
均为 unified diff，`patch < 对应文件 < xxx.patch` 或按内容手动修改。

## fqfilter_v2.pl.patch（zUMIs 2.9.7e）

**问题**：Smart-seq3 bulk 配置（yaml 无 `BC()` 字段）时，所有 reads 的 BC tag 为空 →
BCstats/kept_barcodes 的 XC 为空 → R `fread` 读成 NA → `ScanBamParam tagFilter`
拒收 NA → Counting 阶段全部 worker 报错（`rbindlist: Item 1 is not a data.frame`）。

**补丁**：`fqfilter_v2.pl:268` 附近，空 bcseq 一律赋常量 `CELL1`。

配套要求：yaml 必须 `BarcodeBinning: 0`（单 barcode + binning=1 会让 BCbin 的
setnames 崩溃）。

## zUMIs-BCdetection.R.patch（zUMIs 2.9.7e）

**问题**：`BarcodeBinning: 0` 时不写 `kept_barcodes_binned.txt`，但 Counting 阶段
无条件读取该文件 → 报错。

**补丁**：binning 关闭时也写出 binned 文件。

## umiextract.py.patch（umite 0.1.1）

**问题**：umiextract 对 R1/R2 读名做全等校验；华大 DNBSEQ 读名带 `/1` `/2` 后缀
（如 `E250105976L1C001R00200000053/1`），直接报 `readname mismatch`。
且 umicount 靠 read name 配对，后缀必须去除。

**补丁**：`partition(' ')` 之后剥除末尾 `/1` `/2`。

另注意（非补丁）：`umiextract -d/--output_dir` 只认**已存在的绝对路径**，
相对路径会误报 `Folder does not exist`；`pipelines/umite/scripts/run_umite.sh`
已内置处理。
