# zUMIs (Smart-seq3) 测试样本提交指南

> 目标：学会在 zhllab 集群上提交一个 Slurm 作业，跑通一个样本的 zUMIs UMI 定量。
> 所有路径以本项目为例，换成你自己的目录即可复用。

---

## 0. 背景知识（30 秒版）

- 集群登录节点（你 ssh 上去的机器）**只用来提交和管理作业**，不能跑计算。
- 计算由 **Slurm** 调度到计算节点（cpu01–08 / fat01–02 / gpu01）执行。
- 你写一个 `.sbatch` 脚本（描述"我要多少资源 + 跑什么命令"），用 `sbatch` 提交，
  Slurm 排队、分配节点、执行，日志写到文件里。

## 1. 脚本在哪、做什么

```
/share/home/zhllab_student/yzx/260727_singleworm/scripts/
├── zumis-test-cpu02.sbatch   ← 本次要提交的测试脚本（单样本, cpu02）
├── zumis-batch-fat.sbatch    ← 测试通过后用的全量批量脚本（785 样本, fat 节点）
└── samples.txt               ← 785 个样本名清单（批量版用）
```

测试脚本做三件事：

1. 检查输入文件（原始 fq、STAR 索引、GTF）是否存在
2. 为该样本生成 zUMIs 配置文件 `06zumis/yamls/<sample>.zUMIs.yaml`
   （Smart-seq3 配置：R1 的 12–19 位是 UMI，识别 `ATTGCGCAATG` pattern）
3. 调用 `~/src/zUMIs/zUMIs.sh -c -y <yaml>` 运行，输出到 `06zumis/<sample>/`

## 2. sbatch 脚本头逐行解释

```bash
#!/bin/bash                       # 用 bash 解释执行
#SBATCH --job-name=zumis-test     # 作业名（squeue 里显示）
#SBATCH --partition=compute_cpu   # 提交到哪个分区（队列）
#SBATCH --nodelist=cpu02          # 指定在 cpu02 节点上跑（可省略，让调度器自选）
#SBATCH --time=1-00:00:00         # 最长运行时间 1 天（超时会被杀掉）
#SBATCH --cpus-per-task=32        # 申请 32 个 CPU 核
#SBATCH --mem=90G                 # 申请 90G 内存（cpu02 共 100G）
#SBATCH --output=.../zumis-test.%j.log   # 标准输出+错误写入的日志文件
                                  # %j 会被替换成实际作业号
```

资源申请原则：**够用就好**。申请越多排队越久；cpu02 是 56 核/100G，
这个测试申请 32 核/90G，同节点还能给别人留点余量。

## 3. 提交（在 zhllab 上执行）

```bash
ssh zhllab
sbatch /share/home/zhllab_student/yzx/260727_singleworm/scripts/zumis-test-cpu02.sbatch day11_CF_4
```

- 脚本路径后面的 `day11_CF_4` 是传给脚本的参数（脚本里用 `$1` 接收），即样本名。
- 提交成功会返回：`Submitted batch job 1970XX` —— **记下这个作业号**。

## 4. 监控

```bash
squeue -u zhllab_student
# ST 列: PD=排队中(Pending)  R=运行中(Running)  不见了=已结束

# 实时看日志（<jobid> 换成你的作业号）：
tail -f /share/home/zhllab_student/yzx/260727_singleworm/logs/zumis-test.<jobid>.log
# 按 Ctrl+C 退出 tail，不影响作业运行
```

## 5. 判断成功与失败

```bash
# 作业结束后看最后几行：
tail -20 /share/home/zhllab_student/yzx/260727_singleworm/logs/zumis-test.<jobid>.log
# 成功最后一行是: [DONE] ... day11_CF_4

# 检查产物（最关键的是 UMI counts 和统计）：
ls /share/home/zhllab_student/yzx/260727_singleworm/06zumis/day11_CF_4/zUMIs_output/expression/
# 应看到 day11_CF_4.dgecounts.rds 等文件

# 查看作业实际消耗（以后调资源申请的依据）：
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,AllocCPUS
```

常见失败原因：
- 日志里 `[FAIL] 文件不存在` → 检查样本名拼写、fq 文件路径
- `State=OUT_OF_MEMORY` → 调大 `--mem` 和 yaml 里的 `mem_limit`
- `State=TIMEOUT` → 调大 `--time`

失败了修好直接**重新提交同一条 sbatch 命令**即可；已完成的样本会被自动跳过。

## 6. 测试通过后的全量批量

```bash
sbatch /share/home/zhllab_student/yzx/260727_singleworm/scripts/zumis-batch-fat.sbatch
```

批量版的关键区别（可以打开两个脚本对比学习）：

```bash
#SBATCH --partition=compute_fat    # 换到 fat 分区（160 核/1-2TB 大内存节点）
#SBATCH --array=1-785%2            # 作业数组: 785 个任务, 最多 2 个同时跑
#SBATCH --cpus-per-task=128
#SBATCH --mem=500G
```

- `--array=1-785%2`：Slurm 会启动 785 个"子任务"，每个任务里
  `SLURM_ARRAY_TASK_ID` 分别是 1..785；脚本用这个数字从 `samples.txt`
  里取出第 N 行作为本任务要处理的样本。`%2` 限制最多 2 个并发
  （fat01、fat02 各占一个 128 核任务）。
- 批量监控：`squeue -u zhllab_student` 会显示 `zumis-batch[15-785%2]` 这样的数组作业；
  日志是每个子任务一个文件 `zumis-batch.<主jobid>_<子任务号>.log`。

## 7. 常用命令速查

```bash
sbatch <脚本> [参数...]     # 提交
squeue -u $USER             # 看自己的作业
scancel <jobid>             # 取消作业（数组作业 scancel <主jobid> 全取消）
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # 查历史作业消耗
sinfo -N -o '%N %T %c %m'   # 看节点状态（idle=空闲）
```
