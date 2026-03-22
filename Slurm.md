# SLURM Cluster Configuration

## Account Information
- **Standard account:** `tdlong_lab`
- **GPU account:** `tdlong_lab_gpu`

## Module Loading
```bash
module load python/3.10.2
module load R/4.2.2
```

## CPU Partitions

Got it! Here's a more concise version focused on the most impactful strategy:

---

## CPU Partitions

| Partition | Max memory per core | Default / Max runtime |
|-----------|--------------------|-----------------------|
| standard  | 6 GB               | 2 day / 14 day        |
| highmem   | 10 GB              | 2 day / 14 day        |
| hugemem   | 18 GB              | 2 day / 14 day        |

### Key Rules
1. **Memory:** To get maximum memory per core for a partition, you must explicitly request it with `--mem-per-cpu=XG`
   - Example: `--mem-per-cpu=10G` for highmem partition max
   
2. **Scaling memory:** Once at partition max memory per core, you can only get more total memory by:
   - Adding more cores: `--cpus-per-task=N` (each core gets the per-core memory)
   - Moving to a higher partition (standard → highmem → hugemem)
   
3. **Runtime:** Jobs are killed when they reach the default time (2 days) unless you specify longer with `--time=`
   - Time is rarely an issue for our projects
   - Max is 14 days if needed

### Before Requesting More Memory: Subset Your Data

Memory issues often come from reading entire files when only a subset is needed (e.g., specific replicates).

**Filter while reading:**
```r
# R: Only load the replicates you need
data <- read_tsv("file.tsv") %>% 
  filter(replicate %in% c(1, 5, 10))

# Or select specific columns
data <- read_tsv("file.tsv", col_select = c(replicate, position, score))
```

```python
# Python: Read specific columns
df = pd.read_csv("file.csv", usecols=['replicate', 'position', 'score'])

# Then filter
df = df[df['replicate'].isin([1, 5, 10])]
```

### Resource Request Examples
```bash
# Standard partition: 2 cores × 6 GB = 12 GB total
--mem-per-cpu=6G --cpus-per-task=2

# Highmem partition: 4 cores × 10 GB = 40 GB total  
--mem-per-cpu=10G --cpus-per-task=4

# Hugemem partition: 3 cores × 18 GB = 54 GB total
--mem-per-cpu=18G --cpus-per-task=3
```

---

Better - more focused on the practical reality of your workflow?
### CPU Job Template
```bash
#!/bin/bash
#SBATCH --job-name=cpu_analysis
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=4:00:00

module load python/3.10.2
module load R/4.2.2

# Your commands here
Rscript scripts/analysis.R --input data.tsv --output results.rds
```

### CPU Array Job Template
```bash
#!/bin/bash
#SBATCH --job-name=array_job
#SBATCH -A tdlong_lab
#SBATCH -p standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=4:00:00
#SBATCH --array=1-20

module load python/3.10.2
module load R/4.2.2

# Use SLURM_ARRAY_TASK_ID in your script
python scripts/process.py --task_id ${SLURM_ARRAY_TASK_ID}
```

## GPU Partitions

### GPU Options
- **Paid GPU (recommended for real runs):** 
  - Use `-p gpu` and `-A tdlong_lab_gpu`
  - Not interruptible
  - Cost: ~$0.70/hour per GPU
  - Request: `--gres=gpu:V100:1` (one V100 GPU)
  - Set time limit (e.g., `-t 12:00:00`) to control costs
  - Even 48-hour runs cost only tens of dollars

- **Free GPU:**
  - Use `-p free-gpu`
  - Good for short runs (<5 min)
  - **Can be interrupted** when someone requests paid GPUs
  - Fine for quick tests, NOT for long production runs

- **gpu-debug:**
  - Short limit (~15 min)
  - Same rate as paid GPU
  - Often available instantly
  - Good for testing (e.g., few epochs) before submitting long job
  - No need to set `-t`; system sets the limit

### GPU Resource Notes
- GPU jobs get 1 CPU (3 GB) by default, which may be insufficient for preprocessing
- Request `--cpus-per-task=4-8` and `--mem-per-cpu=6G` if doing CPU preprocessing (R formatting) in the same job
- PyTorch training runs on GPU; preprocessing runs on CPU
- Consider separating preprocessing (CPU/standard partition) from training (GPU partition) to avoid paying for GPU time during CPU work

### PyTorch and CUDA
- Use `module load python/3.10.2`
- Pip-installed PyTorch usually includes CUDA and will use GPU automatically on GPU nodes
- If `torch.cuda.is_available()` returns False on GPU node, load a cluster PyTorch module with CUDA support

### GPU Job Template
```bash
#!/bin/bash
#SBATCH --job-name=train_gpu
#SBATCH -A tdlong_lab_gpu
#SBATCH -p gpu
#SBATCH --gres=gpu:V100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00

module load python/3.10.2
module load R/4.2.2

# Your commands here
python scripts/train_model.py --input data.pkl --epochs 100 --output model.pt
```

## Parameter Passing

### SLURM Parameters vs Script Parameters
Pass SLURM parameters **before** the script name, then script parameters **after**:

```bash
sbatch --array=1-20 slurm_jobs/train_data.sh \
    CHROM=chr3R \
    DIR=fast_train_data
```

Or within the script, call with named parameters:
```bash
python scripts/analysis.py \
    --input data/input.csv \
    --output results/output.pkl \
    --param1 value1
```

## Quick Reference

**Common submission patterns:**
```bash
# Standard CPU job
sbatch -A tdlong_lab -p standard --cpus-per-task=1 --mem-per-cpu=5G slurm_jobs/job.sh

# Higher memory CPU job (highmem)
sbatch -A tdlong_lab -p highmem --cpus-per-task=4 --mem-per-cpu=10G slurm_jobs/job.sh

# Production GPU job
sbatch -A tdlong_lab_gpu -p gpu --gres=gpu:V100:1 --cpus-per-task=8 --mem-per-cpu=6G -t 12:00:00 slurm_jobs/train.sh

# Quick GPU test
sbatch -A tdlong_lab_gpu -p gpu-debug --gres=gpu:V100:1 --cpus-per-task=4 slurm_jobs/test.sh
```
