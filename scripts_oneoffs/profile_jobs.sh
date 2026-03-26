#!/bin/bash
# Profile completed SLURM jobs using seff and report a summary table.
#
# Usage: bash profile_jobs.sh <job_log>
#   job_log : tab-delimited file with columns: label  jobid
#             (one row per job or array job task)
#
# Output: tab-delimited table to stdout with per-job resource stats.
#         Redirect to a file to save: bash profile_jobs.sh jobs.log > jobs.profile.txt
#
# The job_log file can be built as you submit jobs, e.g.:
#   jid=$(sbatch ... | awk '{print $4}')
#   echo -e "founderSites_chrX\t${jid}_1" >> jobs.log
#
# For array jobs, append _<task> to the job ID to profile individual tasks,
# or omit the task suffix to get aggregate stats for the whole array.

job_log=$1

if [[ -z "$job_log" || ! -f "$job_log" ]]; then
  echo "Usage: bash profile_jobs.sh <job_log>"
  echo "job_log format: two tab-separated columns: label<TAB>jobid"
  exit 1
fi

# Print header
printf "%-35s %-12s %-12s %-10s %-8s %-16s %-10s %-16s %-10s\n" \
  "label" "jobid" "state" "wall_time" "cpu_eff" "cpu_utilized" "mem_eff" "mem_utilized" "mem_requested"

printf "%-35s %-12s %-12s %-10s %-8s %-16s %-10s %-16s %-10s\n" \
  "---" "---" "---" "---" "---" "---" "---" "---" "---"

while IFS=$'\t' read -r label jobid; do
  [[ -z "$jobid" || "$label" =~ ^# ]] && continue

  # Run seff and capture output
  seff_out=$(seff "$jobid" 2>/dev/null)

  if [[ -z "$seff_out" ]]; then
    printf "%-35s %-12s %-12s\n" "$label" "$jobid" "NOT_FOUND"
    continue
  fi

  state=$(echo "$seff_out"       | awk '/^State:/         { print $2 }')
  wall=$(echo "$seff_out"        | awk '/^Job Wall-clock/ { print $4 }')
  cpu_eff=$(echo "$seff_out"     | awk '/^CPU Efficiency/ { print $3 }')
  cpu_used=$(echo "$seff_out"    | awk '/^CPU Utilized:/  { print $3 }')
  mem_eff=$(echo "$seff_out"     | awk '/^Memory Efficiency/ { print $3 }')
  mem_used=$(echo "$seff_out"    | awk '/^Memory Utilized:/  { print $3" "$4 }')
  mem_req=$(echo "$seff_out"     | awk '/^Memory Efficiency/ { print $(NF-1)" "$NF }')

  printf "%-35s %-12s %-12s %-10s %-8s %-16s %-10s %-16s %-10s\n" \
    "$label" "$jobid" "$state" "$wall" "$cpu_eff" "$cpu_used" "$mem_eff" "$mem_used" "$mem_req"

done < "$job_log"
