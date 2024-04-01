#!/bin/bash
#Usage: ./runRLIBM.sh platform
echo -e "Performance testing"
if [ $# -ne 1 ]; then
  echo "Usage: ./runRLIBM.sh platform"
  exit
fi

cd ../
if [ $1 = 'k' ]; then
  obj_dir="cmake-build-release-kp"
  log_dir="kp"
else
  obj_dir="cmake-build-release-intel"
  log_dir="intel"
fi

. ./build.sh $1
obj_dir_pref="../../${obj_dir}/performance_test/rlibm"

logs="logs-${log_dir}"
cd performance_test/rlibm || exit
mkdir -p ${logs}

echo "taskset -c 2 perf_test_log"
taskset -c 2 "${obj_dir_pref}/perf_test_log_estrin_fma" > ./${logs}/log_estrin_fma_no_ig.txt

echo "taskset -c 2 ./perf_test_log2"
taskset -c 2 "${obj_dir_pref}/perf_test_log2_estrin_fma" > ./${logs}/log2_estrin_fma_no_ig.txt

echo "taskset -c 2 ./perf_test_log10"
taskset -c 2 "${obj_dir_pref}/perf_test_log10_estrin_fma" > ./${logs}/log10_estrin_fma_no_ig.txt

echo "taskset -c 2 ./perf_test_exp"
taskset -c 2 "${obj_dir_pref}/perf_test_exp_estrin_fma" > ./${logs}/exp_estrin_fma_no_ig.txt

echo "taskset -c 2 ./perf_test_exp2"
taskset -c 2 "${obj_dir_pref}/perf_test_exp2_estrin_fma" > ./${logs}/exp2_estrin_fma_no_ig.txt

echo "taskset -c 2 ./perf_test_exp10"
taskset -c 2 "${obj_dir_pref}/perf_test_exp10_estrin_fma" > ./${logs}/exp10_estrin_fma_no_ig.txt
