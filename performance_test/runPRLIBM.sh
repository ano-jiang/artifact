#!/bin/bash
#Usage: ./runPRLIBM.sh platform
echo -e "Performance testing"
if [ $# -ne 1 ]; then
  echo "Usage: ./runPRLIBM.sh platform"
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
obj_dir_pref="../../${obj_dir}/performance_test/prlibm"

logs="logs-${log_dir}"
cd performance_test/prlibm || exit
mkdir -p ${logs}


echo "taskset -c 2 perf_test_log"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_log" > ./${logs}/log_no_ig.txt

echo "taskset -c 2 ./perf_test_log2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_log2" > ./${logs}/log2_no_ig.txt

echo "taskset -c 2 ./perf_test_log10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_log10" > ./${logs}/log10_no_ig.txt

echo "taskset -c 2 ./perf_test_exp"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_exp" > ./${logs}/exp_no_ig.txt

echo "taskset -c 2 ./perf_test_exp2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_exp2" > ./${logs}/exp2_no_ig.txt

echo "taskset -c 2 ./perf_test_exp10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_exp10" > ./${logs}/exp10_no_ig.txt