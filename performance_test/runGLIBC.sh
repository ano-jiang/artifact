#!/bin/bash
#Usage: . ./runGLIBC.sh
echo -e "Performance testing"


obj_dir="cmake-build-release-intel"
log_dir="intel"

cd ../
echo ". ./build.sh $1"
. ./build.sh $1

logs="logs-${log_dir}"
cd performance_test/glibc/double || exit
mkdir -p ${logs}

obj_dir_pref="../../../${obj_dir}/performance_test/glibc/double"
echo "taskset -c 2 perf_test_log"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_log" > ./${logs}/log_no_ig.txt

echo "taskset -c 2 ./perf_test_log2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_log2" > ./${logs}/log2_no_ig.txt

echo "taskset -c 2 ./perf_test_log10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_log10" > ./${logs}/log10_no_ig.txt

echo "taskset -c 2 ./perf_test_exp"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_exp" > ./${logs}/exp_no_ig.txt

echo "taskset -c 2 ./perf_test_exp2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_exp2" > ./${logs}/exp2_no_ig.txt

echo "taskset -c 2 ./perf_test_exp10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_dle_glibc_exp10" > ./${logs}/exp10_no_ig.txt

cd ../float || exit
mkdir -p ${logs}
obj_dir_pref="../../../${obj_dir}/performance_test/glibc/float"
echo "Performing float test"

echo "taskset -c 2 perf_test_log"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_log" > ./${logs}/log_no_ig.txt

echo "taskset -c 2 ./perf_test_log2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_log2" > ./${logs}/log2_no_ig.txt

echo "taskset -c 2 ./perf_test_log10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_log10" > ./${logs}/log10_no_ig.txt

echo "taskset -c 2 ./perf_test_exp"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_exp" > ./${logs}/exp_no_ig.txt

echo "taskset -c 2 ./perf_test_exp2"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_exp2" > ./${logs}/exp2_no_ig.txt

echo "taskset -c 2 ./perf_test_exp10"
taskset -c 2 "${obj_dir_pref}/perf_test_no_ig_glibc_exp10" > ./${logs}/exp10_no_ig.txt