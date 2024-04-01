#!/bin/bash
#./runPRLIBMTest.sh platform
echo " ------------ Run time: $(date) ------------" >> nohup.out

if [ $# -ne 1 ]; then
  echo "Usage: ./runPRLIBMTest.sh platform"
  exit
fi

platform=$1
obj_dir=""
if [ ${platform} = 'k' ]; then
  obj_dir="cmake-build-release-kp"
else
  obj_dir="cmake-build-release-intel"
fi

exec_pref="../${obj_dir}/correctness_test"

commands=(
"${exec_pref}/correct_test_Exp2 /home/yixin/wxl/rlibm-all/oracles/Exp2Oracle" \
"${exec_pref}/correct_test_Exp /home/yixin/wxl/rlibm-all/oracles/ExpOracle" \
"${exec_pref}/correct_test_Exp10 /home/yixin/wxl/rlibm-all/oracles/Exp10Oracle" \
"${exec_pref}/correct_test_Log2 /home/yixin/wxl/rlibm-all/oracles/Log2Oracle" \
"${exec_pref}/correct_test_Log /home/yixin/wxl/rlibm-all/oracles/LogOracle" \
"${exec_pref}/correct_test_Log10 /home/yixin/wxl/rlibm-all/oracles/Log10Oracle" \
)

run_cmd() {
    cmd=$3
    filename=$(echo $cmd | awk '{print $1}' | cut -d'/' -f4 | cut -d '_' -f 3)
    $cmd > "./logs/${filename}.log" 2>&1
}

export -f run_cmd

# Acquire environmental var.: number of cores
num_cores=-1
if [ -z "${OMP_NUM_THREADS}" ]; then
    echo "Error: OMP_NUM_THREADS is not set！"
    exit
else
    num_cores=${OMP_NUM_THREADS}
    echo "num_cores: $num_cores"
fi

# Use GNU parallel to run commands in parallel
num_jobs=-1
case $num_cores in
"64") num_jobs=2;;
"32") num_jobs=4;;
"16" | "8" | "4") num_jobs=6;;
*)
  echo "Error: Num of cores must be one of numbers in (64, 32, 16, 8, 4)"
  exit;;
esac

echo "num_jobs: $num_jobs"

parallel -j ${num_jobs} run_cmd ::: "${num_cores}" ::: "${i}" ::: "${commands[@]}"

exit
