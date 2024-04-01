#!/bin/bash
#./runGLIBC.sh platform
echo " ------------ Run time: $(date) ------------" >> nohup.out

if [ $# -ne 1 ]; then
  echo "Usage: ./runGLIBC.sh platform"
  exit
fi

platform=$1
obj_dir=""
if [ ${platform} = 'k' ]; then
  obj_dir="cmake-build-release-kp"
else
  obj_dir="cmake-build-release-intel"
fi
# double test
log_dir="glibc/double/logs"
mkdir -p glibc/double/logs
exec_pref="../${obj_dir}/correctness_test/glibc/double"

commands=(
"${exec_pref}/correct_test_glibc_dbl_Exp2 /home/yixin/wxl/rlibm-32-main/oracles/Exp2Oracle" \
"${exec_pref}/correct_test_glibc_dbl_Exp /home/yixin/wxl/rlibm-32-main/oracles/ExpOracle" \
"${exec_pref}/correct_test_glibc_dbl_Exp10 /home/yixin/wxl/rlibm-32-main/oracles/Exp10Oracle" \
"${exec_pref}/correct_test_glibc_dbl_Log2 /home/yixin/wxl/rlibm-32-main/oracles/Log2Oracle" \
"${exec_pref}/correct_test_glibc_dbl_Log /home/yixin/wxl/rlibm-32-main/oracles/LogOracle" \
"${exec_pref}/correct_test_glibc_dbl_Log10 /home/yixin/wxl/rlibm-32-main/oracles/Log10Oracle" \
)

run_cmd() {
    cmd=$3
    log_dir="glibc/double/logs"
    filename=$(echo $cmd | awk '{print $1}' | cut -d'/' -f6 | cut -d '_' -f 5)
    $cmd > "${log_dir}/${filename}.log" 2>&1
}

export -f run_cmd

num_cores=64

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

#parallel -j ${num_jobs} run_cmd ::: "${num_cores}" ::: "${i}" ::: "${commands[@]}"

#float test
echo "testing float"
log_dir="glibc/float/logs"
mkdir -p ${log_dir}
exec_pref="../${obj_dir}/correctness_test/glibc/float"

commands2=(
"${exec_pref}/correct_test_glibc_flt_Exp2 /home/yixin/wxl/rlibm-32-main/oracles/Exp2Oracle" \
"${exec_pref}/correct_test_glibc_flt_Exp /home/yixin/wxl/rlibm-32-main/oracles/ExpOracle" \
"${exec_pref}/correct_test_glibc_flt_Exp10 /home/yixin/wxl/rlibm-32-main/oracles/Exp10Oracle" \
"${exec_pref}/correct_test_glibc_flt_Log2 /home/yixin/wxl/rlibm-32-main/oracles/Log2Oracle" \
"${exec_pref}/correct_test_glibc_flt_Log /home/yixin/wxl/rlibm-32-main/oracles/LogOracle" \
"${exec_pref}/correct_test_glibc_flt_Log10 /home/yixin/wxl/rlibm-32-main/oracles/Log10Oracle" \
)

run_cmd2() {
    cmd=$3
    log_dir="glibc/float/logs"
    filename=$(echo $cmd | awk '{print $1}' | cut -d'/' -f6 | cut -d '_' -f 5)
    $cmd > "${log_dir}/${filename}.log" 2>&1
}

export -f run_cmd2
parallel -j ${num_jobs} run_cmd2 ::: "${num_cores}" ::: "${i}" ::: "${commands2[@]}"

exit
