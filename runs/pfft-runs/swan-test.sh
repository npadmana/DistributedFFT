#PBS -l select=32,place=scatter
#PBS -l walltime=04:00:00
#PBS -j oe
#PBS -N swan-pfft-test

runit() {
    NG=$1
    PROCX=$2
    PROCY=$3
    NLOOPS=$4
    NP=$(( ${PROCX}*${PROCY} ))
    EXEC=./pfft/tests/time_c2c_transposed
    echo --------
    echo Running Ngrid=${NG}, processor grid=${PROCX} ${PROCY}, cpus=${NP}, loops=${NLOOPS}
    aprun -n ${NP} -N 32 ${EXEC} -pfft_opt 1 -pfft_tune 1 -pfft_n ${NG} ${NG} ${NG} -pfft_loops ${NLOOPS} -pfft_np ${PROCX} ${PROCY}
    echo ----
}


cd $PBS_O_WORKDIR

runit 1024 2 32 20
runit 1024 4 32 20
runit 1024 8 32 20
runit 1024 32 8 20
runit 1024 16 16 20
runit 1024 16 32 20
runit 1024 32 16 20
runit 1024 32 32 20

