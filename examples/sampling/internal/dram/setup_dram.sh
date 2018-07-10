#! /bin/bash

BASE_DIR="runs"
DATA_DIR="./data"
DATA_FILE_PREFIX=

EXEC_FILE="../../../../build/dram_theta_internal"

PAR_FILE="./dram.par"
PRIOR_FILE="./priors.par"


Ilist=(1)



if [ ! -d "$BASE_DIR" ]; then
	mkdir $BASE_DIR
fi


for I in "${Ilist[@]}"
do

	II="$(printf "%03d" ${I})"

	RUN_DIR="${BASE_DIR}/run_${II}"

	if [ -d "$RUN_DIR" ]; then
		echo "The directory '${RUN_DIR}' exists. Delete the directory and rerun the script."
		#rm -rf $RUN_DIR
	fi

	mkdir $RUN_DIR


	DATA_FILE="${DATA_DIR}/${DATA_PREFIX}${II}.dat"



	cp ${DATA_FILE} "$RUN_DIR/data.txt"

	cp ${EXEC_FILE}    ${RUN_DIR}
	cp ${PAR_FILE}     ${RUN_DIR}/dram.par
	cp ${PRIOR_FILE}  ${RUN_DIR}/priors.par

done



