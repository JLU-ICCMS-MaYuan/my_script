#!/bin/bash
# iter scheduling remote
#
conda activate base
resource_path=$1
source $resource_path/server.sh
source $resource_path/fun/*

set_server server1
echo_server

fp_calculator="vasp" # support qe, vasp       [DEFAULT: vasp] 

# # # # # # # # # # # # # # # # # # # remote funcions # # # # # # # # # # # # # # # # # # # ##
remove_trailing_slash() {
  local input_string="$1"
  if [[ $input_string == */ ]]; then
    echo "${input_string%/}"
  else
    echo "$input_string"
  fi
}

extract_job_name() {
    local sbatch_line="$1"

    if [[ $sbatch_line =~ --job-name=([a-zA-Z0-9_-]+) ]]; then
        job_name="${BASH_REMATCH[1]}"
        echo "$job_name"
    else
        echo "job-name not found in the string."
        exit 1
    fi
}

submit_remote() {
  # usage: submit_remote "local_task_dir" "submit_task_command"
  # example submit_remote "./fp0000" "sbatch run.sh"
  #
  # res=$(submit_remote "$task" "sbatch run.sh")
  # res="${res##*$'\n'}"
  # res: remote map dir
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

  local local_task_dir="$1"
  local sub_cmd="$2"

  local_task_dir=$(remove_trailing_slash $local_task_dir)
  remote_work_place=$(remove_trailing_slash $remote_work_place)

  #remote_dir=$(od -An -N2 -i /dev/urandom | shasum)  # Mac OS
  remote_dir=$(od -An -N2 -i /dev/urandom | sha1sum)  # Mac OS X
  remote_dir=${remote_dir:0:10}$(date "+%Y%m%d%H%M%S%N")
  remote_dir="tmp."$remote_dir
  remote_dir=$remote_work_place/$remote_dir
  echo "$remote_dir"

  echo "Upload '$local_task_dir' to '$remote_dir'"
  rsync -arPq -e "ssh -i $priv_rsa_path -p $remote_port" $local_task_dir/ $remote_ip:$remote_dir || true
  rsync -arPq -e "ssh -i $priv_rsa_path -p $remote_port" $local_task_dir/ $remote_ip:$remote_dir || true
  ssh -i $priv_rsa_path -p $remote_port $remote_ip "
    cd $remote_dir
    $sub_cmd
  "
  echo "$remote_dir"
}

submit_one(){
  local l_task="$1"
  local submit_cmd="$2"

  qos_limit=100
  while true; do
    if [ "$task_sys" == "slurm" ];then
      queue_count=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "squeue | grep -v USER | wc -l") # slurm
    elif [ "$task_sys" == "pbs" ];then
      sleep 3
      queue_count=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "qstat | grep $task_name | wc -l") # pbs 
    else
      echo "acnn activate learning do not support task_sys $1"
      exit 1
    fi

    if [ "$queue_count" -lt $qos_limit ]; then
      remote_task=$(submit_remote "$l_task" "$submit_cmd")
      remote_task="${remote_task##*$'\n'}"
      echo "$remote_task"
      break
    else
      sleep 60
    fi
  done
}

wait_remote() {
  local task_name="$1"
  if [ "$task_sys" == "slurm" ];then
    jobs_left=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "squeue -n $task_name|grep -v JOBID|wc -l") 
  elif [ "$task_sys" == "pbs" ];then
    sleep 10
    jobs_left=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "qstat | grep $task_name | wc -l")
  else
    echo "acnn activate learning do not support task_sys $1"
    exit 1
  fi

  while [ "$jobs_left" -gt 0 ]; do
    if [ "$jobs_left" -gt 8 ]; then
      sleep 50
    fi
    sleep 10

    if [ "$task_sys" == "slurm" ];then
      jobs_left=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "squeue -n $task_name|grep -v JOBID|wc -l")
    elif [ "$task_sys" == "pbs" ];then
      jobs_left=$(ssh -i $priv_rsa_path -p $remote_port $remote_ip "qstat|grep $task_name|wc -l")
    else
      echo "acnn activate learning do not support task_sys $1"
      exit 1
    fi
    
    echo "waiting remote jobs: $jobs_left"
  done
}

recycling_one(){
  local s_task="$1"
  local l_task="$2"
  echo "recycling remote task ($s_task) to ($l_task)"
  rsync --exclude "WAVECAR" -arPq -e "ssh -i $priv_rsa_path -p $remote_port" $remote_ip:"$s_task"/ "$l_task" || true
  rsync --exclude "WAVECAR" -arPq -e "ssh -i $priv_rsa_path -p $remote_port" $remote_ip:"$s_task"/ "$l_task" || true
  ssh -i $priv_rsa_path -p $remote_port $remote_ip "rm -r $s_task" || true
}

recycling_remote() {
  for ((i = 0; i < ${#remote_task_dirs[@]}; i++)); do
    local s_task=${remote_task_dirs[i]}
    local l_task=${local_task_dirs[i]}
    recycling_one "$s_task" "$l_task"
  done
}
# # # # # # # # # # # # # # # # # # # remote funcions end # # # # # # # # # # # # # # # # # # ##

check_python_version() {
    python_path=$(which python)

    if [ -z "$python_path" ]; then
        echo "Python3 not found."
        return 1
    else
        python_version=$("$python_path" --version 2>&1 | awk '{print $2}')

        major_version=$(echo "$python_version" | cut -d. -f1)
        minor_version=$(echo "$python_version" | cut -d. -f2)

        if [ "$major_version" -lt 3 ]; then
            echo "Python < 3. We need python3."
            return 1
        else
            echo "Python version: $python_version"
            return 0
        fi
    fi
}


wait_task(){
    echo "Waiting...  task_sys: $1  task_name $2"
    if [ "$task_sys" == "donau" ];then
      jobs_cmd="djob -n $2|grep -v JOB_ID|grep $2|wc -l";
    elif [ "$task_sys" == "slurm" ];then
      jobs_cmd="squeue -n $2 | grep -v JOBID | wc -l"
    elif [ "$task_sys" == "pbs" ];then
      jobs_cmd="qstat  | grep $2 | grep -v JOBID | wc -l"
    else
      echo "acnn activate learning do not support task_sys $1"
      exit 1
    fi

    jobs_exec=$(eval "$jobs_cmd")
    wait_time=60
    while( (( $jobs_exec > 0 )) ); do
        echo Jobs left: $jobs_exec
        sleep $wait_time
        if [ "$jobs_exec" -lt 5 ]; then
            wait_time=60
            echo "xxxx"|dconfig
        fi
        jobs_exec=$(eval "$jobs_cmd")
    done
}


terminate_handler() {
    echo "Terminating the task..."
    echo "Note that there may be residual tasks, please kill them manually."
    exit 1
}

handle_error() {
  local exit_code=$?
  local command="$BASH_COMMAND"
  echo "Error occurred in command: '$command', exit code: $exit_code"
  exit $exit_code
}
# End of defined function
# main - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# main - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

set -e
echo "iterate date:  $(date)"
check_python_version
trap 'handle_error' ERR
trap terminate_handler INT
trap terminate_handler SIGTERM SIGINT

print_line="- - - - - - -"

echo "iter_scheduling tast pid: $$"

#- - - - - - - - - - -  - - - - - - - - -  step 0: Prepare
echo
echo "$print_line step 0: Prepare"

last_iter=$(find ./ -maxdepth 1 -type d -regex ".*/[0-9]+" -exec basename {} \; | sort -n | tail -n 1)
if [ -z "$last_iter" ]; then
    iter=0
else
    iter=$((10#$last_iter + 1))
    #iter=$((last_iter))
fi

iter_dir_name=$(printf "%02d" "$iter")
mkdir -p $iter_dir_name

echo "$print_line step 0: over!"

#- - - - - - - - - - -  - - - - - - - - -  step 1: Check duplicates & Make dataset
echo
echo "$print_line step 1: Check duplicates & Make dataset"
mkdir -p "$iter_dir_name/train_dt"
if [ -n "$(ls -A "$iter_dir_name/train_dt")" ]; then
    rm "$iter_dir_name/train_dt"/*
fi

file_count=0
for file in $resource_path/init_dt/*.xsf; do
    if [ -f "$file" ]; then
        cp -n "$file" "$iter_dir_name/train_dt"
        let file_count++ || true
    fi
done
for dir in $(find "$resource_path" -type d -name "*add"); do
    if [ -d "$dir" ]; then
        for file in "$dir"/*.xsf; do
            if [ -f "$file" ]; then
                cp -n "$file" "$iter_dir_name/train_dt"
                let file_count++ || true
            fi
        done
    fi
done

# check same filename
if [ "$file_count" -ne "$(ls -1 "$iter_dir_name/train_dt" | wc -l)" ]; then
    echo "Warning: Copied ($file_count) files mismatch real files in $iter_dir_name/train_dt. Conf. loss."
fi

# check duplicates cfgs.
declare -A file_hashes
find "$iter_dir_name/train_dt" -type f -name "*.xsf" | while read -r file; do
    hash=$(md5sum "$file" | awk '{print $1}')

    if [[ -n ${file_hashes[$hash]} ]]; then
        echo "Warning: same conf.  ${file_hashes[$hash]} and $file"
    else
        file_hashes[$hash]=$file
    fi
done

echo "Train set Total:  $(ls -1 "$iter_dir_name/train_dt" | wc -l)"
echo "$print_line step 1: over!"


#- - - - - - - - - - -  - - - - - - - - -  step 2: Training
echo
echo "$print_line step 2: Training"
cp $resource_path/sub-train.sh $iter_dir_name
cp $resource_path/inv2 $iter_dir_name

cd $iter_dir_name
sed -i "s/traindatapath.*/traindatapath   =   .\/train_dt/" inv2

if [ "$iter" -eq 0 ]; then
  echo "Training from scratch"
elif [ "$iter" -eq 1 ]; then
  last_iter_name=$(printf "%02d" "$((iter - 1))")
  sed -i 's/restart[[:space:]]*= .*/restart             =   true        # restart button/' inv2
  sed -i "s/restart_model[[:space:]]*= .*/restart_model       =   ..\/$last_iter_name\/model\/model-last        # restart model/" inv2
else
  last_iter_name=$(printf "%02d" "$((iter - 1))")
  sed -i 's/restart[[:space:]]*= .*/restart             =   true        # restart button/' inv2
  sed -i "s/restart_model[[:space:]]*= .*/restart_model       =   ..\/$last_iter_name\/model-restart\/model-last        # restart model/" inv2
fi

echo "train step: $(grep ^nbatch inv2|awk '{print $3}')"
echo "save  step: $(grep ^savestep inv2|awk '{print $3}')"

train_step=$(grep ^nbatch inv2 | awk '{print $3}')

echo "Training ... $train_step"
token_get
dsub -s sub-train.sh
while true; do
    sleep 60
    if grep -q "train-finished" 2-train-log; then
        echo "train-finished in 2-train-log"
        break
    #else
    #    echo "Not found model-$train_step in 2-train-log"
    fi
done

cd ..
echo "Training over."
echo "$print_line step2 over!"

#- - - - - - - - - - -  - - - - - - - - -  step 3: MD explore
echo
cp "$resource_path"/sub-select.sh "$iter_dir_name"
echo "$print_line step 3: MD explore"
if [ $iter -eq 0 ]; then
    cd "$iter_dir_name"
    sed -i "s/traindatapath.*/traindatapath   =   .\/train_dt/" inv2
    token_get
    dsub -s sub-select.sh
    while true; do
        sleep 60
        if grep -q "select finished" 4-select-log; then
            echo "select finished in 4-select-log"
            break
        #else
        #    echo "Not found select finished in 4-select-log"
        fi
    done
    cd ..
else
    last_iter_name=$(printf "%02d" "$((iter - 1))")
    cp "$last_iter_name"/Ap_inv* "$iter_dir_name"
fi
#
#
md_task_root="3-md"
mkdir -p "$iter_dir_name"/"$md_task_root"
cp "$resource_path"/md/* "$iter_dir_name"/"$md_task_root"

cd "$iter_dir_name"

last_model_number=0
last_model=""
model_dir="model-restart/model-*"
if [ $iter -eq 0 ];then
  model_dir="model/model-*"
fi

for file in $model_dir; do
    if [ -f "$file" ]; then
        model=$(basename "$file")
        number=${model#"model-"}
        if [ "$number" -gt "$last_model_number" ]; then
            last_model_number=$number
            last_model="$model"
        fi
    fi
done
echo "last model $last_model"

cp model*/"$last_model" $md_task_root
cp Ap_inv* "$md_task_root"
cd $md_task_root
sed -i "s|\(pair_style[[:space:]]*acnn[[:space:]]*\)[^[:space:]]*|\1$last_model|" in.lmp
sed -i -E '/^velocity/ s/[0-9]+\s*$/'"$((RANDOM % 1000))1"'/g' in.lmp

#sed -i "s/model-name/$last_model/g" ares.in
echo "Simulation ..."
# run local >>>>>>>>>>>>>
token_get
dsub -s sub-md.sh || true

while true; do
    sleep 60
    if grep -q "OOS" 3-md-log; then
        echo "OOS file count limit in 3-md-log"
        break
    elif grep -q "Out of sample" 3-md-log; then
        echo "Out of sample in 3-md-log"
        break
    #else
    #    echo "Not found OOS in 3-md-log"
    fi
done
# run local end <<<<<<<<<

cd ..
echo "Simulation over."

cd ..
echo "$print_line step3 over."

#- - - - - - - - - - -  - - - - - - - - -  step 4: Select
echo
echo "$print_line step 4: Select"
cp $resource_path/sub-select.sh $iter_dir_name

cd $iter_dir_name
rm Ap_inv*
mkdir -p pre-select
cp train_dt/* pre-select
cp $md_task_root/oos/* pre-select
sed -i "s/traindatapath.*/traindatapath   =   .\/pre-select/" inv2

echo "Selecting ..."
token_get
dsub -s sub-select.sh
while true; do
    sleep 60
    if grep -q "select finished" 4-select-log; then
        echo "select finished in 4-select-log"
        break
    #else
    #    echo "Not found select finished in 4-select-log"
    fi
done
#tast_name=$(grep "#DSUB -n" sub-select.sh|awk '{print $3}')
#bash task_wait.sh $tast_name donau
# /home/lijx/playground/torchdemo-dev/build/acnn -select ./inv2 > 4-select-log 2>&1
#sh sub-select.sh
echo "Selecting over."

input_string=$(cat ./4-select-log|grep Act)
result=$(echo "$input_string" | sed 's/Act: //')
not_in_a=()
not_in_cnt=0
for file in $result; do
    if [ ! -f "./train_dt/$file" ]; then
        not_in_a+=("$file")
        let not_in_cnt++ || true
    fi
done

echo "Select $not_in_cnt files."
mkdir -p select_label
for file in "${not_in_a[@]}"; do
    cp $md_task_root/oos/$file select_label
done

cd ..
echo "$print_line step4: over."

#- - - - - - - - - - -  - - - - - - - - -  step: 5 FP labeling
echo
echo "$print_line step: 5 FP labeling"
mkdir -p $iter_dir_name/labeling
python $resource_path/convert_xsf_vasp.py $iter_dir_name/select_label $iter_dir_name/labeling

cd $iter_dir_name/labeling
fp_tasks=0
token_get
#--- local run >>> >>>
for dir in $(ls -d */);do
    cp ../../$resource_path/fp/* $dir
    cd $dir
    dsub -s sub-fp.sh
    cd ..
    let fp_tasks++ || true
done
echo "Submit all $fp_tasks fp task."

max_retries=1000
retry_count=0

while true; do
    all_found=true
    for dir in */; do
        if [ -f "$dir/OUTCAR" ]; then
            if ! grep -q "Total CPU time" "$dir/OUTCAR"; then
                all_found=false
                break
            fi
        else
            all_found=false
            break
        fi
    done

    if $all_found; then
        echo "All OUTCAR files have 'Total CPU time'."
        break
    else
        if [ $retry_count -ge $max_retries ]; then
            echo "Reached the maximum retry limit. Deleting folders without 'Total CPU time'..."
            for dir in */; do
                if [ -f "$dir/OUTCAR" ]; then
                    if ! grep -q "Total CPU time" "$dir/OUTCAR"; then
                        echo "Deleting folder: $dir"
                        rm -rf "$dir"
                    fi
                else
                    echo "Deleting folder: $dir"
                    rm -rf "$dir"
                fi
            done
            break
        fi
        retry_count=$((retry_count + 1))
        echo "Not all OUTCAR files have 'Total CPU time'. Sleeping... ($retry_count/$max_retries)"
        sleep 60
    fi
done

cd ../..

echo "$print_line step5: over."

#- - - - - - - - - - -  - - - - - - - - -  step: 6 Collecting
echo
echo "$print_line step: 6 Collecting"
if [ "$fp_calculator" == "vasp" ];then
  python $resource_path/recycling.py $iter_dir_name/labeling $resource_path/"${iter_dir_name}_add" vasp OUTCAR # vasp
elif [ "$fp_calculator" == "qe" ];then
  python $resource_path/recycling.py $iter_dir_name/labeling $resource_path/"${iter_dir_name}_add" qe   # qe
else
  echo "fp_calculator do not suppport $1"
  exit 1
fi 
#python $resource_path/recycling.py $iter_dir_name/labeling $resource_path/"${iter_dir_name}_add"
echo "Collecting $(ls $resource_path/"${iter_dir_name}_add"|wc -l) added file."

echo "$print_line step6: over."

#- - - - - - - - - - -  - - - - - - - - -  well done - - - - - - - -  - -- - - - - - - - - - -- #
echo "AL-Iter $iter_dir_name over normally."
