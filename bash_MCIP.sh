#! /bin/bash
#SBATCH --job-name=MCIP                         # Job name
#SBATCH --output log/MCIP_%J.log                # Output log file
#SBATCH -N 1                      # 作业申请 1 个节点
#SBATCH --ntasks-per-node=1                        # Run on a single CPU
#SBATCH --cpus-per-task=6         # 单任务使用的 CPU 核心数为 4
set -ux

WD=$1 #HCP1200 or HCP_retest
shift
sub_list=$1
shift
label_w=$1
shift
smooth_w=$1
shift
lambda=$1
shift
atlas=$1
#atlas=BN_Atlas
#atlas=Glasser
shift
ts=$1

MESH=32k

if [[ -n $ts ]]; then
  method=MCIP_${atlas}_a${label_w}_b${smooth_w}_l${lambda}_ts${ts}
else
  method=MCIP_${atlas}_a${label_w}_b${smooth_w}_l${lambda}
fi

rfMRI_dir=~/HCP/$WD/rfmri
dti_dir=~/HCP/$WD/dmri

surf_file=../necessary/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
LUT_file=../necessary/${atlas}_LUT_wb.txt

GIFTI=/data0/user/cyli/matlab/GIFTI
CIFTI=/data0/user/cyli/matlab/cifti-matlab
GCO=/data0/user/cyli/indiv_parcel/gco-v3.0/matlab

#matlab_cmd=$(echo "$method" | tr -d '.')_${WD}_$(basename ${sub_list%.*})_$(date +"%Y_%m_%d_%H_%M_%S").m
matlab_cmd=$(echo "$method" | tr -d '.')_${WD}_$(basename ${sub_list%.*}).m
echo "addpath('${GIFTI}');addpath('${CIFTI}');addpath('${GCO}');" > log/$matlab_cmd

for sub in $(cat $sub_list); do
  echo "subject: $sub"

  sub_dir=../indiv_atlas/$WD/$sub
  mkdir -p $sub_dir

  function_ts_list_file=$sub_dir/function_ts_list_file.txt

  if [ ! -f $function_ts_list_file ]; then
  for session in rfMRI_REST1_LR rfMRI_REST1_RL rfMRI_REST2_LR rfMRI_REST2_RL; do
    func_ts_file=$rfMRI_dir/$sub/${session}_Atlas_MSMAll_hp2000_clean.dtseries.nii
    if [ -f $func_ts_file ]; then echo $func_ts_file; fi
  done > $function_ts_list_file
  fi

  for hemi in L R; do
    #    generate individual parcellation
    ana_conn_dot_file=$dti_dir/$sub/ptx_NxN_${MESH}_${hemi}_MSMAll/fdt_matrix3.dot
    ana_conn_mat_file=$dti_dir/$sub/ptx_NxN_${MESH}_${hemi}_MSMAll/fdt_matrix3.mat
    if [ -f $ana_conn_mat_file ]; then ana_conn_file=$ana_conn_mat_file;
    elif [ -f $ana_conn_dot_file ];then ana_conn_file=$ana_conn_dot_file;
    else continue;
    fi

    ref_atlas_file=../necessary/fsaverage.${hemi}.${atlas}.32k_fs_LR.label.gii
    neighbor_file=../necessary/neighbor.${hemi}.${atlas}.32k_fs_LR.mat
    topo_file=../necessary/distance_map.${hemi}.${atlas}.32k_fs_LR.func.gii
    searchlight_file=../necessary/fsaverage.${hemi}.${atlas}_search_light_r20mm.32k_fs_LR.func.gii

    metric_file=$sub_dir/${method}_${hemi}.32k_fs_LR.func.gii
#    metric_file=$sub_dir/${method}_${hemi}_lambdas.mat
    if [ ! -f $metric_file ]; then
      if [[ -n $ts ]]; then
        echo "disp('$sub'); try; MCIP_ts_drive('$hemi', '$function_ts_list_file', '$ana_conn_file', '$ref_atlas_file', '$surf_file', '$neighbor_file', '$topo_file', $label_w, $smooth_w, $lambda, $ts, '$sub_dir', '$method'); end;" >> log/$matlab_cmd
      else
        echo "disp('$sub'); try; MCIP_drive('$hemi', '$function_ts_list_file', '$ana_conn_file', '$ref_atlas_file', '$surf_file', '$neighbor_file', '$topo_file', $label_w, $smooth_w, $lambda, '$sub_dir', '$method'); end;" >> log/$matlab_cmd
      fi
#      echo "disp('$sub'); try; MCIPfast_drive('$hemi', '$function_ts_list_file', '$ana_conn_file', '$ref_atlas_file', '$surf_file', '$neighbor_file', '$topo_file', '$searchlight_file', $label_w, $smooth_w, '$sub_dir', '$method'); end;" >> $matlab_cmd

    fi

  done
done

matlab  -nodisplay -nosplash -r "addpath('log');${matlab_cmd%.*}; exit;"

for sub in $(cat $sub_list); do
  sub_dir=../indiv_atlas/$WD/$sub
  for hemi in L R; do
    metric_file=$sub_dir/${method}_${hemi}.32k_fs_LR.func.gii
    label_file=$sub_dir/${method}_${hemi}.32k_fs_LR.label.gii

    # func.gii -> label.gii
    if [ $hemi == L ]; then CORTEX=CORTEX_LEFT; else CORTEX=CORTEX_RIGHT; fi
    if [ ! -f $label_file ]; then
      wb_command -set-structure $metric_file ${CORTEX}
      wb_command -metric-label-import $metric_file $LUT_file $label_file
    fi
  done
done


