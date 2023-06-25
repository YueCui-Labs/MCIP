#! /bin/bash
#SBATCH -N 1
#SBATCH -p cpu
#SBATCH -J mcip_parc
#SBATCH --mem=4G
set -ux

WD=HCP1200 #or HCP_retest

sub_list=$1

label_w=0.04
smooth_w=0.03

method=mcip_l${label_w}_w${smooth_w}
MESH=32k
DOWNSIZE=5mm

if [ $WD == HCP1200 ]; then
  rfMRI_dir=/n02dat01/users/cyli/HCP/HCP1200/Data
elif [ $WD == HCP_retest ]; then
  rfMRI_dir=/n02dat01/users/cyli/HCP/HCP_retest/rfMRI
else
  echo $WD not agree!
  exit 0
fi
surf_file=necessary/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii
LUT_file=necessary/BN_Atlas_210_LUT_wb.txt
#LUT_file=necessary/Glasser360_LUT_wb.txt


GIFTI=gifti #https://www.gllmflndn.com/software/matlab/gifti/
CIFTI=cifti-matlab #https://github.com/Washington-University/cifti-matlab
GCO=gco-v3.0/matlab #https://github.com/nsubtil/gco-v3.0

sub_num=${SLURM_ARRAY_TASK_ID}
sub=$(sed -n ${sub_num}p $sub_list)
echo "subject: $sub"

sub_dir=$WD/$sub

function_ts_list_file=$sub_dir/function_ts_list_file.txt
if [ -f $function_ts_list_file ];then rm -rf $function_ts_list_file; fi

for session in rfMRI_REST1_LR rfMRI_REST1_RL rfMRI_REST2_LR rfMRI_REST2_RL; do
  if [ $WD == HCP1200 ]; then
    func_ts_file=$rfMRI_dir/$sub/MNINonLinear/Results/$session/${session}_hp2000_clean_filter_smooth_FWHM4.dtseries.nii
#      func_ts_file=$rfMRI_dir/$sub/MNINonLinear/Results/$session/${session}_Atlas_MSMAll_hp2000_clean_filter_smooth_FWHM4.dtseries.nii
  elif [ $WD == HCP_retest ]; then
    func_ts_file=$rfMRI_dir/$sub/$session/${session}_hp2000_clean_filter_smooth_FWHM4.dtseries.nii
#      func_ts_file=$rfMRI_dir/$sub/$session/${session}_Atlas_MSMAll_hp2000_clean_filter_smooth_FWHM4.dtseries.nii
  fi
  if [ -f $func_ts_file ]; then echo $func_ts_file >> $function_ts_list_file; fi
done

out_dir=$sub_dir/indiv_atlas
mkdir -p $out_dir

for hemi in L R; do
  #    generate individual parcellation
  ana_conn_dot_file=$sub_dir/probtrackx_${MESH}_${hemi}_${DOWNSIZE}/fdt_matrix2.dot
  ana_conn_mat_file=$sub_dir/probtrackx_${MESH}_${hemi}_${DOWNSIZE}/fdt_matrix2.mat
  if [ -f $ana_conn_mat_file ]; then ana_conn_file=$ana_conn_mat_file;
  elif [ -f $ana_conn_dot_file ];then ana_conn_file=$ana_conn_dot_file;
  else exit 0;
  fi

  ref_atlas_file=necessary/fsaverage.${hemi}.BN_Atlas.32k_fs_LR.label.gii
#  ref_atlas_file=necessary/fsaverage.${hemi}.Glasser.32k_fs_LR.label.gii

  neighbor_file=necessary/neighbor.${hemi}.BN_Atlas.32k_fs_LR.mat
#  neighbor_file=necessary/neighbor.${hemi}.Glasser.32k_fs_LR.mat

  metric_file=$out_dir/${method}_${hemi}.func.gii
  matlab  -nodisplay -nosplash -r "addpath('${GIFTI}');addpath('${CIFTI}');addpath('${GCO}');mcip_drive('$hemi', '$function_ts_list_file', '$ana_conn_file', '$ref_atlas_file', '$surf_file', '$neighbor_file', $label_w, $smooth_w, '$out_dir', '$method');exit"

  #  func.gii -> label.gii
  if [ $hemi == L ]; then CORTEX=CORTEX_LEFT; else CORTEX=CORTEX_RIGHT; fi
  label_file=$out_dir/${method}_${hemi}.label.gii
  if [ ! -f $label_file ]; then
    wb_command -set-structure $metric_file ${CORTEX}
    wb_command -metric-label-import $metric_file $LUT_file $label_file
  fi

done
