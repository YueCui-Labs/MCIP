# MCIP: Multimodal Connectivity-based Individualized Parcellation
This is the code for paper: Multimodal connectivity-based individualized parcellation and analysis for humans and rhesus monkeys
Now published on IEEE Transcations on Medical Imaging: https://ieeexplore.ieee.org/document/10508267

MCIP simultaneously optimizes a single subject’s within-region homogeneity with the fusion of functional and anatomical connectivity, spatial continuity, and the similarity to a reference atlas.

<img width="500" alt="image" src="https://github.com/YueCui-Labs/MCIP/assets/41955813/89a14cc1-cd61-4e45-921e-f61aeaad0fce">

## How to use it

### MATLAB Dependencies
- Reading GIfTI files in MATLAB: https://www.gllmflndn.com/software/matlab/gifti/
- Reading CIfTI files in MATLAB: https://github.com/Washington-University/cifti-matlab
- GCoptimization - software for energy minimization with graph cuts: https://github.com/nsubtil/gco-v3.0

### Data preparation
Set `rfMRI_dir` as your rfMRI data directory, `dMRI_dir` as your dMRI data directory, and `out_dir` as your output directory in `batch_mcip.sh`.

Revise `ref_atlas_file`, `LUT_file`, and `neighbor_file` in `batch_mcip.sh`, if you want to use another reference atlas. `neighbor_file` can be automatically generated during the program running.

### Using
Run the code with a batch of subjects using Slurm: `sbatch -a 1-${N} batch_mcip.sh $sub_list`, in which `$sub_list` is a .txt file containing subject ID.

If you do not use Slurm, please change `sub_num=...` to be a for loop.

### DATA RELEASE
MCIP-derived individualized parcellations based on Glasser atlas and Brannetome atlas for HCP subjects is available on https://drive.google.com/file/d/1H0mV6Z4icdO9QSda7lPx0TgCO_1c_Dm8/view?usp=drive_link

### :books: Citation
Please cite the [following paper](https://doi.org/10.1109/TMI.2024.3392946) when using MCIP:

Cui Yue, Li Chengyi, Lu Yuheng, et al. Multimodal Connectivity-based Individual Parcellation and Analysis for Humans and Rhesus Monkeys[J]. IEEE Transactions on Medical Imaging, 2024.

