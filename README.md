# MCIP: Multimodal Connectivity-based Individualized Parcellation
### Article: Multimodal connectivity-based individualized parcellation and analysis for humans and rhesus monkeys

MCIP simultaneously optimizes a single subject’s within-region homogeneity with the fusion of functional and anatomical connectivity, spatial continuity, and the similarity to a reference atlas.
<img width="394" alt="image" src="https://github.com/YueCui-Labs/MCIP/assets/41955813/89a14cc1-cd61-4e45-921e-f61aeaad0fce">

1. Download and extract necessary libraries to `./`:

> Reading GIfTI files in MATLAB: https://www.gllmflndn.com/software/matlab/gifti/
>
> Reading CIfTI files in MATLAB: https://github.com/Washington-University/cifti-matlab
>
> GCoptimization - software for energy minimization with graph cuts: https://github.com/nsubtil/gco-v3.0

2. Set `rfMRI_dir` as your rfMRI data directory, `dMRI_dir` as your dMRI data directory, and `out_dir` as your output directory in `batch_mcip.sh`.

3. Revise `ref_atlas_file`, `LUT_file`, and `neighbor_file` in `batch_mcip.sh`, if you want to use another reference atlas. `neighbor_file` can be automatically generated during the program running 

4. Run the code with a batch of subjects using Slurm: `sbatch -a 1-${N} batch_mcip.sh $sub_list`, in which `$sub_list` is a .txt file containing subject ID.
If you do not use Slurm, please change `sub_num=...` to be a for loop.
