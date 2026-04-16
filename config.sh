#### CONFIG FILE -- USER SPECIFIC PATHS

# please add needed file paths compatible with your system to your relevent conda source and environments

### conda setup

# for ortholog mapping -- Hal environment (see READ.ME for installation instructions)
HAL_CONDA_SOURCE="/opt/packages/anaconda3-2024.10-1/etc/profile.d/conda.sh"
HAL_ENV="/jet/home/sturecki/.conda/envs/hal"
HALPER_DIR="/jet/home/sturecki/repos/halLiftover-postprocessing"

# for Gene Ontology -- rgreat environment (see READ.ME for installation instructions)
RGREAT_CONDA_SOURCE="/jet/home/sdhanuka/miniconda3/etc/profile.d/conda.sh"
RGREAT_ENV="/ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env"
RGREAT_LIBPATH="/ocean/projects/bio230007p/sdhanuka/conda_envs/rgreat_env/lib/R/library"

# for enhancer vs. promoter and motif analysis -- HOMER environment (see READ.ME for installation instructions)
HOMER_CONDA_SOURCE="/opt/packages/anaconda3-2024.10-1/etc/profile.d/conda.sh"
HOMER_ENV="/ocean/projects/bio230007p/sturecki/project/conda_envs/homer_env"
