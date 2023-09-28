
module load python
conda activate /global/cfs/cdirs/act/software/iabril/condaenvs/BBenv
cd /global/cfs/cdirs/act/software/iabril/condaenvs/github_reps/BBPower

lmin=2
lmax=502
dell=5
experiment=S4_SAT
nside=256

base_path=/global/cfs/cdirs/act/data/iabril/newr/cells_models
sacc_tracer=sacc_${experiment}_${lmin}_${lmax}_${dell}.fits

cell_path=${base_path}/${sacc_tracer}

fgs='00s d00'

for fg in ${fgs}
do
decorr=0
cross=0

name_config=${fg}_${decorr}_${cross}_${nside}_${lmin}_${lmax}

outdir=${base_path}/saccs/${name_config}/
if [ ! -d ${outdir} ]; then
    mkdir ${outdir}
fi

python -m bbpower BBCompSep --cells_coadded=${cell_path} --cells_coadded_cov=${cell_path} --cells_noise=${cell_path} --cells_fiducial=${cell_path} --config=${base_path}/${name_config}.yml  --output_dir=${outdir}  
done