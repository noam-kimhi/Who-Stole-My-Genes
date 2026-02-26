#!/bin/bash
#SBATCH --job-name=mmseqCBIO
#SBATCH --output=runmmseq.out
#SBATCH --error=runmmseq.err
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#module load cuda/12.2

export PIP_CACHE_DIR="/cs/labs/dina/dorstein/cache_folder"
export CONDA_PKGS_DIRS="/cs/labs/dina/dorstein/cache_folder"
#conda activate clip
#source "/usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6/etc/profile.d/conda.sh"
#conda_path="/usr/local/spack/opt/spack/linux-debian12-x86_64/gcc-12.2.0/miniconda3-24.3.0-iqeknetqo7ngpr57d6gmu3dg4rzlcgk6"
#conda_total_path="$conda_path/etc/profile.d/conda.sh"
#source "$conda_total_path"
#conda activate clip_env
fasta_folder="/cs/76562/davidgersh2004/Miller/Dor/CBIOHACK"
fasta_file="$fasta_folder/combined.fasta"
out="$fasta_folder/combined_out_cluster_mode_1_cov_mode_1"
cd /cs/76562/davidgersh2004/Miller/Dor/CBIOHACK

/sci/labs/fora/share/softwares/mmseqs/bin/mmseqs easy-cluster "$fasta_file" "$out" "$fasta_folder/mmseq_temp" --min-seq-id 0.5 -c 0.80 --cov-mode 1 -e 1e-10 --cluster-mode 1 --alignment-mode 3 --threads ${SLURM_CPUS_PER_TASK} 
#python3 calcEntropy.py
