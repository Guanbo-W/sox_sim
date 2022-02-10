#!/bin/sh
#SBATCH --time=02:0:00
#SBATCH --job-name="500_ABAB"
#SBATCH --account=...
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=15G

cd ...

module --force purge
module load nixpkgs/16.09  
module load gcc/8.3.0
module load StdEnv/2020
module load tktable/2.10
module load gcc/9.3.0
module load r/4.1.2


# Export the nodes names. 
# If all processes are allocated on the same node, NODESLIST contains : node1 node1 node1 node1
# Cut the domain name and keep only the node name
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))
Rscript --vanilla 500_ABAB.R $1 $2
