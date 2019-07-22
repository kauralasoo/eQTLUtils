### Build the Docker container
sudo docker build -t kauralasoo/eqtlutils:<latest_git_commit_hash> .

### To get the hash of the latest git commit, run:
# https://blog.container-solutions.com/tagging-docker-images-the-right-way
git log -1 --pretty=%H

### Push to DockerHub
docker push kauralasoo/eqtlutils

### Build a local copy of the Singularity container
singularity build qtl_norm_qc.img docker://kauralasoo/eqtlutils:latest

### Bind the /gpfs/hpc file system to all singularity containers (helps with testing)
### You can add this line to your ~/.bashrc file
export SINGULARITY_BINDPATH="/gpfs/hpc"
