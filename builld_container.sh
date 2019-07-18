### Build the Docker container
sudo docker build -t kauralasoo/eqtlutils .

### Push to DockerHub
docker push kauralasoo/eqtlutils

### Build a local copy of the Singularity container
singularity build qtl_norm_qc.img docker://kauralasoo/eqtlutils:latest

### Bind the /gpfs/hpc file system to all singularity containers (helps with testing)
### You can add this line to your ~/.bashrc file
export SINGULARITY_BINDPATH="/gpfs/hpc"
