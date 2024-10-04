
# SSH and Setup Instructions

**Connect to the login node:**
```bash
ssh <username>@login00.broadinstitute.org
use UGER
ish -l h_vmem=16G
```

## First-Time Setup:

1. **Transfer Files or Scripts from Cloud:**

   ```bash
   use Google-Cloud-SDK
   gcloud storage cp gs://example .
   ```

2. **Generate SSH Key and Copy to Cluster:**
   ```bash
   ssh-keygen -t rsa -b 4096
   ssh-copy-id <username>@slurm-bits-bigmem-d002
   ssh <username>@slurm-bits-bigmem-d002
   ```

3. **Install and Configure Conda:**
   ```bash
   mkdir -p ~/miniconda3
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
   bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
   rm -rf ~/miniconda3/miniconda.sh
   ~/miniconda3/bin/conda init bash
   source ~/.bashrc
   conda update -n base conda
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba
   conda create --name scpy
   conda update -n base conda
   conda activate scpy
   conda install -c conda-forge scanpy python-igraph leidenalg
   pip install scanpy
   pip install harmonypy
   ```

4. **Install Jupyter Lab and R:**
   ```bash
   conda install -c conda-forge jupyterlab
   # or
   pip install jupyterlab
   pip install jupyter-resource-usage
   conda install conda-forge::r-base r-essentials
   conda install anaconda::ipykernel
   R
   IRkernel::installspec()
   ```



   

## Finally:

1. **Connect to the compute node and Jupyter:**
   ```bash
   ssh <username>@slurm-bits-bigmem-d002
   srun --nodes=1  --mem=128GB  --time=06:00:00 --cpus-per-task=8 --pty /bin/bash
   HOST_IP=$(hostname -i)
   conda activate scpy
   jupyter lab --ip $HOST_IP --port 8790 --no-browser &
   ```

2. **Go to browser and replace IP with what is printed when you run above command:**
   ```bash
   http://10.192.XX.XX:8790/lab
   ```
