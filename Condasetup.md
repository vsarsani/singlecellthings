
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
   conda install conda-forge::r-base
   conda install anaconda::ipykernel
   R
   IRkernel::installspec()
   ```

5. **Configure Jupyter Lab:**
   ```bash
   export HOSTADDR=$(hostname -i)
   jupyter lab --generate-config
   jupyter lab password
   ```

   - Edit the configuration file: `/home/unix/<username>/.jupyter/jupyter_lab_config.py`
     - First line: add `import os`
     - Line 911: uncomment and replace with:
       ```python
       c.ServerApp.ip = os.environ['HOSTADDR']
       ```
     - Line 960: do the same
       ```python
       c.ServerApp.local_hostnames = [os.environ['HOSTADDR']]
       ```

6. **Add a Jupyter Lab Shortcut to `.bashrc`:**
   - Choose ports like 8788, 8789, 8790, etc.
   ```bash
   echo "alias startjupyter='jupyter lab --ip 10.192.4.44 --port 8788 --no-browser'" >> ~/.bashrc
   source ~/.bashrc
   ```

## Finally:

1. **Connect to the compute node:**
   ```bash
   ssh <username>@slurm-bits-bigmem-d002
   ```

2. **Start Jupyter Lab:**
   ```bash
   startjupyter
   ```
