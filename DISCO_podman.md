## Copy your files from gcloud to /broad/macosko/data/


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
   mkdir -p /broad/macosko/data/XXX
   gcloud storage cp gs://example /broad/macosko/data/XXX
   ```



# Instructions for Setting Up RStudio and JupyterLab Environment on DISCO

## 1. Log into DISCO

Use the following command to log into the `login00` node of DISCO:

```bash
ssh login00.broadinstitute.org
```

## 2. Check for Available Resources

Check the available resources for a specific node by using this command:

```bash
scontrol show node "slurm-bits-bigmem-d002" | grep "AllocTRES"| tr '=' '	'|tr ',' '	'| awk  '{cores=64-$3; mem=(3584-$5)/1024; print "Remaining cores: " cores "\nRemaining memory (TB):: " mem}'
```

## 3. Request Resources

To request resources, use the following `srun` command:

```bash
srun -C container --nodes=1 --mem=128GB --cpus-per-task=4 --pty --partition=hpcx_macosko --time=02:00:00 /bin/bash
```

## 4. Navigate to Your Directory

Once the resources are allocated, navigate to your working directory in `/broad/macosko`. For example:

```bash
cd /broad/macosko/data/NPH/
```

## 5. Create and Edit Dockerfile

Create a Dockerfile using `nano`:

```bash
nano Dockerfile
```

Copy the following Dockerfile content into the file:

```dockerfile
# Use rocker/tidyverse as base image
FROM docker.io/rocker/tidyverse:latest

# Set environment variables for username and password
ARG USERNAME=your_username
ARG PASSWORD=your_password

# Switch to root to install system packages
USER root

# Install necessary packages for downloading and installing Miniconda
RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

# Create the user and set the password
RUN useradd -m $USERNAME && echo "$USERNAME:$PASSWORD" | chpasswd

# Install Miniconda
RUN mkdir -p /home/$USERNAME/miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /home/$USERNAME/miniconda3/miniconda.sh && \
    bash /home/$USERNAME/miniconda3/miniconda.sh -b -u -p /home/$USERNAME/miniconda3 && \
    rm -rf /home/$USERNAME/miniconda3/miniconda.sh && \
    /home/$USERNAME/miniconda3/bin/conda init bash && \
    echo ". /home/$USERNAME/miniconda3/etc/profile.d/conda.sh" >> /home/$USERNAME/.bashrc

# Update conda and install dependencies
RUN /home/$USERNAME/miniconda3/bin/conda update -n base conda -y && \
    /home/$USERNAME/miniconda3/bin/conda install -n base conda-libmamba-solver -y && \
    /home/$USERNAME/miniconda3/bin/conda config --set solver libmamba && \
    /home/$USERNAME/miniconda3/bin/conda create --name scpy -y && \
    /home/$USERNAME/miniconda3/bin/conda update -n base conda -y && \
    /home/$USERNAME/miniconda3/bin/conda install -n scpy -c conda-forge scanpy python-igraph leidenalg -y && \
    /home/$USERNAME/miniconda3/bin/pip install scanpy harmonypy

# Install JupyterLab and R-related packages
RUN /home/$USERNAME/miniconda3/bin/conda install -n scpy -c conda-forge jupyterlab

# Install R packages and dependencies
RUN R -e "install.packages('Seurat')" && \
    R -e "install.packages(c('BPCells', 'presto', 'glmGamPoi'))" && \
    R -e "if (!require('remotes')) install.packages('remotes')" && \
    R -e "remotes::install_github('satijalab/seurat-data', quiet = TRUE)" && \
    R -e "remotes::install_github('satijalab/azimuth', quiet = TRUE)" && \
    R -e "remotes::install_github('satijalab/seurat-wrappers', quiet = TRUE)" && \
    R -e "if (!require('BiocManager', quietly = TRUE)) install.packages('BiocManager')" && \
    R -e "BiocManager::install(c('limma', 'variancePartition', 'DESeq2', 'Nebulosa', 'MAST'))" && \
    R -e "library('devtools'); install_github('lme4/lme4', dependencies = TRUE)" && \
    R -e "install.packages('devtools')"

# Set the working directory to the home directory
WORKDIR /home/$USERNAME

# Switch back to the created user
USER $USERNAME

# Expose the port for RStudio
EXPOSE 8787

# Set CMD to start RStudio when the container is run
CMD ["/init"]

```
Save it 
6. Build and run the Dockerfile following your resource allocation instructions.

```bash
podman build -t rstudioconda --build-arg USERNAME=XXXX --build-arg PASSWORD=XXXX -f Dockerfile .
```



## 7. Running RStudio and JupyterLab in Podman

Once you have built and created the image in Podman, follow the steps below to access RStudio and JupyterLab.



```bash
hostname -i
conda activate scpy
podman run -d --name rstudio_container -p 8787:8787 -p 8888:8888 rstudioconda
http://10.192.XX.XX:8787/lab #for Jupyter
http://10.192.XX.XX:8787 # for R-studio

   ```


