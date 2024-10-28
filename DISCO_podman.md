## Copy your files from gcloud to /broad/macosko/data/


**Connect to the login node:**
```bash
ssh <username>@login00.broadinstitute.org
use UGER
ish -l h_vmem=16G
```

## First-Time Setup:

 **Create your directory (if you haven't). Transfer Files or Scripts from Cloud:**

   ```bash
mkdir -p /broad/macosko/data/XXX
use Google-Cloud-SDK
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
cd /broad/macosko/data/XXX/
```

## 5. Create Podman

Create a  Podman:

```bash
podman pull ghcr.io/rocker-org/rstudio:latest
podman run -d --name rstudio-container     -e USER=XXXX   -e PASSWORD=XXXX     -e DEFAULT_USER=XXX     -p 8787:8787     docker.io/rocker/rstudio:latest

```



## 6. Running RStudio and JupyterLab in Podman

Once you have built and created the image in Podman, follow the steps below to access RStudio 

The username is rstudio, and the password is rstudio. 
Change password afterwards


```bash
hostname -i

http://10.192.XX.XX:8787 # for R-studio

   ```


