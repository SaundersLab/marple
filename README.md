![marple-diagnostics logo](https://marple-diagnostics.org/wp-content/uploads/2022/07/MARPLE-logo-1024x364.png)

MARPLE must be in `~/marple`, or the install script will stop. Clone it there:
```
git clone https://github.com/SaundersLab/marple.git ~/marple
cd ~/marple
```
> A downloaded zip extracts to `marple-main`. Rename the folder to `marple` and move it to your home directory.

The install script needs `curl` and `mamba`, and tells you how to install them if they are missing. `mamba` comes with [Miniforge](https://github.com/conda-forge/miniforge#install).

To install the remaining dependencies and create the environment, run this from `~/marple`:
```
. marple-install.sh
```
The script prints a success message when everything is installed, or a list of what is missing.

> This script also adds the functions to run the pipeline in `.bashrc`.

To concatenate and transfer the fastq files created in MinKnow for pgt or pst:
```
transfer-pgt <name of experiment> <barcodeXX>=<name of sample> <barcodeYY>=<name of sample> ...
transfer-pst <name of experiment> <barcodeXX>=<name of sample> <barcodeYY>=<name of sample> ...
```
To run MARPLE:
```
marple
```

> These functions can be executed from anywhere in the system.

> [!IMPORTANT]
> The functions here, assume the reads are in `/var/lib/minknow/data/` and marple in `/home/$USER/`.

Cite this work: \
[1] https://doi.org/10.1186/s12864-025-11428-w \
[2] https://doi.org/10.1186/s12915-019-0684-y 