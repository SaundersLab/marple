![marple-diagnostics logo](https://marple-diagnostics.org/wp-content/uploads/2022/07/MARPLE-logo-1024x364.png)

To install the dependencies and create the environment:
```
. marple-install.sh
```
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

> [!TIP]
> You can also access just the environment to run other commands in `/home/$USER/marple`, using:
> ```
> marple --dev
> ```
