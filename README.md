![marple-diagnostics logo](https://marple-diagnostics.org/wp-content/uploads/2022/07/MARPLE-logo-1024x364.png)

To install the dependencies and create the environment:
```
. marple-install.sh
```
> This also adds the functions to run the pipeline in `.bashrc`.

To concatenate and transfer the fastq files created in MinKnow for pgt or pst:
```
transfer-pgt <name of experiment> <barcodeXX> <name of sample>
transfer-pst <name of experiment> <barcodeXX> <name of sample>
```
To run marple:
```
marple
```
> This can be executed from anywhere in the system.

> [!IMPORTANT]
> These functions assume the reads are in `/var/lib/minknow/data/reads/` and marple in `/home/$USER/`.

> [!TIP]
> You can also access just the environment to run other commands in `/home/$USER/marple`, using:
> ```
> marple --dev
> ```

> [!CAUTION]
> This pipeline is a work in progress, until it gets a numbered version, you may also consider using
> - https://github.com/SaundersLab/marple-pst for Pst
> - https://github.com/SaundersLab/marple-pgt for Pgt
