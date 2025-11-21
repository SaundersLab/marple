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

Cite this work: \
[1] [A portable, nanopore-based genotyping platform for near real-time detection of Puccinia graminis f. sp. tritici lineages and fungicide sensitivity](https://doi.org/10.1186/s12864-025-11428-w) \
[2] [MARPLE, a point-of-care, strain-level disease diagnostics and surveillance tool for complex fungal pathogens](https://doi.org/10.1186/s12915-019-0684-y)
