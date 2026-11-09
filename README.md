![marple-diagnostics logo](https://marple-diagnostics.org/wp-content/uploads/2022/07/MARPLE-logo-1024x364.png)

To install the dependencies and create the environment:
```
. marple-install.sh
```
> This script also adds the functions to run the pipeline in `.bashrc`.

To concatenate and transfer the fastq files created in MinKnow for pgt or pst, use a samplesheet:
```
transfer-pgt <name of experiment> --samplesheet samplesheet.csv
transfer-pst <name of experiment> --samplesheet samplesheet.csv
```
The samplesheet is a csv file with a header row and one row per sample. Excel can create it, and any other columns are ignored:
```
barcode,sample
barcode01,field1
barcode02,field2
```
For a few samples, you can type the pairs instead:
```
transfer-pgt <name of experiment> <barcodeXX>=<name of sample> <barcodeYY>=<name of sample> ...
transfer-pst <name of experiment> <barcodeXX>=<name of sample> <barcodeYY>=<name of sample> ...
```
Sample names can only contain letters, numbers, `.`, `_` and `-`.

The transfer reads `fastq_pass` for live basecalling and `basecalling/pass` for basecalling after the run. If an experiment has several of these directories, it uses the most recent and lists the others.

> If you installed MARPLE before these options existed, run `. marple-install.sh` again and answer `y` to updating the functions.

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