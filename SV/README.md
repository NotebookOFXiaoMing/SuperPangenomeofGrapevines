### structural variation

construct VG graph

dependencies

```
minimap2
syri
SURVIVOR
bcftools
vg
R package tidyverse ggridges
snakemake
```

run

```
snakemake -s minimap2_syri_plotsr.smk --configfiles config.yaml --cores 52 -p
snakemake -s vgGiraffePackCall.smk --configfiles config02.yaml --cores 52 -p
```
