# call INS and DEL ranging in size from larger than 50bp up to 500k using hifi reads with minimap2 and sniffles

```
snakemake -s minimap2sniffles.smk --cores 128 -pn
```

# vg autoindex

```
snakemake -s vgautoindex.smk --cores 128 -pn
```

# vg giraffe 

```
snakemake -s vgGiraffePackCall.smk --cores 128 -pn
```

# paragraph genotyping

```

```
