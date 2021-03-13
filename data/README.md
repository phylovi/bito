# Random notes about data files

## `hello_rooted.nwk` and `hello_rooted_two_trees.nwk`

    iqtree -m JC -s hello.fasta

gives

    (mars:0.2072560544,saturn:0.0694244266,jupiter:0.0694247559);

which when we reroot in two places gives the trees in the files.


## 7-taxon-slice-of-ds1.fasta

Made by

```
sc --cut 500:1000 ds1.fasta data/7-taxon-slice-of-ds1.fasta
```

Then truncating to 7 taxa and renumbering things to be 0-indexed.
