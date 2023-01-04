## Running MINUUR

To do a dry-run, enter the command 

`snakemake -n`.

To run the workflow on 3 cores with conda installing all required packages as you go, run: 

`snakemake --cores 3 --use-conda`

I would also recommend tacking on the following command: 

`snakemake --cores 3 --use-conda --keep-going` 

A known error I run into is where MINUUR can't find any bins from your assemblies. This is common when trying to recover MAGs using this approach since you may not have the required read depth. The `--keep-going` command ensures the pipeline still runs even when you encounter this error. Just know that subsequent steps after binning (CheckM) won't run since there will be no bins to work on.
