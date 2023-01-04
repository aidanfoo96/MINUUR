# Inputting Files To MINUUR
Inputting your files to MINUUR can be done in one of two ways: 

##### Automatic 
You can input a paired fastq file into the `resources` directory of MINUUR. Simply `mv` or `cp` some raw paired gzipped fastq files into this directory. The caveat, however, is this uses a strict naming scheme, where your files must be named in the following format:

`{samplename}_1.fastq.gz` and `{samplename}_2.fastq.gz`

Your `{samplename}` must then be specified in the sample_list.tsv file located in the `config` directory and look a little something like this: 


```{figure} pics/Sample_table_texample_2.png

```

##### Manual 
Maybe you have an alternative naming scheme or files deposited in a directory you don't want them moved from? MINUUR offers a more flexible approach.

Here, input your sample name and relative path to your data in the `samples_table.tsv` file located in the `config` directory. Similar to above, please also specify your sample names in the `sample_list.tsv` file. Both files should look a little like this: 

```{figure} pics/sample_list_example.png
```

### Putting It All Together
By now you should have downloaded the taxonomic databases you might be interested in, configured your data inputs using either Manual or Automatic Options and placed your host reference genome someone (maybe in the resources directory of MINUUR), but that doesn't really matter. 

If you think all is well, you can try running `snakemake -np` from the `working` directory of MINUUR. This will perform a dry run of your analysis and print a load of inputs and outputs that you're intending to run. Have a look at some of the steps and you can see what might be happening under the hood. Cool. However, if you tried running this for real, you'd likely run into a string of errors. That's because you haven't *configured* your workflow yet. We need to specify some paths towards certain databases and steps you want to run / don't want to run and the behaviour of certain parameters. Juicy. Lets talk about this next. 