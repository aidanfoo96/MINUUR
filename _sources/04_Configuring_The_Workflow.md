## Configuring the Pipeline
Navigate to the `config` directory of MINUUR and open the file `config.yaml`. You will come across a list of options of how you want the pipeline to run. Let's go through these options.

###### Pre-Processing

```
samples: ../config/samples_list.tsv

input_type: 
  fastq: True
  bam: False

input_config: 
  automatic: False
  manual: ../config/sample_table.tsv
```
**Samples**

Do not edit this - this will point to your `sample_list.tsv` file in which you'll have listed your sample names.

**input_type**

Choose between fastq or bam. Most likely your input will be a fastq file. Some of the samples I processed were bam files so I needed this option to start the run at this point. Note that if you have a mixture of bam and fastq files, MINUUR will not be able to process both simultanesouly. Please do them one at a time. 

**input_config**

Choose how you want to specify your files. `automatic` or `manual`. To understand the difference between the two, look back at the [previous](03_Inputting_Files.md) section. 

**kraken_db_location**

```
kraken_db_location: /home/db/Kraken2_db/
```

Provide a path to where you've put your indexed Kraken2 database. 

###### Quality Control and Host Separation

**QC**
```
QC: 
  Activate: True
  CutadaptParams: "--minimum-length 50 -q 30"
  CutadaptThreads: 10
```

`Activate:` Insert `True` if you would like to perform QC on your samples, `False` if you would like to skip. 

`CutadaptParams`: Insert specific cutadapt parameters. For my reads I trimmed each read <50 bp and <30 phred score

`CutadaptThreads`: specify the number of threads you'd like to use

**Remove Host From Fastq Gz**

```
RemoveHostFromFastqGz: 
  Threads: 15
  # Give path to directory - allows MINUUR to check there is an existing directory 
  ContaminantDirectory: /home/db/User_dbs/aidan_foo/bt_db/
  # Give path to Bowtie2 indexed host genome
  ContaminantIndex: /home/db/User_dbs/aidan_foo/bt_db/aegypti_db
  AlignmentSensitivty: "--sensitive-local"
```

`Threads`: Specify the number of threads for bowtie2 to use to separate your indexed host reference 

`ContaminantDirectory`: Specify the directory your indexed reference database is stored. This is used as a CheckPoint for the pipeline to confirm you have a host to separate

`ContaminantIndex`: Specify the path to your bowtie2 indexed host to separate. You will have a list of `host.1.bt1, host.2.bt2, host.3.bt2, host.4.bt2, host.rev.1.bt2, host.rev.2.bt2`. When specify this path, do so like `/path/to/your/host` - exclude the extension. 

`AlignmentSensitivity`: specify the sensitivity and alignment type to align your reads. Options for this include

```
- "--very-fast-local"
- "--fast-local"
- "--sensitive-local"
- "--very-sensitive-local" 
- "--very-fast-global" 
- "--fast-global"
- "--sensitive-global"
- "--very-sensitive-global"
```

```
ProcessBam:
  Activate: True
  # Activate if input was Fastq (prior). If Input is bam then FromFastq = False
  FromFastq: True
```

`Activate`: Pick `True` if you would like to generate some statistics of your alignments. 

`FromFastq`: Pick `True` if you input fastq files. If you had bam files choose `False`

###### Read Classification

After QC, you will have some gzipped fastq files pertaining to reads that didn't map to your host in question. Cool! Lets see how much microbial information we can get from this. One path you can take is to classify those reads against the Kraken2 database you downloaded. 


**Kraken Classification**

```
KrakenClassification:
  Activate: True
  Threads: 20
  ConfidenceScore: 0
```

`Activate` : Pick `True` if you would like to do some taxonomic classifications with Kraken2

`Threads`: Choose the number of threads Kraken should use for classification 

`ConfidenceScore`: Choose a number between 0 - 1 which sets how stringent the Kraken2 classifications will be. I picked 0, the default. 1 = very stringent. 

**Kraken Summaries**

```
KrakenSummaries: 
  Activate: True 
  GenusReadThreshold: 6000
  SpeciesReadThreshold: 10000
  StratThreshold: 10000
  KrakenDbStandard: True
```

`Activate`: Use `True` to get some summary data of your Kraken2 classifications. These can be more easily parsed and wrangled in R - especially in the Tidyverse. 

`GenusReadThreshold`: MINUUR will attempt to make some plots of genus level classifications for you. Set a threshold for your plots so they don't appear to crowded. 

`SpeciesReadThreshold`: Same as above but for species level classification.

`StratThreshold`: Again, a read threshold when making stratified plots of all your classifications.

**Extract Kraken Classified Reads**

This option is useful if you are interested in a specific taxon or just want to filter out bacterial / archael / viral reads. I've found this option useful when I wanted to do some metagenome assemblies but found a lot of reads unclassified. This could mean two things 1. There are reads left over from your host or 2. There are potentially mirobes present in your unmapped reads not present in the database. You can choose if you'd like to filter out those bacterial reads and use the resultant paired fastq files for metagenome assembly and binning in the later steps

```
ExtractKrakenTaxa: 
  Activate: True
  taxon_choice: "2"
```

`Activate`: Use `True` to activate this option 

`taxon_choice`: Use "2" to filter all bacteria. This is the taxonomy ID for [bacteria](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=2&lvl=3&lin=f&keep=1&srchmode=1&unlock)

**Bracken Reestimation**

Bracken is used to redistribute Kraken2 classified reads to estimate taxonomic abundance in your sample. This will help remove some of the over/under estimation of classified reads at species and genus level. 

```
BrackenReestimation:
  Activate: True
  BrakenDb: /home/db/Kraken2_db/
  ClassificationLvl: 'S'
  DistributionThresh: 10
  PlotThreshold: 10000
```

`Activate`: Pick `True` if you would like to use Bracken to reestimate reads 

`BrackenDb`: Give a path to your Bracken Database. FYI - if you downloaded the Kraken2 index from the link given [prior](02_Installation.md#installing-a-database), a bracken database is also included here, so specify the same path to your Kraken2 database.

`ClassificationLvl`: choose `S` for species or `G` for genus. Bracken will estimate taxonomic abundance at these levels accordingly.

`PlotThreshold`: MINUUR will try and make some plots for you - set a read threshold so the plots don't appear to crowded. 

**MetaPhlAn3 Classification**

An alternative to Kraken2 is MetaPhlAn3. MetaPhlAn3 uses a marker gene appraoch for read classification. I found MetaPhlAn3 to perform poorly compared to Kraken2 for my mosquito metagenomic data. 

```
MetaphlanClassification:
  Activate: True
  Database: /home/db/User_db/metaphlan_db
  CleanMetaphlanReport: True
  NProc: 8
```

`Activate`: Choose `True` to activate MetaPhlAn

`Database`: Give a path to your metaphlan database

`CleanMetaphlanReport`: Produce a summary file of your MetaPhlAn classifications to easily parse 

`NProc`: Set the number of threads

###### Metagenome Assembly 

The third path MINUUR takes is to try and do some metagenome assemblies, followed by binning of the resultant contigs and metagenome assembled genome QC. This was of most interest to me (and maybe you too!) as it meant I could recover MAGs associated to mosquitoes for further analysis. However...this won't always work. Given that the reads we're trying to assemble were never intended for metagenome assemblies it is possible you'll never get a set of high-quality MAGs. Nevertheless it's still worth trying! And I've had good success using this approach, especailly when the number of classified reads are high enough for a particular taxon. Around ~200,000 worked for me.

**Metagenome Assembly**

```
MetagenomeAssm: 
  Activate: True
  UseKrakenExtracted: False
  Threads: 14
  Memory: 0.5
```

`Activate`: If `True` will perform metagenome assemblies with MegaHit on your paried unmapped fastq files.

`UseKrakenExtracted`: Another way to do your assemblies. If you extracted some taxa of interest, or extracted all bacterial classified reads, you can choose to do this option and perform metagenome assemblies on those reads. I found that assemblies with Kraken2 classified reads result in fewer MAGs but had less contamination and higher-completeness. Check out my [paper](https://www.biorxiv.org/content/biorxiv/early/2022/08/11/2022.08.09.5) for more info. 

**Binning** 

Binning is where we can pull some actual genomes from our assembeld contigs. MINUUR uses MetaBat2 to do this. 

```
MetagenomeBinning: 
  Activate: True
  Threads: 20
  MinimumContigLength: 1500
```

`Activate`: Use `True` if you'd like to bin your genomes 

`Threads`: set the number of threads to use 

`MinimumContigLength`: Don't use contigs below a given threshold. The minimum here is 1500bp. If you go lower than this, you will raise an error. 

**Bin Quality Assurance**

Check the quality of those bins! MINUUR will use CheckM to assess bin quality 

```
CheckmBinQA: 
  Activate: True
  Threads: 15 
```

`Activate`: Use `True` to activate 

`Threads`: specify the number of threads 

Optionally you can also use BUSCO to assess bin quality 

```
BUSCO: 
  Activate: True
```
`Activate`: Use `True` to activate

This can help identify whether low quality MAGs pertain to any eukaryotic type assemblies

###### Woohoo, I've configured my pipeline, now let me run it!

By now you should a configured pipeline, with all the things you want to run and paths specified to all your databases. Lets run it!