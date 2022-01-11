##### Download Required MINUUR Databases #####
conda create --name install_db -c bioconda metaphlan humann
conda activate install_db

# Install Standard Kraken2 + Bracken Indexed Database
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz -P ../../resources/ 
mkdir ../../resources/kraken_db
mv ../../resources/k2_standard_20210517.tar.gz ../../resources/kraken_db
tar xvzf --verbose ../../resources/kraken_db/k2_standard_20210517.tar.gz

# Install metaphlan_db
metaphlan --install --bowtie2db ../../resources/

# Install humann_databases 
humann_databases --download chocophlan full ../../resources/
humann_databases --download uniref uniref90_diamond ../../resources/
