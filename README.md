# metasite_coverage


This tool generates mean coverage plot around list of sites.

### Requirements:

* Python3
* pip
* samtools


#### Python modules:

* matplotlib
* numpy
* pandas
* pysam

### Installation:

```
git clone https://github.com/BioInfoData/metasite_coverage
cd metasite_coverage
pip install -r requirements.txt

```

### Input files:
1. BAM file
2. List of sites in BED6 format


### Output files:
1. plot of the coverage around the sites +- range (png image)
2. Additional plot focus in range +- 5 around sites (.png image)
3. Table with raw coverage data for each site (.csv file)


### Usage:

```
python metasite_coverage.py

options:                                                                                                                                              
  -h, --help            show this help message and exit                                                                                               
  --bam BAM             path to bam file                                                                                                              
  --sites_file SITES_FILE                                                                                                                             
                        bed file with sites                                                                                                           
  --range RANGE         range to plot around the site (+-range)must be >=5                                                                            
  --sample SAMPLE       sample name
  --min_score MIN_SCORE
                        min score value (column 5 at bed file). Not required
                        ```
  --min_coverage MIN_COVERAGE
                        min coverage value to include in site in analysis. Not required, default is 0
```

### Steps of analysis:
1. For each site calculate coverage around the site (using pysam).
2. Filter out sites with 0 coverage (or other min is given by the user).
3. For each site normalize to max value.
4. Sum all coverage from all sites (normalised) to one vector.
4. plot data +- range around sites.
5. plot second  focus plot in range +- 5 around sites.

