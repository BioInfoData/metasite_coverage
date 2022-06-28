#import matplotlib.pyplot as plt
#import pysam
import numpy as np
import pandas as pd
import logging
import argparse
import os
import sys


logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
streamHeandeler = logging.StreamHandler()
logger.addHandler(streamHeandeler)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHeandeler.setFormatter(formatter)

""""
Check if input file exists
"""

def check_file(file_name):
    if os.path.isfile(file_name):
        return
    else:
        logger.critical("file {} not found".format(file_name))
        sys.exit(1)


""""
Check if input range is acceptable (>=5)
"""

def check_range(myRange):
    if int(myRange) >=5:
        return
    else:
        logger.critical("Range must be >=5")
        sys.exit(1)


"""
Get pysam.AlignmentFile and positions for one site.
Returns coverage for one site
"""
def one_area_covearge(samfile, chr, start, end, strand):

    coverage = samfile.count_coverage(
                      contig = chr,
                      start = start,
                      stop = end,)

    cov = np.asarray(coverage)
    cov_all = cov[0]
    for i in range(1,len(cov)):
        cov_all = np.add(cov_all, coverage[i])
    if strand == "-":
        cov_all = np.flip(cov_all)

    return(cov_all)


"""
Get bam file and a sites file. 
For each site calculate coverage around the site (using pysam)
Collect data from all sites.

Returns data frame with raw data (counts from all sites in the input file)
"""
def area_coverage_all(bam,sites_file, myRange, score_cutoff = None):
    sites = pd.read_csv(sites_file,sep = "\t", header= None)
    if score_cutoff:
        sites = sites[sites[4] >= int(score_cutoff)]
    samfile = pysam.AlignmentFile(bam, "rb")
    dict_res = {}
    for i,row in sites.iterrows():
        chr = row[0]
        strand = row[5]
        one_site = row[2]
        if strand == "-":
            start = one_site - (myRange +1)
            end = one_site + (myRange -1)
        else:
            start = one_site - myRange
            end = one_site + myRange

        id = chr + ":" + str(one_site) + strand
        dict_res[id] = one_area_covearge(samfile, chr, start, end, strand)
    df_raw = pd.DataFrame.from_dict(dict_res, orient='index')
    samfile.close()
    return(df_raw)

"""
Gets data frame with raw data (counts from all sites in the input file).
filer out lines by min coverage value (default = 0)
"""

def filter_data(df_raw, min_val):
    df_filt = df_raw.loc[df_raw.apply(np.min, axis=1) >= min_val]
    return df_filt
"""
Gets data frame with raw data (counts from all sites in the input file, after filtering low coverage).
normalise each row (=site) to the max value in the row
"""
def normalyse_data(df_filt):
    df_norm = df_filt.apply(lambda x: round(x / np.max(x), 2), axis=1)
    return df_norm

"""
Gets data frame with counts data after filtering low coverage and after normalisation.
sum all rows and: 1. plot the data. 2. write data to csv file
"""

def plot_data(df_norm,sample_name, plot_type, myRange):
    meta_data = df_norm.sum(axis = 0)
    meta_data.index =  list(range(-1*myRange + 1, (myRange+1)))
    meta_data.to_csv("{}_plotData_{}_range{}.csv".format(sample_name, plot_type,myRange))
    if plot_type == "bar":
        # site is norm to 1
        meta_data = meta_data/meta_data.loc[0]
        meta_data.plot.bar(width = 1.0, edgecolor = 'black')
        plt.axvline(x=myRange-1, color="red")
    elif plot_type == "line":
        meta_data.plot()
        plt.axvline(x = 0, color="red")
    plt.ylim(0, np.max(meta_data) + np.max(meta_data)*0.1)
    plt.title(sample_name)
    plt.ylabel("Coverage (normalised)")
    plt.xlabel("Relative position")
    plt.savefig("{}_coverageSite_{}_range{}.png".format(sample_name, plot_type,myRange))

def write_raw_data_to_file(df_raw,myRange,sample_name):
    df_raw.to_csv("{}_coverageSite_range{}.tsv".format(sample_name,myRange), sep = "\t")

"""
Get bam file and a sites file and plot the coverage around the sites +- range

1. For each site calculate coverage around the site (using pysam)
2. Filter out sites with 0 coverage (or other min is given by the user)
3. For each site normalise to max value
4. Sum all coverage from all sites (normalised) to one vector
4. plot data +- range around sites

Returns data frame with normalised data (all sites +- range)
"""

def meta_site_coverage(bam,sites_file, myRange, score_cutoff, sample_name, min_val,write_res):
    # For each site calculate coverage around the site (using pysam)
    df_raw = area_coverage_all(bam,sites_file, myRange, score_cutoff)
    # Filter out sites with 0 coverage (or other min is given by the user)
    df_filt = filter_data(df_raw, min_val)
    # For each site normalise to max value
    df_norm = normalyse_data(df_filt)
    logger.info("Plotting results full range")
    # plot data
    plot_data(df_norm, sample_name, "line", myRange)
    if(write_res):
        write_raw_data_to_file(df_raw,myRange, sample_name)
    return(df_norm)



"""
plot focus plot in range +- 5 around sites
"""
def focus_plot(df_norm,midRange,myRange, sample_name):
    logger.info("Plotting results in range {}".format(midRange))
    df_focus = df_norm.iloc[:, myRange - midRange:myRange + midRange]
    plt.clf()
    plot_data(df_focus, sample_name, "bar", midRange)


"""
Main function.
Get bam file and a sites file:
1. plot the coverage around the sites +- range.
2. plot second  focus plot in range +- 5 around sites.
3. write csv files with data.

Steps of analysis:
1. For each site calculate coverage around the site (using pysam)
2. Filter out sites with 0 coverage (or other min is given by the user)
3. For each site normalise to max value
4. Sum all coverage from all sites (ormalised) to one vector
4. plot data +- range around sites
5. plot second  focus plot in range +- 5 around sites

Output are png files of plots in the running directory + csv files of the data
"""
def meta_site_coverage_main(bam,sites_file, myRange, sample_name,score_cutoff = None,  min_val =0 ):
    # check that input files exists
    check_file(bam)
    check_file(sites_file)
    check_range(myRange)
    logger.info("running  bam file: {}".format(bam))
    # get coverage for each position, filter by min coverage, normalise to max to get the same scale and plot results
    df_norm = meta_site_coverage(bam,sites_file, myRange, score_cutoff, sample_name, min_val, write_res = True)
    # plot a second graph only +-5 bases from the site
    midRange = 5
    focus_plot(df_norm, midRange, myRange, sample_name)
    logger.info("finished successfully")



def parse_user_data():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True,
                        help='path to bam file')
    parser.add_argument('--sites_file', required=True,
                        help='bed file with sites')
    parser.add_argument('--range', required=True,
                        help='range to plot around the site (+-range)'
                             'must be >=5')
    parser.add_argument('--sample', required=True,
                        help='sample name')
    parser.add_argument('--min_score', required=False,
                        help='min score value (column 5 at bed file). '
                             'Not required', default= None)
    parser.add_argument('--min_coverage', required=False,
                        help='min coverage value to include in site in analysis. '
                             'Not required, default is 0', default=0)

    args = parser.parse_args()

    meta_site_coverage_main(args.bam, args.sites_file, int(args.range), args.sample, args.min_score, int(args.min_coverage))

if __name__ == '__main__':
    parse_user_data()






