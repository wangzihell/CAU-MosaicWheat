#!/usr/bin/env python

import numpy as np
import pybedtools
from optparse import OptionParser


def structure_2bed (bedfile1, bedfile2, BINSIZE):
    x = pybedtools.bedtool.BedTool(bedfile1)
    y = pybedtools.bedtool.BedTool(bedfile2)
    xy = x.cat(y)
    total_bed = xy.merge()

    str_bed_list =[]
    for i in range(0, len(total_bed)):
        interval = total_bed[i]
        CHR = interval.chrom
        START = interval.start
        END = interval.end
        sample_name_x = []
        sample_name_y = []
        for start in range(START, END, BINSIZE):
            tmp_bed = pybedtools.BedTool(CHR+" "+str(start)+" "+str(start+BINSIZE), from_string=True)
            laped_bed_x = x.intersect(tmp_bed)
            laped_bed_y = y.intersect(tmp_bed)
            if len(laped_bed_x) != 0 and len(laped_bed_y) != 0:

                for j in range(0, len(laped_bed_x)):
                    sample_name_x.append(laped_bed_x[j].name)

                for j in range(0, len(laped_bed_y)):
                    sample_name_y.append(laped_bed_y[j].name)

                str_bed_list.append([CHR, start, start+BINSIZE, sample_name_x, sample_name_y])

            sample_name_x = []
            sample_name_y = []

    return str_bed_list


def main():
    usage = "Usage: %prog [-i <input file>] [-o <output file>] [-s <column number of first sample>]\n"
    #
    parser = OptionParser(usage)
    parser.add_option("-d", dest="dist_folder",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", default = 5000000)
    parser.add_option("-l", dest="BINSIZE",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", default = 5000000)
    parser.add_option("-a", dest="bedfile1",
                  help="Input file, use STDIN if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    parser.add_option("-b", dest="bedfile2",
                  help="column id for chromosome [default: %default]")
    parser.add_option("-s", dest="orderfile",
                  help="column id for chromosome [default: %default]")
    parser.add_option("-o", dest="outfile",
                  help="Output file, use STDOUT if omit; "
                       "gzip file end with \".gz\".", metavar="FILE")
    #
    (options, args) = parser.parse_args()
    #
    WD = options.dist_folder

    bed_matrix = structure_2bed(options.bedfile1, options.bedfile2, options.BINSIZE)

    ORDER_IN = open(options.orderfile, 'r')
    OUTFILE = open(options.outfile, "w")
    #
    SM_order = []
    for line in ORDER_IN :
        tokens = line.strip()
        SM_order.append(tokens)

    for i in range(0, len(bed_matrix)):
        CHR = bed_matrix[i][0]
        START = bed_matrix[i][1]
        END = bed_matrix[i][2]
        sample_list_x = bed_matrix[i][3]
        sample_list_y = bed_matrix[i][4]
        data_matrix = np.loadtxt(WD+CHR+"."+str(int(END/options.BINSIZE))+".5M.dist.txt", "i2")
        sample_series_x = [SM_order.index(j) for j in sample_list_x]
        sample_series_y = [SM_order.index(j) for j in sample_list_y]
        for n in range(0, len(sample_series_x)):
            for m in range(0, len(sample_series_y)):
                OUTFILE.write(str(data_matrix[sample_series_x[n], sample_series_y[m]])+"\n")
                # OUTFILE.write("\t".join([CHR, str(START), str(END), str(data_matrix[sample_series_x[n], sample_series_y[m]]), sample_list_x[n], sample_list_y[m]])+"\n")
# ===========================================
if __name__ == "__main__":
    main()
