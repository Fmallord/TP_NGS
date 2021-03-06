## program used to calculate sigma squared and Va, and that calculates their ratio afterwards.

#!/usr/bin/env python3

import numpy as np
import argparse
import os


def variance_additive(array): #sum of the snp variances
    n=array.shape[0] #nb ligne = nb de snps
    m=array.shape[1]  #nb colonnes = nb d'individus
    X=0
    for i in range (n):
        X+=np.var(array[i,:])
    return X


def sigma_squared(array):   #variance of the summed individuals
    n=array.shape[0] #nb ligne = nb de snps
    m=array.shape[1]  #nb colonnes = nb d'individus
    X=[]
    for j in range (m):
        X+=[np.sum(array[:,j])]
    return np.var(X)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--vcf', required=True, type=str,
                        dest="v", metavar="<vcf>",
                        help="The relative name of the .vcf file")
    parser.add_argument('-c', '--cutoff', required=False, type=int,
                        dest="c", metavar="<cutoff>", default=1,
                        help="Cut-off for minor allele")
    args = parser.parse_args()

    variant_list = []

    path = os.getcwd()

    vcf_file = open("{0}/{1}".format(path, args.v), 'r')
    for line in vcf_file:
        if line[0] != "#":
            split_line = line.replace('\n', '').split("\t")[9:]
            sparse_line = np.zeros(len(split_line), dtype=np.int8)
            for i, genotype in enumerate(split_line):
                value = genotype.count("1")
                if value > 0:
                    sparse_line[i] = value

            if args.c >= np.sum(sparse_line) > 0:
                variant_list.append(sparse_line)
    vcf_file.close()

    sparse_array = np.array(variant_list, dtype=np.int8)
    nbr_cols = sparse_array.shape[0]
    nbr_rows = sparse_array.shape[1]
    print("{0} SNPs".format(nbr_cols))
    print("{0} individuals".format(nbr_rows))

    tsv_file = open("{0}/{1}.analysis.tsv".format(path, args.v[:-4]), 'w')
    tsv_file.write('\t'.join(['SNP', 'sample_size', 's^2', 'Va', 's^2/Va']) + '\n')

    va = variance_additive(sparse_array)

    print("Va={0}".format(va))

    s2 = sigma_squared(sparse_array)

    print("\sigma^2={0}".format(s2))
    print("\sigma^2/Va={0}".format(s2 / va))
    tsv_file.write('\t'.join([str(nbr_cols), str(nbr_rows), str(s2), str(va), str(s2 / va)]) + '\n')
    tsv_file.close()
    print("Analysis performed")
