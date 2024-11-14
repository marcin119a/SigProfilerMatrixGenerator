from SigProfilerMatrixGenerator.scripts.utils import *
import os
from SigProfilerMatrixGenerator.scripts.one_generated import catalogue_generator_single


def read_vcf(vcf_path, sample):
    """
    Wczytaj dane z pliku VCF (bez użycia zewnętrznych bibliotek) i przekształć je na listę zawierającą dane
    chromosomu, pozycję oraz zmiany (ref, alt).

    :param vcf_path: Ścieżka do pliku VCF.
    :return: Lista zawierająca dane chromosomu, pozycji oraz ref i alt dla każdej mutacji.
    """
    lines = []
    with open(vcf_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Pomijanie linii nagłówka

            columns = line.strip().split('\t')
            chrom = columns[0]
            pos = int(columns[1])
            ref = columns[3]
            alt = columns[4].split(',')  # Obsługa wariantów z wieloma ALT

            for alt_allele in alt:
                lines.append([sample, chrom, pos, ref, alt_allele])

    return lines

import pandas as pd
if __name__ == "__main__":
    mutation_pd = {}
    mutation_dict = {}
    print(mut_types)
    samples = ['GRCh37_bench']
    mutation_dict["6144"] = pd.DataFrame([[[] for _ in samples] for _ in mut_types], index=mut_types, columns=samples)
    chrom_list = [str(i) for i in range(1, 23)] + ['X', 'Y']

    for chrom in chrom_list:
        lines = read_vcf('SigProfilerMatrixGenerator/references/vcf_files/GRCh37_bench/GRCh37_bench.vcf', sample=samples[0])
        chrom_lines = [line for line in lines if line[0] == chrom]
        chrom = chrom

        vcf_path = 'SigProfilerMatrixGenerator/references/tests/WGS/GRCh37/output/temp/test_01270cc6-d80d-4b2f-b4d1-066e6a46953b/SNV/'
        chrom_path = 'SigProfilerMatrixGenerator/references/chromosomes/tsb/GRCh37/'
        output_matrix = 'SigProfilerMatrixGenerator/references/tests/WGS/GRCh37/output/'
        exome = False
        mutation_dinuc_pd_all = pd.DataFrame([[[] for _ in samples] for _ in mut_types], index=mut_types, columns=samples)
        functionFlag = True
        bed = False

        transcript_path = 'SigProfilerMatrixGenerator/references/chromosomes/transcripts/GRCh37/'
        seqInfo = False
        gs = False
        log_file = 'log'


        context = "6144"
        (
            mutation_pd,
            skipped_mut,
            total,
            total_DINUC,
            mutation_dinuc_pd_all,
        ) = catalogue_generator_single(
            lines,
            chrom,
            mutation_dict,
            mutation_dinuc_pd_all,
            mutation_types_tsb_context,
            vcf_path,
            chrom_path,
            output_matrix,
            exome,
            functionFlag,
            bed,
            tsb_ref,
            transcript_path,
            seqInfo,
            gs,
            log_file
        )
    mutation_dict["6144"].to_csv('output_mutations.csv')
