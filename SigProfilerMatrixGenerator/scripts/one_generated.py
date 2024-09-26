import os
from MutationMatrixGenerator import gene_range
from scipy import stats
import statsmodels.stats.multitest as sm

def catalogue_generator_single(
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
    log_file,
):
    """
    Generates the mutational matrix for 96, 1536, 384, and 6144 context using a single vcf file with all samples of interest.

    Parameters:
                                    vcf_path  -> path to vcf file of interest
               vcf_path_original  -> path to the original vcf files
                               vcf_files  -> actual vcf file
                       bed_file_path  -> path to the user-provided bed files
                              chrom_path  -> path to chromosome reference files. The chromosomes are saved as strings witht the following
                                                            file name: '1.txt', '2.txt', etc.
                                     project  -> unique name given to the set of samples (ex. 'BRCA')
                       output_matrix  -> path where the final mutational matrix is stored
                                     context  -> desired context (ex. 96, 1536, 384, 6144)
                                       exome  -> flag that generates the matrix in only the exome
                                      genome  -> the desired reference genome
                              ncbi_chrom  -> dictionary that allows for the converstion of ncbi chromosome names to standard format
                                                            for the mm10 assembly.
                            functionFlag  -> flag that is used when calling this function from an alternate script
                                             bed  -> parameter used to filter the mutations on a user-provided BED file
                              bed_ranges  -> dictionary that contains all of the ranges for each chromosome dictated by the user's input BED file
                             chrom_based  -> flag that generates the matrices on a chromosome basis
                                            plot  -> flag that plots the matrices after they are generated
                                    tsb_stat  -> performs a transcriptional strand bias test for the 24, 384, and 6144 contexts. The output is
                                                             saved into the output/TSB directory

                                     tsb_ref  -> dictionary that allows for switching between binary and biologically relevant strings
                     transcript_path  -> path to the transcript files
                                              gs  -> flag that generates a file for the strand bias on a gene basis.
                                    log_file  -> path to the output log file

    Returns:
            If called as a function, returns a nested dictionary of panda data frames for each matrix

    Outputs:
            Outputs a mutational matrix for the desired context within [user-provided input path]/output/[mut_type]/

    """
    out = open(log_file, "a")
    mnv = 5
    # Small functions to provide reverse complements of TSB and sequence info:
    revcompl = lambda x: "".join(
        [
            {
                "A": "T",
                "C": "G",
                "G": "C",
                "T": "A",
                "N": "N",
                "[": "[",
                "]": "]",
                ">": ">",
            }[B]
            for B in x
        ][::-1]
    )
    revbias = lambda x: "".join(
        [
            {
                "0": "0",
                "3": "3",
                "1": "2",
                "2": "1",
                "U": "T",
                "T": "U",
                "B": "B",
                "N": "N",
            }[B]
            for B in x
        ][::-1]
    )

    # Provides the sorting order for the TSB matrices
    bias_sort = {"T": 0, "U": 1, "N": 3, "B": 2}
    tsb = ["T", "U", "N", "B"]
    bases = ["A", "C", "G", "T"]

    # Instantiates all relevant variables
    types = []
    flag = True
    i = 0
    sample_start = None
    gene_counts = {}
    skipped_count = 0
    total_analyzed = 0
    total_analyzed_DINUC = 0
    sequence = ""

    # Instantiates the necessary variables/data structures for DINUCs
    if exome:
        exome_temp_file = "exome_temp.txt"
        exome_file = open(vcf_path + exome_temp_file, "a")
        exome_temp_file_context_tsb_DINUC = "exome_temp_context_tsb_DINUC.txt"
        exome_file_context_tsb_DINUC = open(
            vcf_path + exome_temp_file_context_tsb_DINUC, "a"
        )

    if bed:
        bed_temp_file = "bed_temp.txt"
        bed_file = open(vcf_path + bed_temp_file, "a")
        bed_temp_file_context_tsb_DINUC = "bed_temp_context_tsb_DINUC.txt"
        bed_file_context_tsb_DINUC = open(
            vcf_path + bed_temp_file_context_tsb_DINUC, "a"
        )

    dinuc_tsb_ref = ["CC", "CT", "TC", "TT", "AA", "AG", "GA", "GG"]
    chrom_start = chrom

    # Opens the input vcf file
    # tutaj masz każdy chromosome
    with open(chrom_path + chrom_start + ".txt", "rb") as f2:
        chrom_string = f2.read()

        dinuc_sub = [int(y[2]) - int(x[2]) for x, y in zip(lines, lines[1:])]
        dinuc_index_temp = [i for i, x in enumerate(dinuc_sub) if x == 1]
        mnv_index_temp = [i for i, x in enumerate(dinuc_sub) if x <= mnv and x > 0]

        if seqInfo:
            seqOut_path = output_matrix + "vcf_files/SNV/"
            seqOut_path_dinuc = output_matrix + "vcf_files/DBS/"
            seqOut_path_mns = output_matrix + "vcf_files/MNS/"
            if not os.path.exists(seqOut_path):
                os.makedirs(seqOut_path)
            if not os.path.exists(seqOut_path_dinuc):
                os.makedirs(seqOut_path_dinuc)
            if not os.path.exists(seqOut_path_mns):
                os.makedirs(seqOut_path_mns)
            seqOut = open(seqOut_path + chrom_start + "_seqinfo.txt", "w")
            seqOut_dinuc = open(seqOut_path_dinuc + chrom_start + "_seqinfo.txt", "w")
            seqOut_mns = open(seqOut_path_mns + chrom_start + "_seqinfo.txt", "w")

            mnv_index = []
            if len(mnv_index_temp) > 1:
                if mnv_index_temp[0] + 1 == mnv_index_temp[1]:
                    if lines[mnv_index_temp[0]][0] == lines[mnv_index_temp[1]][0]:
                        mnv_index.append(mnv_index_temp[0])
                elif mnv_index_temp[0] not in dinuc_index_temp:
                    mnv_index.append(mnv_index_temp[0])
                else:
                    pass
                for i in range(1, len(mnv_index_temp) - 1, 1):
                    if mnv_index_temp[i] + 1 == mnv_index_temp[i + 1]:
                        if (
                            lines[mnv_index_temp[i]][0]
                            == lines[mnv_index_temp[i + 1]][0]
                        ):
                            mnv_index.append(mnv_index_temp[i])
                    elif mnv_index_temp[i] - 1 == mnv_index_temp[i - 1]:
                        if (
                            lines[mnv_index_temp[i]][0]
                            == lines[mnv_index_temp[i - 1]][0]
                        ):
                            mnv_index.append(mnv_index_temp[i])
                    elif mnv_index_temp[i] not in dinuc_index_temp:
                        mnv_index.append(mnv_index_temp[i])
                if mnv_index_temp[-1] - 1 == mnv_index_temp[-2]:
                    if lines[mnv_index_temp[-1]][0] == lines[mnv_index_temp[-2]][0]:
                        mnv_index.append(mnv_index_temp[-1])
                elif mnv_index_temp[-1] not in dinuc_index_temp:
                    mnv_index.append(mnv_index_temp[-1])

                dinuc_index = [x for x in dinuc_index_temp if x not in mnv_index]

            else:
                dinuc_index = dinuc_index_temp
                mnv_index = [x for x in mnv_index_temp if x not in dinuc_index]

            if mnv_index:
                i = 0
                mnv_seq = ""
                mut_seq = ""
                save_start = None
                change_start = True
                while i < len(mnv_index):
                    line1 = lines[mnv_index[i]]
                    line2 = lines[mnv_index[i] + 1]
                    start1 = int(line1[2])
                    start2 = int(line2[2])
                    previous_ref = line1[3]
                    ref = line2[3]
                    previous_mut = line1[4]
                    mut = line2[4]
                    sample = line1[0]
                    sample2 = line2[0]

                    if sample != sample2:
                        mnv_seq = ""
                        mut_seq = ""

                    else:
                        if change_start:
                            save_start = start1
                            change_start = False
                        mnv_seq += previous_ref
                        mut_seq += previous_mut

                        for l in range(start1 + 1, start2, 1):
                            mnv_seq += tsb_ref[chrom_string[l - 1]][1]
                            mut_seq += tsb_ref[chrom_string[l - 1]][1]

                        if i < len(mnv_index) - 1:
                            if mnv_index[i] + 1 != mnv_index[i + 1]:
                                mnv_seq += tsb_ref[chrom_string[start2 - 1]][1]
                                mut_seq += mut
                                print(
                                    "\t".join(
                                        [
                                            sample,
                                            chrom,
                                            str(save_start),
                                            mnv_seq + ">" + mut_seq,
                                        ]
                                    ),
                                    file=seqOut_mns,
                                )
                                mnv_seq = ""
                                mut_seq = ""
                                change_start = True
                        elif i == len(mnv_index) - 1:
                            mnv_seq += tsb_ref[chrom_string[start2 - 1]][1]
                            mut_seq += mut
                            print(
                                "\t".join(
                                    [
                                        sample,
                                        chrom,
                                        str(save_start),
                                        mnv_seq + ">" + mut_seq,
                                    ]
                                ),
                                file=seqOut_mns,
                            )
                    i += 1

        else:
            dinuc_index = dinuc_index_temp

        for x in dinuc_index:
            strand = "1"
            line1 = lines[x]
            line2 = lines[x + 1]
            previous_ref = line1[3]
            ref = line2[3]
            previous_mut = line1[4]
            mut = line2[4]
            dinuc = "".join([previous_ref, ref, ">", previous_mut, mut])
            sample = line1[0]
            sample2 = line2[0]
            start = int(line1[2])


            if sample != sample2:
                continue

            try:
                dinuc_seq = "".join(
                    [
                        tsb_ref[chrom_string[start - 2]][1],
                        "[",
                        dinuc,
                        "]",
                        tsb_ref[chrom_string[start + 1]][1],
                    ]
                )
                bias = tsb_ref[chrom_string[start - 1]][0]
                if "N" in dinuc_seq:
                    print(
                        "The position is out of range. Skipping this mutation: "
                        + chrom
                        + " "
                        + str(start)
                        + " "
                        + ref
                        + " "
                        + mut,
                        file=out,
                    )
                    skipped_count += 1
                    continue

            except:
                print(
                    "The position is out of range. Skipping this mutation: "
                    + chrom
                    + " "
                    + str(start)
                    + " "
                    + ref
                    + " "
                    + mut,
                    file=out,
                )
                skipped_count += 1
                continue

            dinuc_seq_tsb = bias + ":" + dinuc_seq
            if dinuc[:2] not in dinuc_tsb_ref:
                dinuc_seq_tsb = "Q" + dinuc_seq_tsb[1:]
                strand = "0"

                if dinuc_seq_tsb in mutation_types_tsb_context:
                    mutation_dinuc_pd_all.at[dinuc_seq_tsb, sample].extend(f'{chrom_start}')

                else:
                    dinuc_seq_tsb = "".join(
                        [
                            "Q:",
                            revcompl(dinuc_seq_tsb[-1]),
                            "[",
                            revcompl(dinuc_seq_tsb[4:6]),
                            ">",
                            revcompl(dinuc_seq_tsb[7:9]),
                            "]",
                            revcompl(dinuc_seq_tsb[2]),
                        ]
                    )
                    mutation_dinuc_pd_all.at[dinuc_seq_tsb, sample].extend(f'{chrom_start}')

            else:
                if dinuc_seq_tsb in mutation_types_tsb_context:
                    mutation_dinuc_pd_all.at[dinuc_seq_tsb, sample].extend(f'{chrom_start}')

                else:
                    strand = "-1"
                    dinuc_seq_tsb = "".join(
                        [
                            revbias(dinuc_seq_tsb[0]),
                            ":",
                            revcompl(dinuc_seq_tsb[-1]),
                            "[",
                            revcompl(dinuc_seq_tsb[4:6]),
                            ">",
                            revcompl(dinuc_seq_tsb[7:9]),
                            "]",
                            revcompl(dinuc_seq_tsb[2]),
                        ]
                    )
                    mutation_dinuc_pd_all.at[dinuc_seq_tsb, sample].extend(f'{chrom_start}')

            if seqInfo:
                print(
                    "\t".join([sample, chrom, str(start), dinuc_seq_tsb, strand]),
                    file=seqOut_dinuc,
                )

            # Saves the DINUC into temporary files for exome sorting
            if exome:
                exome_file_context_tsb_DINUC.write(
                    sample
                    + "\t"
                    + chrom
                    + "\t"
                    + str(start)
                    + "\t"
                    + dinuc_seq_tsb
                    + "\t"
                    + ref
                    + "\t"
                    + mut
                    + "\n"
                )

            # Saves the DINUC into temporary files for region sorting
            if bed:
                bed_file_context_tsb_DINUC.write(
                    sample
                    + "\t"
                    + chrom
                    + "\t"
                    + str(start)
                    + "\t"
                    + dinuc_seq_tsb
                    + "\t"
                    + ref
                    + "\t"
                    + mut
                    + "\n"
                )

            total_analyzed_DINUC += 1

        for line in lines:
            range_index = 0
            if gs:
                out = open(output_matrix + "gene_strand_bias_counts_SNV.txt", "w")
                out_hot = open(
                    output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt", "w"
                )
                (
                    gene_ranges,
                    gene_counts,
                    gene_names,
                    sample_mut_counts_per_gene,
                    sample_mut_counts_per_mut_type,
                ) = gene_range(transcript_path)

            try:
                sample = line[0]
                start = int(line[2])
                ref = line[3][0].upper()
                mut = line[4][0].upper()

                # Pulls out the relevant sequence depending on the context
                try:
                    sequence = "".join(
                        [
                            tsb_ref[chrom_string[start - 3]][1],
                            tsb_ref[chrom_string[start - 2]][1],
                            tsb_ref[chrom_string[start - 1]][1],
                            tsb_ref[chrom_string[start]][1],
                            tsb_ref[chrom_string[start + 1]][1],
                        ]
                    )
                except:
                    print(
                        "The position is out of range. Skipping this mutation: "
                        + chrom
                        + " "
                        + str(start)
                        + " "
                        + ref
                        + " "
                        + mut,
                        file=out,
                    )
                    out.flush()
                    skipped_count += 1
                    continue

                bias = tsb_ref[chrom_string[start]][0]
                char = sequence[int(len(sequence) / 2)]

                # Prints the sequence and position if the pulled sequence doesn't match
                # the variant from the file
                if char != ref:  # and revcompl(char) != ref:
                    print(
                        "The reference base does not match the reference genome. Skipping this mutation: "
                        + chrom
                        + " "
                        + str(start)
                        + " "
                        + ref
                        + " "
                        + mut,
                        file=out,
                    )
                    out.flush()
                    skipped_count += 1
                    continue

                # Saves the sequence/mutation type if it matched the reference/reverse strand
                else:
                    strand = "1"
                    if ref == "A" or ref == "G":
                        strand = "-1"
                        bias = revbias(bias)
                        ref = revcompl(ref)
                        mut = revcompl(mut)
                        sequence = revcompl(sequence)

                    # Performs the gene strand bias test if desired
                    if gs:
                        for ranges in gene_ranges[chrom_start][range_index:]:
                            dict_key = ref + ">" + mut
                            if start < ranges[0]:
                                break
                            if ranges[0] <= start <= ranges[1]:
                                gene_index = gene_ranges[chrom_start].index(ranges)
                                gene = gene_names[chrom_start][gene_index]
                                if int(strand) + int(ranges[2]) == 0:
                                    dict_key = "T:" + dict_key
                                    gene_counts[gene][dict_key] += 1
                                    if sample not in gene_counts[gene]["samples"]:
                                        gene_counts[gene]["samples"].append(sample)
                                        sample_mut_counts_per_gene[gene][sample] = 1
                                        sample_mut_counts_per_mut_type[gene][sample] = {
                                            "T:C>A": 0,
                                            "T:C>G": 0,
                                            "T:C>T": 0,
                                            "T:T>A": 0,
                                            "T:T>C": 0,
                                            "T:T>G": 0,
                                            "U:C>A": 0,
                                            "U:C>G": 0,
                                            "U:C>T": 0,
                                            "U:T>A": 0,
                                            "U:T>C": 0,
                                            "U:T>G": 0,
                                        }
                                        sample_mut_counts_per_mut_type[gene][sample][
                                            dict_key
                                        ] += 1

                                    else:
                                        sample_mut_counts_per_gene[gene][sample] += 1
                                        sample_mut_counts_per_mut_type[gene][sample][
                                            dict_key
                                        ] += 1
                                elif strand == ranges[2]:
                                    dict_key = "U:" + dict_key
                                    gene_counts[gene][dict_key] += 1
                                    if sample not in gene_counts[gene]["samples"]:
                                        gene_counts[gene]["samples"].append(sample)
                                        sample_mut_counts_per_gene[gene][sample] = 1
                                        sample_mut_counts_per_mut_type[gene][sample] = {
                                            "T:C>A": 0,
                                            "T:C>G": 0,
                                            "T:C>T": 0,
                                            "T:T>A": 0,
                                            "T:T>C": 0,
                                            "T:T>G": 0,
                                            "U:C>A": 0,
                                            "U:C>G": 0,
                                            "U:C>T": 0,
                                            "U:T>A": 0,
                                            "U:T>C": 0,
                                            "U:T>G": 0,
                                        }
                                        sample_mut_counts_per_mut_type[gene][sample][
                                            dict_key
                                        ] += 1
                                    else:
                                        sample_mut_counts_per_gene[gene][sample] += 1
                                        sample_mut_counts_per_mut_type[gene][sample][
                                            dict_key
                                        ] += 1

                    # Saves the mutation key for the current variant
                    mut_key = "".join(
                        [
                            bias,
                            ":",
                            sequence[0 : int(len(sequence) / 2)],
                            "[",
                            ref,
                            ">",
                            mut,
                            "]",
                            sequence[int(len(sequence) / 2 + 1) :],
                        ]
                    )
                    mutation_dict["6144"].at[mut_key, sample].extend(f'{chrom_start}')
                    total_analyzed += 1

                    # If exome is specified, it will write the variant to a temporary exome file.
                    if exome:
                        exome_file.write(
                            sample
                            + "\t"
                            + chrom
                            + "\t"
                            + str(start)
                            + "\t"
                            + mut_key
                            + "\t"
                            + ref
                            + "\t"
                            + mut
                            + "\n"
                        )
                    if bed:
                        bed_file.write(
                            sample
                            + "\t"
                            + chrom
                            + "\t"
                            + str(start)
                            + "\t"
                            + mut_key
                            + "\t"
                            + ref
                            + "\t"
                            + mut
                            + "\n"
                        )
                    if seqInfo:
                        print(
                            "\t".join([sample, chrom, str(start), mut_key, strand]),
                            file=seqOut,
                        )

            except:
                print(
                    "There appears to be an error in this line. Skipping this mutation: "
                    + chrom
                    + " "
                    + str(start)
                    + " "
                    + ref
                    + " "
                    + mut,
                    file=out,
                )
                out.flush()
                skipped_count += 1

            # Once all variants are accounted for, complete the gene strand bias test/output to the final file.
            if gs:
                pvals = []
                qvals = []
                pvals_hot = []
                qvals_hot = []
                hotspots = {}
                gene_bias = []
                gene_bias_hotspots = []
                for gene in gene_counts:
                    if gene not in gene_bias:
                        gene_bias.append(gene)
                    total_count = sum(sample_mut_counts_per_gene[gene].values())
                    for sample in sample_mut_counts_per_gene[gene]:
                        mut_count = sample_mut_counts_per_gene[gene][sample]
                        if mut_count > 3:  # and mut_count/total_count > 0.5:
                            if gene not in gene_bias_hotspots:
                                gene_bias_hotspots.append(gene)
                            if gene not in hotspots:
                                hotspots[gene] = {}
                                for mut, count in sample_mut_counts_per_mut_type[gene][
                                    sample
                                ].items():
                                    hotspots[gene][mut] = count
                                hotspots[gene]["samples"] = [sample]
                                for mut, count in sample_mut_counts_per_mut_type[gene][
                                    sample
                                ].items():
                                    gene_counts[gene][mut] -= count
                                gene_counts[gene]["samples"].remove(sample)
                            else:
                                for mut, count in sample_mut_counts_per_mut_type[gene][
                                    sample
                                ].items():
                                    hotspots[gene][mut] += count
                                    gene_counts[gene][mut] -= count
                                gene_counts[gene]["samples"].remove(sample)
                                hotspots[gene]["samples"].append(sample)

                    sum_tran = 0
                    sum_untran = 0
                    for mut, counts in gene_counts[gene].items():
                        if mut[0] == "T":
                            sum_tran += counts
                        elif mut[0] == "U":
                            sum_untran += counts
                    n = sum_tran + sum_untran
                    if n > 0:
                        pval = stats.binomtest(k=sum_tran, n=n).pvalue
                    else:
                        pval = 1
                    pvals.append(pval)

                    sum_tran_hot = 0
                    sum_untran_hot = 0
                    if gene in hotspots:
                        for mut, counts in hotspots[gene].items():
                            if mut[0] == "T":
                                sum_tran_hot += counts
                            elif mut[0] == "U":
                                sum_untran_hot += counts
                    n_hot = sum_tran_hot + sum_untran_hot
                    if n_hot > 0:
                        pval_hot = stats.binomtest(k=sum_tran_hot, n=n_hot).pvalue
                    else:
                        pval_hot = 1
                    pvals_hot.append(pval_hot)

                qvals = sm.fdrcorrection(pvals)[1]
                qvals_hot = sm.fdrcorrection(pvals_hot)[1]
                ind = pvals.index("BMP7")
                ind2 = pvals_hot.index("BMP7")

                gene_ind = 0
                for gene in gene_bias:
                    gene_counts[gene]["samples"] = len(gene_counts[gene]["samples"])
                    print(gene, end="", file=out, flush=False)
                    sum_tran = 0
                    sum_untran = 0
                    for mut, counts in gene_counts[gene].items():
                        if mut[0] == "T":
                            sum_tran += counts
                        elif mut[0] == "U":
                            sum_untran += counts
                        print("\t" + str(counts), end="", file=out, flush=False)
                    print(
                        "\t"
                        + str(sum_tran)
                        + "\t"
                        + str(sum_untran)
                        + "\t"
                        + str(qvals[gene_ind]),
                        flush=False,
                        file=out,
                    )
                    gene_ind += 1
                out.close()
                with open(output_matrix + "gene_strand_bias_counts_SNV.txt") as f2:
                    lines = [line.strip().split() for line in f2]
                output = open(output_matrix + "gene_strand_bias_counts_SNV.txt", "w")
                print(
                    "GENE\tT:C>A\tT:C>G\tT:C>T\tT:T>A\tT:T>C\tT:T>G\tU:C>A\tU:C>G\tU:C>T\tU:T>A\tU:T>C\tU:T>G\tSampleCount\tTranscribed_total\tUntranscribedTotal\tq_value",
                    file=output,
                )
                for line in sorted(lines, key=lambda x: (float(x[-1])), reverse=False):
                    print("\t".join(line), file=output)
                output.close()

                # Gene strand bias test for hot spot samples.
                gene_ind = 0
                for gene in gene_bias_hotspots:
                    hotspots[gene]["samples"] = len(hotspots[gene]["samples"])
                    print(gene, end="", file=out_hot, flush=False)
                    sum_tran_hot = 0
                    sum_untran_hot = 0
                    for mut, counts in hotspots[gene].items():
                        if mut[0] == "T":
                            sum_tran_hot += counts
                        elif mut[0] == "U":
                            sum_untran_hot += counts
                        print("\t" + str(counts), end="", file=out_hot, flush=False)
                    print(
                        "\t"
                        + str(sum_tran_hot)
                        + "\t"
                        + str(sum_untran_hot)
                        + "\t"
                        + str(qvals_hot[gene_ind]),
                        flush=False,
                        file=out_hot,
                    )
                    gene_ind += 1
                out_hot.close()
                with open(
                    output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt"
                ) as f2:
                    lines = [line.strip().split() for line in f2]
                output = open(
                    output_matrix + "gene_strand_bias_counts_hotspots_SNV.txt", "w"
                )
                print(
                    "GENE\tT:C>A\tT:C>G\tT:C>T\tT:T>A\tT:T>C\tT:T>G\tU:C>A\tU:C>G\tU:C>T\tU:T>A\tU:T>C\tU:T>G\tSampleCount\tTranscribed_total\tUntranscribedTotal\tq_value",
                    file=output,
                )
                for line in sorted(lines, key=lambda x: (float(x[-1])), reverse=False):
                    print("\t".join(line), file=output)
                output.close()

    # Organizes the required dictionaries for the final matrix generation.
    print("Chromosome " + chrom_start + " done", file=out)
    out.flush()

    if seqInfo:
        seqOut.close()
        seqOut_dinuc.close()

    if exome:
        exome_file.close()
    if bed:
        bed_file.close()
    # Calls the function to generate the final matrix
    if functionFlag:
        chrom_start = None
        chrom_start = None
        out.close()

        return (
            mutation_dict,
            skipped_count,
            total_analyzed,
            total_analyzed_DINUC,
            mutation_dinuc_pd_all,
        )

    out.close()

from utils import *



import pandas as pd

if __name__ == "__main__":
    lines = [['GRCh37_bench', '6', '1424677', 'C', 'T'], ['GRCh37_bench', '7', '1424677', 'C', 'T'], ['GRCh37_bench', '6', '1566963', 'C', 'T'], ['GRCh37_bench', '6', '2222541', 'G', 'A'], ['GRCh37_bench', '6', '5085791', 'C', 'T'], ['GRCh37_bench', '6', '6380833', 'C', 'T'], ['GRCh37_bench', '6', '8650024', 'A', 'C'], ['GRCh37_bench', '6', '10971947', 'C', 'G'], ['GRCh37_bench', '6', '11223198', 'C', 'T'], ['GRCh37_bench', '6', '12336416', 'A', 'G'], ['GRCh37_bench', '6', '12928920', 'T', 'C'], ['GRCh37_bench', '6', '13494076', 'T', 'C'], ['GRCh37_bench', '6', '14373449', 'C', 'A'], ['GRCh37_bench', '6', '15331808', 'G', 'A'], ['GRCh37_bench', '6', '15898454', 'A', 'T'], ['GRCh37_bench', '6', '20385749', 'G', 'T'], ['GRCh37_bench', '6', '20771293', 'A', 'T'], ['GRCh37_bench', '6', '20793220', 'G', 'C'], ['GRCh37_bench', '6', '23100985', 'G', 'A'], ['GRCh37_bench', '6', '23452006', 'C', 'T'], ['GRCh37_bench', '6', '23481169', 'G', 'T'], ['GRCh37_bench', '6', '23539997', 'T', 'C'], ['GRCh37_bench', '6', '24710650', 'G', 'A'], ['GRCh37_bench', '6', '24989682', 'C', 'T'], ['GRCh37_bench', '6', '25098946', 'C', 'T'], ['GRCh37_bench', '6', '26598110', 'G', 'A'], ['GRCh37_bench', '6', '26733134', 'G', 'A'], ['GRCh37_bench', '6', '28729315', 'G', 'C'], ['GRCh37_bench', '6', '29148965', 'T', 'A'], ['GRCh37_bench', '6', '29899045', 'T', 'C'], ['GRCh37_bench', '6', '30411401', 'T', 'C'], ['GRCh37_bench', '6', '33231045', 'A', 'C'], ['GRCh37_bench', '6', '35117416', 'C', 'T'], ['GRCh37_bench', '6', '36096542', 'G', 'A'], ['GRCh37_bench', '6', '38479776', 'T', 'C'], ['GRCh37_bench', '6', '39108123', 'G', 'C'], ['GRCh37_bench', '6', '42241105', 'G', 'T'], ['GRCh37_bench', '6', '44154833', 'G', 'C'], ['GRCh37_bench', '6', '44906435', 'T', 'A'], ['GRCh37_bench', '6', '45731245', 'T', 'C'], ['GRCh37_bench', '6', '45894455', 'C', 'T'], ['GRCh37_bench', '6', '47297438', 'G', 'T'], ['GRCh37_bench', '6', '47432353', 'A', 'C'], ['GRCh37_bench', '6', '47454303', 'G', 'A'], ['GRCh37_bench', '6', '47601051', 'A', 'C'], ['GRCh37_bench', '6', '48839562', 'G', 'T'], ['GRCh37_bench', '6', '49309438', 'G', 'T'], ['GRCh37_bench', '6', '50769775', 'C', 'T'], ['GRCh37_bench', '6', '51234636', 'C', 'T'], ['GRCh37_bench', '6', '54461223', 'T', 'A'], ['GRCh37_bench', '6', '57232114', 'A', 'C'], ['GRCh37_bench', '6', '62361201', 'T', 'G'], ['GRCh37_bench', '6', '63472499', 'G', 'A'], ['GRCh37_bench', '6', '64488078', 'A', 'G'], ['GRCh37_bench', '6', '66643616', 'A', 'T'], ['GRCh37_bench', '6', '69041761', 'G', 'C'], ['GRCh37_bench', '6', '70590745', 'C', 'T'], ['GRCh37_bench', '6', '71790610', 'T', 'C'], ['GRCh37_bench', '6', '72548442', 'G', 'A'], ['GRCh37_bench', '6', '73592597', 'C', 'T'], ['GRCh37_bench', '6', '74231441', 'A', 'T'], ['GRCh37_bench', '6', '75976167', 'T', 'C'], ['GRCh37_bench', '6', '76803662', 'G', 'C'], ['GRCh37_bench', '6', '78859462', 'T', 'C'], ['GRCh37_bench', '6', '82964214', 'C', 'T'], ['GRCh37_bench', '6', '83037520', 'C', 'T'], ['GRCh37_bench', '6', '83672705', 'G', 'A'], ['GRCh37_bench', '6', '84217067', 'A', 'C'], ['GRCh37_bench', '6', '84683378', 'A', 'C'], ['GRCh37_bench', '6', '85492515', 'G', 'T'], ['GRCh37_bench', '6', '86331800', 'G', 'A'], ['GRCh37_bench', '6', '86553739', 'C', 'A'], ['GRCh37_bench', '6', '88515530', 'C', 'T'], ['GRCh37_bench', '6', '89181059', 'C', 'T'], ['GRCh37_bench', '6', '89357052', 'T', 'A'], ['GRCh37_bench', '6', '89927824', 'C', 'A'], ['GRCh37_bench', '6', '90672658', 'G', 'T'], ['GRCh37_bench', '6', '93784840', 'T', 'C'], ['GRCh37_bench', '6', '94266218', 'T', 'C'], ['GRCh37_bench', '6', '95512534', 'G', 'T'], ['GRCh37_bench', '6', '95619388', 'C', 'T'], ['GRCh37_bench', '6', '96582320', 'A', 'G'], ['GRCh37_bench', '6', '99005884', 'G', 'C'], ['GRCh37_bench', '6', '100639519', 'G', 'C'], ['GRCh37_bench', '6', '101461619', 'C', 'A'], ['GRCh37_bench', '6', '104348931', 'A', 'G'], ['GRCh37_bench', '6', '104621981', 'A', 'C'], ['GRCh37_bench', '6', '104799872', 'G', 'T'], ['GRCh37_bench', '6', '105670440', 'G', 'A'], ['GRCh37_bench', '6', '105814291', 'T', 'A'], ['GRCh37_bench', '6', '106763538', 'A', 'T'], ['GRCh37_bench', '6', '108085369', 'C', 'G'], ['GRCh37_bench', '6', '109275827', 'C', 'T'], ['GRCh37_bench', '6', '111495092', 'C', 'G'], ['GRCh37_bench', '6', '111786718', 'C', 'T'], ['GRCh37_bench', '6', '113181941', 'C', 'G'], ['GRCh37_bench', '6', '113585649', 'A', 'G'], ['GRCh37_bench', '6', '114611981', 'G', 'C'], ['GRCh37_bench', '6', '114805944', 'G', 'C'], ['GRCh37_bench', '6', '116244699', 'T', 'C'], ['GRCh37_bench', '6', '116888880', 'A', 'T'], ['GRCh37_bench', '6', '116930038', 'G', 'C'], ['GRCh37_bench', '6', '119181329', 'C', 'A'], ['GRCh37_bench', '6', '119812194', 'A', 'G'], ['GRCh37_bench', '6', '120173588', 'G', 'A'], ['GRCh37_bench', '6', '124588315', 'C', 'A'], ['GRCh37_bench', '6', '129603734', 'G', 'T'], ['GRCh37_bench', '6', '129636671', 'G', 'A'], ['GRCh37_bench', '6', '131564152', 'C', 'G'], ['GRCh37_bench', '6', '131980500', 'T', 'C'], ['GRCh37_bench', '6', '133053189', 'G', 'A'], ['GRCh37_bench', '6', '133074924', 'C', 'G'], ['GRCh37_bench', '6', '133583570', 'C', 'G'], ['GRCh37_bench', '6', '135475450', 'G', 'T'], ['GRCh37_bench', '6', '136338301', 'G', 'C'], ['GRCh37_bench', '6', '137823689', 'A', 'G'], ['GRCh37_bench', '6', '138391127', 'T', 'A'], ['GRCh37_bench', '6', '139476382', 'G', 'T'], ['GRCh37_bench', '6', '139505793', 'G', 'C'], ['GRCh37_bench', '6', '143228370', 'T', 'C'], ['GRCh37_bench', '6', '143538639', 'G', 'A'], ['GRCh37_bench', '6', '144127480', 'A', 'G'], ['GRCh37_bench', '6', '146239442', 'C', 'G'], ['GRCh37_bench', '6', '148087706', 'G', 'A'], ['GRCh37_bench', '6', '148368887', 'G', 'C'], ['GRCh37_bench', '6', '149062107', 'G', 'T'], ['GRCh37_bench', '6', '149298849', 'C', 'T'], ['GRCh37_bench', '6', '150076322', 'G', 'A'], ['GRCh37_bench', '6', '150402465', 'G', 'T'], ['GRCh37_bench', '6', '152386309', 'G', 'A'], ['GRCh37_bench', '6', '154310740', 'G', 'A'], ['GRCh37_bench', '6', '155986300', 'G', 'T'], ['GRCh37_bench', '6', '156095629', 'G', 'A'], ['GRCh37_bench', '6', '156777418', 'T', 'C'], ['GRCh37_bench', '6', '156874159', 'G', 'A'], ['GRCh37_bench', '6', '157425083', 'C', 'G'], ['GRCh37_bench', '6', '157796791', 'T', 'A'], ['GRCh37_bench', '6', '158665410', 'G', 'C'], ['GRCh37_bench', '6', '160854898', 'G', 'A'], ['GRCh37_bench', '6', '161814872', 'C', 'T'], ['GRCh37_bench', '6', '162742734', 'A', 'T'], ['GRCh37_bench', '6', '162858850', 'A', 'C'], ['GRCh37_bench', '6', '165361916', 'C', 'T'], ['GRCh37_bench', '6', '166950179', 'A', 'C'], ['GRCh37_bench', '6', '167625069', 'A', 'T'], ['GRCh37_bench', '6', '167908110', 'A', 'T'], ['GRCh37_bench', '6', '168907182', 'C', 'G'], ['GRCh37_bench', '6', '169835224', 'C', 'T']]
    chrom = '7'
    mutation_pd = {}

    samples = ['GRCh37_bench']
    mutation_dict = {}
    print(mut_types)
    mutation_dict["6144"] = pd.DataFrame([[[] for _ in samples] for _ in mut_types], index=mut_types, columns=samples)

    vcf_path = '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/tests/WGS/GRCh37/output/temp/test_01270cc6-d80d-4b2f-b4d1-066e6a46953b/SNV/'
    chrom_path = '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/chromosomes/tsb/GRCh37/'
    output_matrix = '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/tests/WGS/GRCh37/output/'
    exome = False



    mutation_dinuc_pd_all = pd.DataFrame([[[] for _ in samples] for _ in mut_types], index=mut_types, columns=samples)

    functionFlag = True
    bed = False
    tsb_ref = {
        0: ["N", "A"],
        1: ["N", "C"],
        2: ["N", "G"],
        3: ["N", "T"],
        4: ["T", "A"],
        5: ["T", "C"],
        6: ["T", "G"],
        7: ["T", "T"],
        8: ["U", "A"],
        9: ["U", "C"],
        10: ["U", "G"],
        11: ["U", "T"],
        12: ["B", "A"],
        13: ["B", "C"],
        14: ["B", "G"],
        15: ["B", "T"],
        16: ["N", "N"],
        17: ["T", "N"],
        18: ["U", "N"],
        19: ["B", "N"],
    }
    transcript_path = '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/chromosomes/transcripts/GRCh37/'
    seqInfo = False
    gs = False
    log_file = '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/tests/WGS/GRCh37/logs/SigProfilerMatrixGenerator_test_GRCh372024-09-26.out'
    print(mutation_dict["6144"].to_csv('test.csv'))

    for file in ['/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/vcf_files/GRCh37_bench/GRCh37_bench.vcf', '/home/mw/SigProfilerMatrixGenerator/SigProfilerMatrixGenerator/references/vcf_files/GRCh37_bench/GRCh37_bench2.vcf']:
        if not os.path.exists(chrom_path + chrom + ".txt"):
            continue
        with open(vcf_path + file) as f:
            lines = [line.strip().split() for line in f]
        lines = sorted(lines, key=lambda x: (x[0], int(x[2])))

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