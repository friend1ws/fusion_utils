#! /usr/bin/env python

import os, subprocess
import process
import pysam
import filter
import annotation

def comp_main(args):

    if os.path.getsize(args.fusion1) == 0:
        hOUT = open(args.output, 'w')
        hOUT.close()
        return
    
    if args.type1 in ["fusionfusion", "star_fusion", "genomon_fusion", "mapsplice2", "tophat_fusion"] and args.type2 == "genomonSV":
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 500000, 10, args.type1)
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 10, 10, "genomonSV")
    elif args.type2 in ["fusionfusion", "star_fusion", "genomon_fusion", "mapsplice2", "tophat_fusion"] and args.type1 == "genomonSV":
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 10, 10, "genomonSV")
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 500000, 10, args.type2)
    else:
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 10, 10, args.type1)
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 10, 10, args.type2)


    hOUT = open(args.output + ".fusion_comp.bedpe", 'w')
    subprocess.call([args.bedtools_path + "/bedtools", "pairtopair", "-a", args.output + ".fusion1.bedpe", "-b", args.output + ".fusion2.bedpe"], stdout = hOUT)
    hOUT.close()

    # create dictionary
    fusion_comp = {}
    hIN = open(args.output + ".fusion_comp.bedpe", 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        fusion_comp[F[6]] = F[16]
    
    hIN.close()

    # add SV annotation to fusion
    hIN = open(args.fusion1, 'r')
    hOUT = open(args.output, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        chr1, pos1, dir1, chr2, pos2, dir2 = process.get_position(F, args.type1)
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        SV_info = fusion_comp[ID] if ID in fusion_comp else "---"
        print >> hOUT, '\t'.join(F) + '\t' + SV_info

    hIN.close()
    hOUT.close()

    # remove intermediate files
    subprocess.call(["rm", "-rf", args.output + ".fusion1.bedpe"])
    subprocess.call(["rm", "-rf", args.output + ".fusion2.bedpe"])
    subprocess.call(["rm", "-rf", args.output + ".fusion_comp.bedpe"])


def rmdup_main(args):

    if os.path.getsize(args.fusion) == 0:
        hOUT = open(args.output, 'w')
        hOUT.close()
        return

    process.convert_to_bedpe(args.fusion, args.output + ".fusion.bedpe", 10, 10, args.type)

    hOUT = open(args.output + ".fusion_comp.bedpe", 'w')
    subprocess.call([args.bedtools_path + "/bedtools", "pairtopair", "-a", args.output + ".fusion.bedpe", "-b", args.output + ".fusion.bedpe"], stdout = hOUT)
    hOUT.close()

    # create dictionary
    fusion_comp = {}
    hIN = open(args.output + ".fusion_comp.bedpe", 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if F[6] != F[16]:
            if int(F[7]) < int(F[17]) or (int(F[7]) == int(F[17]) and F[6] not in fusion_comp):
                fusion_comp[F[6]] = F[16]

    hIN.close()

    hIN = open(args.fusion, 'r')
    hOUT = open(args.output, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        chr1, pos1, dir1, chr2, pos2, dir2 = process.get_position(F, args.type)
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        if ID not in fusion_comp:
            print >> hOUT, '\t'.join(F)


    hIN.close()
    hOUT.close()

    # remove intermediate files
    subprocess.call(["rm", "-rf", args.output + ".fusion.bedpe"])
    subprocess.call(["rm", "-rf", args.output + ".fusion_comp.bedpe"])



def filt_main(args):

    annotation_dir = args.db_dir
    ref_gene_bed = annotation_dir + "/refGene.bed.gz"
    ref_exon_bed = annotation_dir + "/refExon.bed.gz"
    ens_gene_bed = annotation_dir + "/ensGene.bed.gz"
    ens_exon_bed = annotation_dir + "/ensExon.bed.gz"
    grch2ucsc_file = annotation_dir + "/grch2ucsc.txt"

    ref_gene_tb = pysam.TabixFile(ref_gene_bed)
    ref_exon_tb = pysam.TabixFile(ref_exon_bed)
    ens_gene_tb = pysam.TabixFile(ens_gene_bed)
    ens_exon_tb = pysam.TabixFile(ens_exon_bed)

    # relationship between CRCh and UCSC chromosome names
    grch2ucsc = {}
    with open(grch2ucsc_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            grch2ucsc[F[0]] = F[1]

    hIN = open(args.fusion, 'r')
    hOUT = open(args.output, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        chr1, pos1, dir1, chr2, pos2, dir2 = process.get_position(F, args.type)
        read_num = process.get_read_num(F, args.type)

        if read_num < args.thres: continue
 
        if chr1 in ["MT", "chrM", "hs37d5"] or chr2 in ["MT", "chrM", "hs37d5"]: continue
        if chr1.startswith("GL0") or chr2.startswith("GL0"): continue

        chr1_ucsc = grch2ucsc[chr1] if chr1 in grch2ucsc else chr1 
        chr2_ucsc = grch2ucsc[chr2] if chr2 in grch2ucsc else chr2
        gene1 = annotation.get_gene_info(chr1_ucsc, pos1, ref_gene_tb, ens_gene_tb)
        gene2 = annotation.get_gene_info(chr2_ucsc, pos2, ref_gene_tb, ens_gene_tb)
        junction1 = annotation.get_junc_info(chr1_ucsc, pos1, ref_exon_tb, ens_exon_tb, 5)
        junction2 = annotation.get_junc_info(chr2_ucsc, pos2, ref_exon_tb, ens_exon_tb, 5)

        junc_test = 0
        if (';'.join(junction1).find("start") != -1 and ';'.join(junction2).find("end") != -1) or (';'.join(junction1).find("end") != -1 and ';'.join(junction2).find("start") != -1):
            junc_test = 1
        
        if junc_test == 0: continue

        print >> hOUT, '\t'.join(F)

    hIN.close()
    hOUT.close()

