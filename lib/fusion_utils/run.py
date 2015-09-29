#! /usr/bin/env python

import os, subprocess
import process
import filter

def comp_main(args):

    if os.path.getsize(args.fusion1) == 0:
        hOUT = open(args.output, 'w')
        hOUT.close()
        return
    
    if args.type1 in ["fusionfusion", "star_fusion", "genomon_fusion"] and args.type2 == "genomonSV":
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 500000, 10, args.type1)
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 10, 10, "genomonSV")
    elif args.type2 in ["fusionfusion", "star_fusion", "genomon_fusion"] and args.type1 == "genomonSV":
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


def filt_main(args):

    if args.type == "fusionfusion":
        filter.filter_fusionfusion(args.fusion, args.output, args.thres)
    elif args.type == "genomon_fusion":
        filter.filter_genomon_fusion(args.fusion, args.output, args.thres)
    elif args.type == "star_fusion":
        filter.filter_star_fusion(args.fusion, args.output, args.thres)
    else:
        raise ValueError("the input type should be fusionfusion, genomon_fusion or star_fusion")


