#! /usr/bin/env python

import subprocess
import process

def main(args):
    
    fusion1 = args.fusion1
    type1 = args.type1
    fusion2 = args.fusion2
    type2 = args.type2
    output = ar


    if args.type1 in ["fusionfusion", "star_fusion", "genomon_fusion"] and args.type2 == "genomonSV":
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 500000, 10, args.type1)
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 10, 10, "genomonSV")
    elif args.type2 in ["fusionfusion", "star_fusion", "genomon_fusion"] and args.type1 == "genomonSV":
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 10, 10, "genomonSV")
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 500000, 10, args.type2)
    else:
        process.convert_to_bedpe(args.fusion1, args.output + ".fusion1.bedpe", 10, 10, args.type1)
        process.convert_to_bedpe(args.fusion2, args.output + ".fusion2.bedpe", 10, 10, args.type2)


    hOUT = open(output_file + ".fusion_comp.bedpe", 'w')
    subprocess.call([args.bedtools_path + "/bedtools", "pairtopair", "-a", args.output + ".fusion1.bedpe", "-b", args.output + ".fusion2.bedpe"], stdout = hOUT)
    hOUT.close()

    # create dictionary
    fusion_comp = {}
    hIN = open(output_file + ".fusion_comp.bedpe", 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        fusion_comp[F[6]] = F[16]
    
    hIN.close()

    # add SV annotation to fusion
    hIN = open(args.fusion1, 'r')
    hOUT = open(args.output, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        chr1, pos1, dir1, chr2, pos2, dir2 = process.get_position(F, method)
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        SV_info = fusion_to_SV[ID] if ID in fusion_to_SV else "---"
        print >> hOUT, '\t'.join(F) + '\t' + SV_info

    hIN.close()
    hOUT.close()

    # remove intermediate files
    subprocess.call(["rm", "-rf", output_file + ".fusion1.bedpe"])
    subprocess.call(["rm", "-rf", output_file + ".fusion2.bedpe"])
    subprocess.call(["rm", "-rf", output_file + ".fusion_comp.bedpe"])

