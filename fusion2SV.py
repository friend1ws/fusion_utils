#! /usr/bin/env python

import subprocess, re

ReReg = re.compile(r'([^ \t\n\r\f\v,]+):([\+\-])(\d+)\-([^ \t\n\r\f\v,]+):([\+\-])(\d+)')

def convert_to_bedpe(input_file, output_file, margin_major, margin_minor, method):

    if method not in ["fusionfusion", "genomonSV", "star_fusion", "genomon_fusion"]:
        raise ValueError("the argument method should be fusionfusion, genomonSV, star_fusion")

    hIN = open(input_file, 'r')
    hOUT = open(output_file, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')

        if method == "star_fusion" and F[0] == "#fusion_name": continue

        if method in ["fusionfusion", "genomonSV"]:
            chr1, pos1, dir1, chr2, pos2, dir2 = F[0], F[1], F[2], F[3], F[4], F[5]
        elif method == "star_fusion":
            chr1, pos1, dir1 = F[4].split(':')
            chr2, pos2, dir2 = F[7].split(':')
            dir2 = '-' if dir2 == '+' else '-'
        elif method == "genomon_fusion":
            keyMatch = ReReg.match(F[0])
            chr1, dir1, pos1, chr2, dir2, pos2 = keyMatch.group(1), keyMatch.group(2), keyMatch.group(3), keyMatch.group(4), keyMatch.group(5), keyMatch.group(6)

        if chr1 > chr2 or chr1 == chr2 and int(pos1) > int(pos2):
            chr1, chr2, pos1, pos2, dir1, dir2 = chr2, chr1, pos2, pos1, dir2, dir1

        start1, end1, start2, end2 = pos1, pos1, pos2, pos2
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        chr1 = chr1.replace("chr", "")
        chr2 = chr2.replace("chr", "")
       
        if method == "fusionfusion":
            read_num = F[7]
        elif method == "genomonSV":
            read_num = F[13]
        elif method == "star_fusion":
            read_num = F[1]
        elif method == "genomon_fusion":
            read_num = F[20]
 
        if dir1 == '+':
            start1 = str(int(start1) - int(margin_minor))
            end1 = str(int(end1) + int(margin_major))
        else:
            start1 = str(int(start1) - int(margin_major))
            end1 = str(int(end1) + int(margin_minor))

        if dir2 == '+':
            start2 = str(int(start2) - int(margin_minor))
            end2 = str(int(end2) + int(margin_major))
        else:
            start2 = str(int(start2) - int(margin_major))
            end2 = str(int(end2) + int(margin_minor))

        print >> hOUT, '\t'.join([chr1, start1, end1, chr2, start2, end2, ID, read_num, dir1, dir2])

    hIN.close()
    hOUT.close()


def fusion2SV(input_fusion, input_SV, output_file, method, bedtools_path):

    convert_to_bedpe(input_fusion, output_file + ".fusion.bedpe", 500000, 50, method)
    convert_to_bedpe(input_SV, output_file + ".SV.bedpe", 5, 5, "genomonSV")

    hOUT = open(output_file + ".fusion_SV.bedpe", 'w')
    subprocess.call([bedtools_path + "/bedtools", "pairtopair", "-a", output_file + ".fusion.bedpe", "-b", output_file + ".SV.bedpe"], stdout = hOUT)
    hOUT.close()

    # create dictionary
    fusion_to_SV = {}
    hIN = open(output_file + ".fusion_SV.bedpe", 'r')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        fusion_to_SV[F[6]] = F[16]
    
    hIN.close()

    # add SV annotation to fusion
    hIN = open(input_fusion, 'r')
    hOUT = open(output_file, 'w')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        if method == "star_fusion" and F[0] == "#fusion_name": continue
        if method in ["fusionfusion", "genomonSV"]:
            chr1, pos1, dir1, chr2, pos2, dir2 = F[0], F[1], F[2], F[3], F[4], F[5]
        elif method == "star_fusion":
            chr1, pos1, dir1 = F[4].split(':')
            chr2, pos2, dir2 = F[7].split(':')
            dir2 = '-' if dir2 == '+' else '-'
        elif method == "genomon_fusion":
            keyMatch = ReReg.match(F[0])
            chr1, dir1, pos1, chr2, dir2, pos2 = keyMatch.group(1), keyMatch.group(2), keyMatch.group(3), keyMatch.group(4), keyMatch.group(5), keyMatch.group(6)

        if chr1 > chr2 or chr1 == chr2 and int(pos1) > int(pos2):
            chr1, chr2, pos1, pos2, dir1, dir2 = chr2, chr1, pos2, pos1, dir2, dir1

        start1, end1, start2, end2 = pos1, pos1, pos2, pos2
        ID = chr1 + ':' + dir1 + pos1 + '-' + chr2 + ':' + dir2 + pos2

        SV_info = fusion_to_SV[ID] if ID in fusion_to_SV else "---"
        print >> hOUT, '\t'.join(F) + '\t' + SV_info

    hIN.close()
    hOUT.close()

    # remove intermediate files
    # subprocess.call(["rm", "-rf", output_file + ".fusion.bedpe"])
    # subprocess.call(["rm", "-rf", output_file + ".SV.bedpe"])
    # subprocess.call(["rm", "-rf", output_file + ".fusion_SV.bedpe"])



if __name__ == "__main__":
    import sys
    fusion2SV(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
 
