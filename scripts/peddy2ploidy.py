#!/usr/bin/env python3
import os,sys

def get_gender(peddyfile):
    with open(peddyfile) as infile:
        l1 = infile.readline().strip().split(',')
        l2 = infile.readline().strip().split(',')
    data = {i:j for (i,j) in zip(l1,l2)}
    s = data['sample_id']
    p = data['ped_sex']
    g = data['predicted_sex']
    if p != 'unknown':
        print("Sample %s's predicted gender is: %s"%(s, g))
        print('Using predicted gender.')
        return s, p
    else:
        print("Sample %s's predicted gender is: %s"%(s, g))
        print("Sample %s's listed gender is: %s"%(s, p))
        if g != p:
            print('#'*100 + '\n[[[WARNING]]] Predicted gender != listed gender\n'+'#'*100)
        if p != 'unknown':
            print('Using listed gender.')
            return s, g
        else:
            print('Using predicted gender.')
            return s, p

def main():
    try:
        peddyfile = sys.argv[1]
        output = sys.argv[2]
    except:
        sys.exit(sys.argv[0] + ' [Peddy sex_check.csv file] [Output ploidy vcf file]')
        
    sample, gender = get_gender(peddyfile)
    
    with open(output, 'w') as savefile:
        savefile.write('''##fileformat=VCFv4.1
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n'''%(sample))
        if gender == 'male':
            savefile.write('''chrX\t0\t.\tN\t<CNV>\t.\tPASS\tEND=10000\tCN\t1
chrX\t2781479\t.\tN\t<CNV>\t.\tPASS\tEND=155701382\tCN\t1
chrX\t156030895\t.\tN\t<CNV>\t.\tPASS\tEND=156040895\tCN\t1
chrY\t0\t.\tN\t<CNV>\t.\tPASS\tEND=57227415\tCN\t1\n''')
        else:
            savefile.write('''chrY\t0\t.\tN\t<CNV>\t.\tPASS\tEND=57227415\tCN\t0\n''')

if __name__ == '__main__':
    main()