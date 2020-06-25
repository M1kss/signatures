import os
import gzip
import string

GTRD_slice_path="/home/abramov/PARAMETERS/Master-lines.tsv"
alignments_path="/home/abramov/Alignments/"
dbsnp_path = "/home/abramov/ParamsForSign/00-All.vcf.gz"
sorted_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
out_folder = "/home/abramov/Signatures/CellTypes/"

Nucleotides = {'A', 'T', 'G', 'C'}


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def remove_punctuation(x):
    table = str.maketrans({key: "_" for key in string.punctuation if key not in {'-', '+'}})
    return x.translate(table).replace(" ", "_")


def create_path_from_gtrd_function(line, for_what, ctrl=False):
    end = ""
    if for_what == "vcf":
        end = ".vcf.gz"
    if for_what == "annotated_table":
        end = "_table_annotated.txt"
    if for_what == "p-value_table":
        end = "_table_p.txt"
    if ctrl:
        return alignments_path + "CTRL/" + line[10] + "/" + line[14] + end
    else:
        return alignments_path + "EXP/" + line[1] + "/" + line[0] + "/" + line[6] + end


def unpack_line(line):
    a = line.split()
    return a[0], int(a[1]), a[2], a[3:]  # chr, pos, id


def annotate_vcf(opened_vcf, out_path):
    annotated = 0
    with gzip.open(dbsnp_path, "rt") as snps, open(out_path, 'w') as out:
        db_line = snps.readline()
        while db_line.startswith("#"):
            db_line = snps.readline()
        db_chr, db_pos, db_id, db_args = unpack_line(db_line)
        for vcf_line in opened_vcf:
            try:
                vcf_chr, vcf_pos, vcf_id, vcf_args = unpack_line(vcf_line)
            except ValueError:
                print('Wrong line: {} {}'.format(opened_vcf.name, vcf_line))
                out.write(vcf_line)
                continue
            if vcf_chr not in sorted_chromosomes:
                continue
            if not len(vcf_line[3]) == 1 or not len(vcf_line[4]) == 1:
                continue
            if vcf_line[3] not in Nucleotides or vcf_line[4] not in Nucleotides:
                continue
            while vcf_chr != db_chr or vcf_pos > db_pos:
                db_line = snps.readline()
                db_chr, db_pos, db_id, db_args = unpack_line(db_line)
            if vcf_pos == db_pos:
                if vcf_id and vcf_id != db_id:
                    print('Mismatch! {} {} {} {} {}'.format(opened_vcf.name, vcf_chr, vcf_pos, vcf_id, db_id))
                vcf_id = db_id
                annotated += 1
            out.write(pack([vcf_chr, vcf_pos, vcf_id] + vcf_args))
    return annotated


def read_vcfs():
    annotated = 0
    counted_controls = set()
    with open(GTRD_slice_path, "r") as master_list:
        for line in master_list:
            if line[0] == "#":
                continue
            split_line = line.strip('\n').split("\t")
            vcf_path = create_path_from_gtrd_function(split_line, for_what="vcf")
            if os.path.isfile(vcf_path):
                with gzip.open(vcf_path, "rt") as vcf_buffer:
                    name = os.path.join(out_folder, remove_punctuation(split_line[4]))
                    if not os.path.isdir(name):
                        os.mkdir(name)
                    annotated += annotate_vcf(vcf_buffer, os.path.join(name, split_line[6] + '.vcf'))
            if len(split_line) > 10:
                vcf_path = create_path_from_gtrd_function(split_line, for_what="vcf", ctrl=True)
                if vcf_path in counted_controls:
                    continue
                if os.path.isfile(vcf_path):
                    with gzip.open(vcf_path, "rt") as vcf_buffer:
                        name = os.path.join(out_folder, remove_punctuation(split_line[12]))
                        if not os.path.isdir(name):
                            os.mkdir(name)
                        annotated += annotate_vcf(vcf_buffer, os.path.join(name, split_line[14] + '.vcf'))
                counted_controls.add(vcf_path)
    return annotated


if __name__ == '__main__':
    print('Total: {}'.format(read_vcfs()))
