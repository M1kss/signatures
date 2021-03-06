import os
import gzip
import string

GTRD_slice_path="/home/abramov/PARAMETERS/Master-lines.tsv"
alignments_path="/home/abramov/Alignments/"
dbsnp_path = "/home/abramov/ParamsForSign/00-All.sorted.vcf"
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
    a = line.strip("\n").split("\t")
    return a[0], int(a[1]), a[2], a[3:]  # chr, pos, id


def annotate_vcf(opened_vcf, name, out_path):
    annotated = 0
    inv_chr_set = set()
    db_ended = False
    with open(dbsnp_path, "r") as snps, open(out_path, 'w') as out:
        db_line = snps.readline()
        while db_line.startswith("#"):
            db_line = snps.readline()
        if not db_line:
            return annotated
        db_chr, db_pos, db_id, db_args = unpack_line(db_line)
        for vcf_line in opened_vcf:
            if not db_ended:
                try:
                    vcf_chr, vcf_pos, vcf_id, vcf_args = unpack_line(vcf_line)
                except (ValueError, IndexError):
                    if not vcf_line.startswith('#'):
                        print('Wrong line: {} {}'.format(name, vcf_line))
                    out.write(vcf_line)
                    continue
                if vcf_chr not in sorted_chromosomes:
                    if vcf_chr not in inv_chr_set:
                        print('Invalid chromosome: {}')
                        inv_chr_set.add(vcf_chr)
                    continue
                if not len(vcf_args[0]) == 1 or not len(vcf_args[1]) == 1:
                    continue
                if vcf_args[0] not in Nucleotides or vcf_args[1] not in Nucleotides:
                    continue
                while vcf_chr != 'chr' + db_chr or vcf_pos > db_pos:
                    db_line = snps.readline()
                    while db_line.startswith("#"):
                        db_line = snps.readline()
                    if not db_line:
                        db_ended = True
                        out.write(vcf_line)
                        continue
                    db_chr, db_pos, db_id, db_args = unpack_line(db_line)
                print("Now doing {}, pos {}, name {}".format(vcf_chr, vcf_pos, name))
                if vcf_pos == db_pos and db_args[1] == vcf_args[1] and db_args[0] == vcf_args[0]:
                    if vcf_id != "." and vcf_id != db_id:
                        print('Mismatch! {} {} {} {} {}'.format(name, vcf_chr, vcf_pos, vcf_id, db_id))
                    vcf_id = db_id
                    annotated += 1
                out.write(pack([vcf_chr, vcf_pos, vcf_id] + vcf_args))
            else:
                out.write(vcf_line)
    return annotated


def sort_vcf(vcf_buffer):
    sorted_lines = []
    lines_to_remember = []
    for line in vcf_buffer:
        if line.startswith("#"):
            lines_to_remember.append(line)
        else:
            if not line.startswith("chr"):
                print("WTF {}".format(line))
            sorted_lines.append(line)
    return lines_to_remember + sorted(sorted_lines, key=lambda x: x[:x.find("\t")])


def read_vcfs():
    annotated = 0
    counted_controls = set()
    counter = 0
    with open(GTRD_slice_path, "r") as master_list:
        for i, line in enumerate(master_list):
            if line[0] == "#":
                continue
            split_line = line.strip('\n').split("\t")
            vcf_path = create_path_from_gtrd_function(split_line, for_what="vcf")
            if os.path.isfile(vcf_path):
                with gzip.open(vcf_path, "rt") as vcf_buffer:
                    name = os.path.join(out_folder, remove_punctuation(split_line[4]))
                    if not os.path.isdir(name):
                        os.mkdir(name)
                    if i > 965:
                        annotated += annotate_vcf(sort_vcf(vcf_buffer),
                                                  vcf_buffer.name, os.path.join(name, split_line[6] + '.vcf'))
                        counter += 1
                        if counter % 10 == 0:
                            print('Made {} vcfs, annotated: {}'.format(counter, annotated))
            if len(split_line) > 10:
                vcf_path = create_path_from_gtrd_function(split_line, for_what="vcf", ctrl=True)
                if vcf_path in counted_controls:
                    continue
                if os.path.isfile(vcf_path):
                    with gzip.open(vcf_path, "rt") as vcf_buffer:
                        name = os.path.join(out_folder, remove_punctuation(split_line[12]))
                        if not os.path.isdir(name):
                            os.mkdir(name)
                        if i > 965:
                            annotated += annotate_vcf(sort_vcf(vcf_buffer),
                                                      vcf_buffer.name, os.path.join(name, split_line[14] + '.vcf'))
                            counter += 1
                            if counter % 10 == 0:
                                print('Made {} vcfs, annotated: {}'.format(counter, annotated))
                counted_controls.add(vcf_path)
    return annotated


if __name__ == '__main__':
    print('Total: {}'.format(read_vcfs()))
