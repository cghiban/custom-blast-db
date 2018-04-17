#!/usr/bin/env python3

from Bio import SeqIO
import time
import sys
import os
import sqlite3


def usage():

    me = sys.argv[0]
    print("""\nUsage:
    %s <sqlite.db> <accesstion_list_file> [<rRNA product filter>]

    %s test2.db accns.txt
    %s test2.db new-accns.txt "16S ribosomal RNA"
    %s test2.db new-accns.txt "cytochrome oxidase subunit I"
    """ % (me, me, me, me))


def search(db, accns):
    data = {}
    conn = sqlite3.connect(db)
    c = conn.cursor()
    res = c.execute("""
            SELECT file,version
            FROM entries WHERE version in (%s)
        """ % ','.join(['?' for a in accns]), accns)

    for r in res:
        if r[0] in data:
            data[r[0]].append(r[1])
        else:
            data[r[0]] = [r[1]]

    c.close()
    return data


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# --------------------------------------------------------------
# __main__
#

if len(sys.argv) < 3:
    usage()
    sys.exit(1)

db_file = sys.argv[1]

# this is the file with a list of accesstions
input_file = sys.argv[2]

# only look for certain products, eg:
#  cytochrome oxidase subunit I
#  16S ribosomal RNA
filter_by_product = None
if len(sys.argv) > 3:
    filter_by_product = sys.argv[3]

#eprint(sys.argv)

#cwd = cwd();

accns = []
#daccns = {}
with open(input_file) as f:
    for line in f:
        id = line.rstrip()
        accns.append(id)
        #daccns[id] = 1

#eprint(len(accns))
found = 0
x = 0

data = search(db_file, accns)
for gb in data:
    eprint(gb, len(data[gb]))
    #x += len(data[gb])

    if not os.path.exists(gb):
        continue

    accns = {a:1 for a in data[gb]}
    #print(accns)

#if True:
    #gb = '300-01963-new.gb'
    #gb = '300-01892-new.gb'
    #gb = 'aaa.gb'
    t0 = time.time()
    for seq in SeqIO.parse(gb, "genbank"):
        if seq.id not in accns:
            continue
        #if seq.id in daccns:
        #    found += 1
        #eprint(seq.id, len(seq.features))
        features = [
                f for f in seq.features if f.type in ['source', 'rRNA']
            ]
        info = {}
        for ft in features:
            #eprint("\t", ft.type, str(ft.location))
            if ft.type  == 'source':
                #print(ft.qualifiers)
                for q in ft.qualifiers:
                    qval = ','.join(ft.qualifiers[q])
                    if (qval != ''):
                        info[q] = qval

            # if this filter is set, we reject any seq that doesn't mach it
            elif filter_by_product is not None and ft.type == 'rRNA':
                if 'product' in ft.qualifiers and filter_by_product in ft.qualifiers['product']:
                    start = ft.location.start.position
                    end = ft.location.end.position
                    info['seq'] = str(ft.extract(seq.seq))


        # since we're not filtering by a gene/product, we'll get the entire seq
        if 'seq' not in info and filter_by_product is None:
            info['seq'] = str(seq.seq)
        #eprint(info)

        seq_len = len(info['seq']) if 'seq' in info else 0
        if (seq_len < 150 or seq_len > 9000):
            continue

        classification = ','.join(seq.annotations['taxonomy'])
        seq_id = '; '.join([
                seq.id,
                '='.join(['organism', seq.annotations['organism']]),
                '='.join(['taxon', info['db_xref'][6:]]),
                '='.join(['classification', classification]),
                seq.description
            ])
        print(">%s\n%s\n" % (seq_id, info['seq']))
        #break
    ss = time.time() - t0
    eprint("BioPython:", ss, "seconds")

eprint("found:", found)

