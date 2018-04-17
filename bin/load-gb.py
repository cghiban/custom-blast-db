#!/usr/bin/env python3

from Bio import (SeqIO, BiopythonParserWarning)
import time
import sys
import os
import sqlite3
import random

# with a few ideas gotten from:
# http://charlesleifer.com/blog/going-fast-with-sqlite-and-python/
#

def usage():

    me = sys.argv[0]
    print("""\nUsage:
    %s <sqlite.db> <genbank_file>

    %s test2.db ./fasta-20180410/300-00001.gb
    """ % (me, me))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_db(path):
    
    conn = sqlite3.connect(path, isolation_level=None)
    #version, = conn.execute('select sqlite_version();').fetchone()
    #eprint("** Using sqlite v%s" % version)

    #zz = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='entries';").fetchone()
    #eprint(zz)

    conn.execute('pragma journal_mode=wal;')
    conn.execute("""
        CREATE TABLE IF NOT EXISTS entries (
            id rowid,
            accession text not null,
            version text not null,
            organism text not null,
            modified_on date not null,
            file text not null,
            UNIQUE (version)
        )
    """);

    return conn

#@transaction
def search(conn, accns):
    data = {}
    c = conn.cursor()
    res = c.execute(""" SELECT version FROM entries WHERE version in (%s)""" \
            % ','.join(['?' for a in accns]),
            accns
        )

    data = { r[0]: 1 for r in res }
    c.close()
    return data

# --------------------------------------------------------------
# __main__
#

if len(sys.argv) < 3:
    usage()
    sys.exit(1)

db_file = sys.argv[1]

# this is the file with a list of accesstions
input_file = os.path.relpath(sys.argv[2])

conn = get_db(db_file)

#num = 0
# new data to add into db
new_entries = []
gbiter = SeqIO.parse(input_file, "genbank")
while True:
    try:
        seq = next(gbiter)
        new_entries.append({
            'accession': seq.name,
            'version': seq.id,
            'organism': seq.annotations['organism'],
            'modified_on': seq.annotations['date'],
            'file': input_file,
        })
    except ValueError as ex:
        eprint('E: ingnoring entry.', ex)
    except StopIteration as ex:
        break

if len(new_entries) == 0:
    eprint("** No new entries in [%s] Bailing out.." % input_file)
    sys.exit(0)

#eprint('** Done parsing..')
#sys.exit(0)

# make sure these don't already exist in the db
sentries = search(conn, [e['version'] for e in new_entries])

random.seed()
time.sleep(random.random())  # sleep under 1 second

# start transaction
conn.execute('BEGIN')
try:
    c = conn.cursor()
    for ne in new_entries:
        if ne['version'] not in sentries:
            #eprint('** about to add:', ne['version'])
            c.execute("""
                INSERT INTO entries (accession, version, organism, modified_on, file) 
                VALUES (?, ?, ?, ?, ?)
            """, [ ne['accession'], ne['version'], ne['organism'], ne['modified_on'], ne['file'] ])
    c.close()
except sqlite3.IntegrityError as ex:
    conn.rollback()  # Roll back all changes if an exception occurs.
    #print(ex)
    raise
else:
    conn.commit()

num, = conn.execute('SELECT total_changes()').fetchone()
if num != 0:
    eprint('** %s\tadded %d entries' % (input_file, num))

conn.close()

