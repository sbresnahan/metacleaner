#!/usr/bin/env python

import sys
import getopt
import datetime
import os.path
import shutil
import subprocess
import os
import time
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO

def myfunc(argv):
    global arg_query
    arg_query = ""
    global arg_outdir
    arg_outdir = os.getcwd()
    global arg_blastdbdir
    arg_blastdbdir = ""
    global arg_badblastdb
    arg_badblastdb = ""
    global arg_goodblastdb
    arg_goodblastdb = ""
    global arg_badblastdbinput
    arg_badblastdbinput = "NONE"
    global arg_goodblastdbinput
    arg_goodblastdbinput = "NONE"
    global arg_chunks
    arg_chunks = "100"
    global arg_pident
    arg_pident = "100"
    global arg_qcovs
    arg_qcovs = "100"
    global arg_threads
    arg_threads = "1"
    global arg_sortonly
    arg_sortonly = "F"
    global arg_filterlevel
    arg_filterlevel = "genus"
    global arg_taxadb
    arg_taxadb = ""
    global arg_addfilter
    arg_addfilter = "NONE"
    arg_help = "\nmetacleaner options:\
    \n-q, --query <(required) path to query fasta file>\
    \n-o, --outdir <path to output directory>\
    \n-e, --pident <(default=100) percent identity threshold for filtering>\
    \n-v, --qcovs <(default=100) query cover threshold for filtering>\
    \n-s, --sortonly <(default=F), T/F: start at sorting step (requires previously generated blastn output files in output directory specified by -o)>\
    \n-l, --filterlevel <(default=genus) one of superkingdom, phylum, class, order, family, genus, or species>\
    \n-b, --blastdbdir <(required) path to badblastdb directory>\
    \n-x, --badblastdb <(required) badblastdb name>\
    \n-f, --badblastdbinput <path to fasta file for badblastdb (required if badblastdb does not already exist in directory specified by -b)>\
    \n-y, --goodblastdb <(required) goodblastdb name>\
    \n-g, --goodblastdbinput <path to fasta file for goodblastdb (required if goodblastdb does not already exist in directory specified by -b)>\
    \n-d, --taxadb <(required unless constructing from scratch) path to directory containing accessionTaxa.sql file for taxonomizr>\
    \n-c, --chunks <(default=100) number of chunks to split query file into (higher values may increase speed for larger query files)>\
    \n-t, --threads <(default=1) number of threads for blastn>\
    \n-w, --addfilter <Additional filter level, one of taxa levels as in -l (--filterlevel) and filter value, separated by a comma;\
    \n hits at pident and qcovs with taxa info other than in addfilter will be removed>".format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hq:b:x:y:f:g:o:c:e:v:t:s:l:d:w:", ["help", "query=", "blastdbdir=", "badblastdb=", \
        "goodblastdb=", "badblastdbinput=", "goodblastdbinput=", "outdir=", "chunks=", "pident=", \
        "qcovs=", "threads=", "sortonly=", "filterlevel==","taxadb==","addfilter=="])
    except:
        print(arg_help)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-q", "--query"):
            arg_query = arg
        elif opt in ("-b", "--blastdbdir"):
            arg_blastdbdir = arg
        elif opt in ("-x", "--badblastdb"):
            arg_badblastdb = arg
        elif opt in ("-y", "--goodblastdb"):
            arg_goodblastdb = arg
        elif opt in ("-f", "--badblastdbinput"):
            arg_badblastdbinput = arg
        elif opt in ("-g", "--goodblastdbinput"):
            arg_goodblastdbinput = arg
        elif opt in ("-o", "--outdir"):
            arg_outdir = arg
        elif opt in ("-c", "--chunks"):
            arg_chunks = str(arg)
        elif opt in ("-e", "--pident"):
            arg_pident = str(arg)
        elif opt in ("-v", "--qcovs"):
            arg_qcovs = str(arg)
        elif opt in ("-t", "--threads"):
            arg_threads = str(arg)
        elif opt in ("-s", "--sortonly"):
            arg_sortonly = arg
        elif opt in ("-l", "--filterlevel"):
            arg_filterlevel = arg
        elif opt in ("-d", "--taxadb"):
            arg_taxadb = arg
        elif opt in ("-w", "--addfilter"):
            arg_addfilter = arg

    print('----------------------------------')
    if not os.stat(arg_query).st_size == 0:
        print(datetime.datetime.now())
        print('> Runnning metacleaner with options:')
        print('query =', arg_query)
        print('blastdbdir =', arg_blastdbdir)
        print('badblastdb =', arg_badblastdb)
        print('badblastdbinput =', arg_badblastdbinput)
        print('goodblastdb =', arg_goodblastdb)
        print('goodblastdbinput =', arg_goodblastdbinput)
        print('outdir =', arg_outdir)
        print('chunks =', arg_chunks)
        print('pident =', arg_pident)
        print('qcovs =', arg_qcovs)
        print('threads =', arg_threads)
        print('sortonly =', arg_sortonly)
        print('filterlevel =', arg_filterlevel)
        print('taxadb =', arg_taxadb)
        print('addfilter =', arg_addfilter)
        print('----------------------------------')
    else:
        print("Error: query file specified by -q does not exist.")
        sys.exit(2)

def check_badblastdb():
    # badblastdb
    fname = arg_blastdbdir + "/" + arg_badblastdb + ".nhd"
    if os.path.isfile(fname):
        print(arg_badblastdb + " exists in directory specified by -b")
    else:
        if shutil.which("makeblastdb") is not None:
            badblastdb = arg_blastdbdir + "/" + arg_badblastdb
            if os.path.isfile(arg_badblastdbinput):
                print("badblastdb "+ arg_badblastdb + " does not exist in directory specified by -b. Attempting to construct:")
                subprocess.call(['makeblastdb', '-in', arg_badblastdbinput, '-out', badblastdb, '-dbtype',\
                'nucl', '-hash_index'])
            else:
                print("badblastdb " + arg_badblastdb + " does not exist in directory specified by -b, but fasta file specified by -f does not exist.")
                sys.exit(2)
        else:
            print("Error: BLAST command line tools are not installed.")
            sys.exit(2)

def check_goodblastdb():
    fname = arg_blastdbdir + "/" + arg_goodblastdb + ".nhd"
    if os.path.isfile(fname):
        print(arg_goodblastdb + " exists in directory specified by -b")
    else:
        if shutil.which("makeblastdb") is not None:
            goodblastdb = arg_blastdbdir + "/" + arg_goodblastdb
            if os.path.isfile(arg_goodblastdbinput):
                print("goodblastdb "+ arg_goodblastdb + " does not exist in directory specified by -b. Attempting to construct:")
                subprocess.call(['makeblastdb', '-in', arg_goodblastdbinput, '-out', goodblastdb, '-dbtype', \
                'nucl', '-hash_index'])
            else:
                print("goodblastdb " + arg_goodblastdb + " does not exist in directory specified by -b, but fasta file specified by -f does not exist.")
                sys.exit(2)
        else:
            print("Error: BLAST command line tools are not installed.")
            sys.exit(2)

def split_fasta(file_in,file_out):
    print("> Preprocessing query file")
    if shutil.which("pyfasta") is not None:
        print('Splitting query fasta file into ' + arg_chunks + " chunks")
        subprocess.call(['cp', file_in, file_out])
        subprocess.call(['pyfasta', 'split', '-n', arg_chunks, file_out], \
        stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        subprocess.call(['rm', file_out])
    else:
        print('Error: pyfasta must be installed. Do \'pip install pyfasta\'.')
        sys.exit(2)

def badblastn():
    extlen = len(arg_query.split('.', 1))
    ext = arg_query.split('.', 1)[extlen-1]
    fpath = tempdir + "/*." + ext
    print("> Searching for query sequences in "+arg_badblastdb)
    flist = glob.glob(fpath)
    chunk = 0
    for i in flist:
        if not os.stat(i).st_size == 0:
            chunk = chunk + 1
            chunkout = tempdir + "/blastn_chunk" + str(chunk) + ".out1"
            print("Running blastn search of query chunk " + str(chunk)+" | "+str(datetime.datetime.now()))
            blastdbpath = arg_blastdbdir+"/"+arg_badblastdb
            subprocess.call(['blastn', '-db', blastdbpath, '-query', i, '-out', chunkout, \
            '-outfmt', '6 qseqid pident qcovs sseqid', '-num_threads', arg_threads, '-max_target_seqs', "1"],stderr=subprocess.DEVNULL)
        else:
            subprocess.call(['rm', i])
    filenames = glob.glob(tempdir+"/*.out1")
    with open(arg_outdir+"/blastn_"+arg_badblastdb+"_output.tsv", 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def badsort():
    print("> Filtering accessions from query due to hits with "+arg_pident+" identity and "\
    +arg_qcovs+" coverage on sequences in "+arg_badblastdb+":")
    badseqids = []
    tophit = []
    with pd.read_csv(arg_outdir+"/blastn_"+arg_badblastdb+"_output.tsv", sep='\t', header=None, chunksize=10000) as reader:
        for chunk in reader:
            chunk.columns = ['qseqid', 'pident', 'qcovs', 'sseqid']
            qseqids = chunk.qseqid.unique()
            for i in qseqids:
                chunk_sub = chunk[chunk["qseqid"].isin([i])].reset_index(drop=True)
                if chunk_sub[(chunk_sub["pident"]>=float(arg_pident)) & \
                (chunk_sub["qcovs"]>=float(arg_qcovs))].shape[0] > 0:
                    badseqids.append(i)
                    tophit.append(chunk_sub.iloc[0]['sseqid'])
    print(">>> "+str(len(badseqids))+" accession(s) filtered")
    badseqout = pd.DataFrame({'badseqid' : badseqids, 'tophit' : tophit, 'reason' : ["in "+arg_badblastdb] * len(badseqids)})
    badseqout.to_csv(tempdir+"/"+fileout+"_badseqids.txt", sep='\t', header=True, index=False)
    return badseqids

def filter_fasta_badseqids(badseqids):
    headers = []
    with open(arg_query, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            headers.append(record.description)
    goodseqids = set([x for x in headers if x not in badseqids])
    print("> Saving temporary fasta file: of "+str(len(headers))+" query accession(s), "\
    +str(len(badseqids))+" accession(s) filtered, "+str(len(goodseqids))+" accession(s) retained")
    query_filtered = (r for r in SeqIO.parse(arg_query, "fasta") if r.id in goodseqids)
    SeqIO.write(query_filtered,tempdir+"/"+fileout+"_nobadseqids.fasta","fasta")

def goodblastn():
    extlen = len(arg_query.split('.', 1))
    ext = arg_query.split('.', 1)[extlen-1]
    fpath = tempdir2 + "/*." + ext
    print("> Searching for query sequences in "+arg_goodblastdb)
    flist = glob.glob(fpath)
    chunk = 0
    for i in flist:
        if not os.stat(i).st_size == 0:
            chunk = chunk + 1
            chunkout = tempdir2 + "/blastn_chunk" + str(chunk) + ".out2"
            print("Running blastn search of query chunk " + str(chunk)+" | "+str(datetime.datetime.now()))
            blastdbpath = arg_blastdbdir+"/"+arg_goodblastdb
            subprocess.call(['blastn', '-db', blastdbpath, '-query', i, '-out', chunkout, \
            '-outfmt', '6 qseqid pident qcovs sseqid', '-num_threads', arg_threads, '-max_target_seqs', "10"],stderr=subprocess.DEVNULL)
        else:
            subprocess.call(['rm', i])
    filenames = glob.glob(tempdir2+"/*.out2")
    with open(arg_outdir+"/blastn_"+arg_goodblastdb+"_output.tsv", 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def goodsort():
    print("> Filtering accessions from query due to hits with < "+arg_pident+" identity and < "\
    +arg_qcovs+" coverage on sequences in "+arg_goodblastdb+":")
    badseqids = []
    badtophit = []
    goodseqids = []
    goodtophit = []
    with pd.read_csv(arg_outdir+"/blastn_"+arg_goodblastdb+"_output.tsv", sep='\t', header=None, chunksize=10000) as reader:
        for chunk in reader:
            chunk.columns = ['qseqid', 'pident', 'qcovs', 'sseqid']
            qseqids = list(set(chunk.qseqid.unique()) - set(badseqids + goodseqids))
            for i in qseqids:
                chunk_sub = chunk[chunk["qseqid"].isin([i])].reset_index(drop=True)
                if (chunk_sub[(chunk_sub['pident']>=float(arg_pident)) & \
                (chunk_sub["qcovs"]>=float(arg_qcovs))].shape[0] > 0) or \
                (chunk_sub.iloc[0]['qseqid']==chunk_sub.iloc[0]['sseqid']):
                    goodseqids = goodseqids + list(chunk_sub[(chunk_sub['pident']>=float(arg_pident)) & \
                    (chunk_sub["qcovs"]>=float(arg_qcovs))]['qseqid'])
                    goodtophit = goodtophit + list(chunk_sub[(chunk_sub['pident']>=float(arg_pident)) & \
                    (chunk_sub["qcovs"]>=float(arg_qcovs))]['sseqid'])
                else:
                    badseqids.append(i)
                    badtophit.append(chunk_sub.iloc[0]['sseqid'])
    print(">>> "+str(len(badseqids))+" accession(s) flagged")
    badseqout = pd.DataFrame({'badseqid' : badseqids, 'tophit' : badtophit, 'reason' : ["not in "+arg_goodblastdb] * len(badseqids)})
    goodseqout = pd.DataFrame({'goodseqid' : goodseqids, 'tophit' : goodtophit})
    oldbadseqids = pd.read_csv(tempdir+"/"+fileout+"_badseqids.txt", sep='\t')
    newbadseqids = pd.concat([oldbadseqids,badseqout])
    newbadseqids.to_csv(tempdir+"/"+fileout+"_badseqids.txt", sep='\t', header=True, index=False)
    goodseqout.to_csv(tempdir+"/"+fileout+"_goodseqids.txt", sep='\t', header=True, index=False)

def taxafilter():
    scriptdir = os.path.realpath(os.path.dirname(__file__))
    print("> Collating taxonomy information and filtering sequences with mislabeled "+arg_filterlevel)
    print("taxafilter.R is in "+scriptdir)
    subprocess.call(['Rscript', scriptdir+'/taxafilter.R', arg_taxadb+"/accessionTaxa.sql",\
    tempdir+"/"+fileout+"_badseqids.txt", tempdir+"/"+fileout+"_goodseqids.txt",\
    arg_outdir, arg_filterlevel, fileout, arg_badblastdb, arg_goodblastdb, arg_addfilter])

def filter_fasta_goodseqids():
    headers = []
    with open(arg_query, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            headers.append(record.description)
    badseqids = list(pd.read_csv(tempdir+"/"+fileout+"_badseqids.txt", sep='\t', header=0)['badseqid'])
    goodseqids = set([x for x in headers if x not in badseqids])
    print("> Saving cleaned fasta file: of "+str(len(headers))+" query accession(s), "\
    +str(len(badseqids))+" accession(s) filtered, "+str(len(goodseqids))+" accession(s) retained")
    query_filtered = (r for r in SeqIO.parse(arg_query, "fasta") if r.id in goodseqids)
    SeqIO.write(query_filtered,arg_outdir+"/"+fileout+"_clean.fasta","fasta")
    newheaders = []
    with open(arg_outdir+"/"+fileout+"_clean.fasta", "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            newheaders.append(record.description)
    cleantax = pd.read_csv(arg_outdir+"/"+fileout+"_clean.tax", header=None)
    cleantax['qseqidcat'] = pd.Categorical(cleantax[1], categories=newheaders, ordered=True)
    cleantax.sort_values('qseqidcat', inplace=True)
    cleantax = cleantax.drop('qseqidcat', axis=1)
    cleantax.to_csv(arg_outdir+"/"+fileout+"_clean.tax", header=None, index=None)


def cleanup():
    subprocess.call(['cp', tempdir+"/"+fileout+"_badseqids.txt", arg_outdir+"/"+fileout+"_badseqids.txt"])
    if arg_sortonly=="F":
        subprocess.call(['rm', '-r', tempdir])
    else:
        subprocess.call(['rm', tempdir+"/"+fileout+"_goodseqids.txt"])
        subprocess.call(['rm', '-r', tempdir2])


if __name__ == "__main__":

    myfunc(sys.argv)

    global tempdir
    global tempdir2
    global fileout
    head, tail = os.path.split(arg_query)
    fileout = tail.split('.', 1)[0]

    if arg_sortonly=="T":
        tempdir = arg_outdir

    if arg_sortonly=="F":
        tempdir = arg_outdir + "/" + "temp_" + str(time.time()).split('.', 1)[0]
        print('Creating temporary directory at: ' + tempdir)
        os.makedirs(tempdir)
        check_badblastdb()
        check_goodblastdb()
        split_fasta(arg_query,tempdir+"/"+fileout+".fasta")
        badblastn()

    badseqids = badsort()
    filter_fasta_badseqids(badseqids)
    tempdir2 = tempdir + "/" + "temp_" + str(time.time()).split('.', 1)[0]
    os.makedirs(tempdir2)

    if arg_sortonly=="F":
        split_fasta(tempdir+"/"+fileout+"_nobadseqids.fasta",tempdir2+"/"+fileout+"_nobadseqids.fasta")
        goodblastn()
    goodsort()

    taxafilter()
    filter_fasta_goodseqids()
    cleanup()
