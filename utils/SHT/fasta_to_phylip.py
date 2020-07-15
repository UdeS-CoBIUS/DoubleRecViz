#!/usr/bin/env python

"""
Convert fasta alignemnts to relaxed phylip ones in constant memory.
Written by Lucas Sinclair.
Kopimi.

You can use this script from the shell like this::
$ fasta_to_phylip seqs.fasta seqs.phylip
"""
# The libraries we need #
import sys, os, random, string, re
import argparse

def build_arg_parser(path):
    parser = argparse.ArgumentParser(description="Label exons")   
    parser.add_argument('-p', '--clusters_path',  required=True)
    return parser



###############################################################################
class Sequence(object):
    """The Sequence object has a string *header* and
    various representations."""

    def __init__(self, header, seq):
        self.header = re.findall('^>(\S+)', header)[0]
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    @property
    def phylip(self):
        return self.header + " " + self.seq.replace('.','-') + "\n"

    @property
    def fasta(self):
        return ">" + self.header + "\n" + self.seq + "\n"

def fasta_parse(path):
    """Reads the file at *path* and yields
       Sequence objects in a lazy fashion"""
    header = ''
    seq = ''
    with open(path) as f:
        for line in f:
            line = line.strip('\n')
            if line.startswith('>'):
                if header: yield Sequence(header, seq)
                header = line
                seq = ''
                continue
            seq += line
    yield Sequence(header, seq)


def main_fasta_to_phylip_one(fa_path, ph_path):
    # Check that the path is valid #
    if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
    # Use our two functions #
    seqs = fasta_parse(fa_path)
    # Write the output to temporary file #
    tm_path = ph_path + '.' + ''.join(random.choice(string.ascii_letters) for i in range(10))
    # Count the sequences #
    count = 0
    with open(tm_path, 'w') as f:
        for seq in seqs:
            f.write(seq.phylip)
            count += 1
    # Add number of entries and length at the top #
    with open(tm_path, 'r') as old, open(ph_path, 'w') as new:
        new.write(" " + str(count) + " " + str(len(seq)) + "\n")
        new.writelines(old)
    # Clean up #
    command = "phyml -i "  + ph_path
    os.system(command)
        
    os.remove(tm_path)


def main_fasta_to_phylip(clusters_path):
    for root, dirs, files in os.walk(clusters_path):  
        for filename in  files:
            if filename.endswith(".fasta"):
                newname = filename.split(".")[0]
                main_fasta_to_phylip_one(clusters_path  +"/" +filename, clusters_path + "/phylipFormat/" +newname+".phylip")                            

def main_fasta_to_phylip_aux(clusters_path):
    path =  sys.path[0]      
    save_path  = path + "/clusters/sequencesClusters/"
    return main_fasta_to_phylip(clusters_path)


if __name__ == "__main__":    
    path =  sys.path[0]
    parser = build_arg_parser(path)
    arg = parser.parse_args()        
    clusters_path  = arg.clusters_path 
    main_fasta_to_phylip(clusters_path)
