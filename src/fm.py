import argparse
import sys
import re
import ast

def suffixArray(x: str) -> list:
    """Given x return suffix array SA(x). 
       We use Python's sorted function here 
       for simplycity, but we can do better.

    Args:
        x (str): input string.
    """  
    satups = sorted([(x[i:], i) for i in range(len(x))])
    
    return list(map(lambda t: t[1], satups))


def count_to_bucket(count: str) -> dict:
    '''
    >>> count_to_bucket("$iiiimppss")
    bucket = {$ : 0, i : 1, m : 5, p : 6, s : 8}
    '''
    C = {}
    for i,c in enumerate(count):
        if c in C:
            continue
        else:
            C[c] = i

    return C
    
def bwt_C_O(x: str) -> tuple():
    
    x += '$'
    sa = suffixArray(x)
    col = len(x)-1 # column at the back

    bwt = ''.join([x[(i + col)%len(x)] for i in sa])
    count = ''.join([x[i] for i in sa])
    C = count_to_bucket(count) # dict with cumulative counts
    O = calc_O(bwt, C) # dict (table) with offsets
    return sa, C , O

def calc_O(bwt: str, C: dict) -> dict:
    '''
    >>>calc_O('aaba$')
    O = { '$' : [0, 0, 0, 1, 1, 1], 'a' : [0, 1, 1, 1, 2, 3], 'b' : [0, 0, 1, 1, 1, 1]}
    '''
    O = C.copy()
    for k in O.keys():
        O[k] = [0]
    
    for char in bwt:
        for k in O.keys():
            if k == char:
                O[k].append(O[k][-1]+1)
            else:
                O[k].append(O[k][-1])
    
    return O        

def fasta_func(fastafile: str) -> dict:
    '''Function that can take file or list of strings and return dictionary
    with fasta sequence coupled with its sequence name'''

    sequence = []
    name = ''
    fasta_dict = {}
    for line in fastafile:
        if type(line) == list:
            line = line[0]
        if line.startswith('>'):
            if name != '':
                fasta_dict[name] = ''.join(sequence)
                sequence = []
            name = line[1:].strip()
        else:
            sequence.append(line.strip())

    if name != '':
        fasta_dict[name] = ''.join(sequence)

    return fasta_dict

def fastq_func(fastqfile: str) -> dict: 
    read = []
    name = ''
    fastq_dict = {}
    for line in fastqfile:
        if line.startswith('@'):
            if name != '':
                fastq_dict[name] = ''.join(read)
                read = []
            name = line[1:].strip()
        else:
            read.append(line.strip())

    if name != '':
        fastq_dict[name] = ''.join(read)

    return fastq_dict
    
def process_file(fasta_dict: dict, filename: str):
    '''Create a file containing the name of each string 
    with their suffix array, their bucket dict and their O table'''
    filename = filename.split('.')[0]
    file = filename + '_prepro.txt'

    with open(file, 'w') as f:
        final = ''
        for k, v in fasta_dict.items():
            final += '>' + k + '\n'
            sa, C, O = bwt_C_O(v)
            final += '#' + str(sa) + '\n'
            final += '+' + str(C) + '\n'
            final += '@' + str(O) + '\n'
        f.writelines(final)

def fm_search(prepro_file: str, reads: str) -> str:

    # give room for muliple fasta sequences
    fastanames, sa_list, C_list, O_list = [], [], [], []

    with open(prepro_file, 'r') as f:    
        lines = f.readlines()
        for line in lines:
            if line.startswith('>'):
                name = line[1:].strip()
                fastanames.append(name)
            elif line.startswith('#'):
                l = ast.literal_eval(line[1:].strip())
                sa_list.append(l)
            elif line.startswith('+'):
                d = ast.literal_eval(line[1:].strip())
                C_list.append(d) # convert back to dict
            elif line.startswith('@'):
                d = ast.literal_eval(line[1:].strip())
                O_list.append(d) # convert back to dict
    
    fastq_dict = fastq_func(reads)

    res = []
    for readname, read in fastq_dict.items():
        
        if len(sa_list) == len(C_list) == len(O_list):
            # For each fasta sequence
            for i in range(len(sa_list)):
                genomename = fastanames[i]

                sa, O, C = sa_list[i], O_list[i], C_list[i]
                L, R = 0, len(sa)
                
                for char in reversed(read):
                    if L == R or char not in C:
                        L = R # for char not in C
                        break
                    else:
                        L = C[char] + O[char][L]
                        R = C[char] + O[char][R]
                
                for a in range(L,R):
                    match = sa[a]+1
                    res.append('\t'.join([readname, genomename, str(match), f'{str(len(read))}M', read]))
    
    return '\n'.join(res)

def main():

    argparser = argparse.ArgumentParser(
        description="FM-index exact pattern matching",
        usage="\n\tfm -p genome\n\tfm genome reads"
    )
    argparser.add_argument(
        "-p", action="store_true",
        help="preprocess the genome."
    )
    argparser.add_argument(
        "genome",
        help="Simple-FASTA file containing the genome.",
        type=argparse.FileType('r')
    )
    argparser.add_argument(
        "reads", nargs="?",
        help="Simple-FASTQ file containing the reads.",
        type=argparse.FileType('r')
    )
    args = argparser.parse_args()

    if args.p:
        print(f"Preprocess {args.genome}")
        fasta_dict = fasta_func(args.genome)

        process_file(fasta_dict, args.genome.name)
    else:
        # here we need the optional argument reads
        if args.reads is None:
            argparser.print_help()
            sys.exit(1)
        print(f"Search {args.genome} for {args.reads}")
        
        prepro_file = args.genome.name.split('.')[0]+'_prepro.txt'
        fm_search(prepro_file, args.reads)



if __name__ == '__main__':
    main()
