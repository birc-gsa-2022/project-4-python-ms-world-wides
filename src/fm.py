import argparse
import sys
import re

def suffixArray(t):
    """Given T return suffix array SA(T). 
       We use Python's sorted function here 
       for simplycity, but we can do better.

    Args:
        t (str): input string.
    """  
    satups = sorted([(t[i:], i) for i in range(len(t))])
    
    return list(map(lambda x: x[1], satups))

def cigar_to_dict(cigar: str):
    lst = [(int(i), op) for i, op in re.findall(r"(\d+)([^\d]+)", cigar)]
    if(cigar!=''):
        dict = {lst[i][1]: lst[i][0] for i in range(0, len(lst))}
    return dict

def dict_to_cigar(dict):
    cigar = []
    for key in dict:
        cigar.append(dict[key])
        cigar.append(key)
    return ''.join(str(e) for e in cigar)

def count_to_bucket(count: str):
    '''
    >>> cigar_to_edits("1$4i1m2p4s")
    '0$1i5m6p8s'
    '''
    C = cigar_to_dict(count)
    letter = list(C.keys())
    bucket = {}
    bucket[letter[0]] = 0
    for i in range(1,len(letter)):
        bucket[letter[i]] = bucket[letter[i-1]] + C[letter[i-1]]

    return dict_to_cigar(bucket)
    

#BWT
def bwt_C(x):
    x += '$'
    sa = suffixArray(x)
    col = len(x)-1

    bwt = ''.join([x[(i + col)%len(x)] for i in sa])
    count = ''.join([x[i] for i in sa])
    C = count_to_bucket(rle(count))         #  C is the bucket list
    #O = calc_O(sa, C)
    return bwt#, C, O

def calc_O(sa, C):
    
    for k in C.keys():
        C[k] = [0]
    
    O = C
    for s in sa:
        for k in O.keys():
            if k == s:
                O[k].append(O[k][-1]+1)
            else:
                O[k].append(O[k][-1])
    
    return O        

def fasta_func(fastafile):
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

def rle(s):
    rle_bwt = ""
    i = 0
    while (i <= len(s)-1):
        cnt = 1
        chr = s[i]
        j = i
        while (j < len(s)-1):
            if (s[j] == s[j+1]):
                cnt += 1
                j += 1
            else:
                break
        rle_bwt += str(cnt)+chr
        i = j+1
    return rle_bwt

def process_file(fasta_dict, filename):
    filename = filename.split('.')[0]
    file = filename + '_rle.txt'

    with open(file, 'w') as f:
        for k, v in fasta_dict.items():
            final = '>' + k + '/n'
            bwt, C, O = bwt_C(v)
            comp_bwt = rle(bwt)
            final += '#' + comp_bwt + '/n'
            comp_C = C
            final += '+' + comp_C + '/n'
            # add O
        f.writelines(final)


def main():

    # x = 'aaaaaaa'
    # f = 'genome'

    # save_as_file(x, f)

    # argparser = argparse.ArgumentParser(
    #     description="FM-index exact pattern matching",
    #     usage="\n\tfm -p genome\n\tfm genome reads"
    # )
    # argparser.add_argument(
    #     "-p", action="store_true",
    #     help="preprocess the genome."
    # )
    # argparser.add_argument(
    #     "genome",
    #     help="Simple-FASTA file containing the genome.",
    #     type=argparse.FileType('r')
    # )
    # argparser.add_argument(
    #     "reads", nargs="?",
    #     help="Simple-FASTQ file containing the reads.",
    #     type=argparse.FileType('r')
    # )
    # args = argparser.parse_args()

    # if args.p:
    #     print(args.genome.name)
    #     print(f"Preprocess {args.genome}")
    #     fasta_dict = fasta_func(args.genome)
    #     # for k, v in fasta_dict.items():
    #     #     save_as_file(v, k)
    #     process_file(fasta_dict, args.genome.name)
    # else:
    #     # here we need the optional argument reads
    #     if args.reads is None:
    #         argparser.print_help()
    #         sys.exit(1)
    #     print(f"Search {args.genome} for {args.reads}")

    print(count_to_bucket('1$4i1m2p4s'))


if __name__ == '__main__':
    main()
