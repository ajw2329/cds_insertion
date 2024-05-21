# Count kmers and calculate GC-content over the transcriptome.
import argparse
import sys
from Bio import SeqIO

parser = argparse.ArgumentParser()
	parser.add_argument("-feat", type = str, help = "Feature CSV file with required columns: [tx_len, cds_len, utr5_len, utr3_len].")
	parser.add_argument("-fa", type = str, help = "Transcriptome FASTA file.")
	parser.add_argument("-o", type = str, help = "Output file name.")
	parser.add_argument("-kmin", type = int, help = "Min kmer length.")
	parser.add_argument("-kmax", type = int, help = "Max kmer length.")
	
args = parser.parse_args(args)

# Produce kmers.
text = 'ACGT'
kmers = sum((list(product(text,repeat=i)) for i in range(args.kmin, args.kmax + 1)),[])
kmers = [''.join(kmer) for kmer in kmers]

# Load transcript information.
df = pd.read_csv(args.feat)[['V1', 'tx_len', 'cds_len', 'utr5_len', 'utr3_len']]
# Make into dict: {transcript: [tx_length, CDS_length, UTR5_length, UTR3_length], ...}
dfDict = df.set_index(['V1']).T.to_dict('list')

# Load transcript sequences.
txIDs1 = []
seqs1 = []

fasta_sequences = SeqIO.parse(open(args.fa),'fasta')
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    txIDs1.append(name)
    seqs1.append(sequence)
    
# Filter out transcripts from the FASTA that aren't in the final dataset.
txIDs = []
seqs = []
for i in range(len(txIDs1)):
    if txIDs1[i] in dfDict.keys():
        txIDs.append(txIDs1[i])
        seqs.append(seqs1[i])

print('initial transcript list length: ', len(txIDs1), 
      '\nfiltered transcript list length: ', len(txIDs))

# Establish empty frequency dictionaries.
seqDict = {new_list: [] for new_list in txIDs}
seqDictGC = {new_list: [] for new_list in txIDs}
cdsDictGC = {new_list: [] for new_list in txIDs}
utr5DictGC = {new_list: [] for new_list in txIDs}
utr3DictGC = {new_list: [] for new_list in txIDs}
sequenceDict = {new_list: [] for new_list in txIDs}

# Count kmers for each transcript.
for i in range(len(seqs)):
    txID = txIDs[i]
    seq = seqs[i]
    L = len(seq)
    cds = seq[dfDict[txID][2] + 1: dfDict[txID][0] - dfDict[txID][3]]
    utr5 = seq[:dfDict[txID][2] + 1]
    utr3 = seq[dfDict[txID][1] + dfDict[txID][2]:]
    seqDictGC[txID].append((seq.count('C') + seq.count('G')) / L)
    sequenceDict[txID].append(seq)
    if len(cds) > 0:
        cdsDictGC[txID].append((cds.count('C') + cds.count('G')) / len(cds))
        for kmer in kmers:
            l = len(kmer)
            if (len(cds) - l + 1) > 0:
                seqDict[txID].append(cds.count(kmer) / (len(cds) - l + 1))
    else:
        cdsDictGC[txID].append(0)
        seqDict[txID].append(0)
    if len(utr5) > 0:
        utr5DictGC[txID].append((utr5.count('C') + utr5.count('G')) / len(utr5))
    else:
        utr5DictGC[txID].append(0)
    if len(utr3) > 0:
        utr3DictGC[txID].append((utr3.count('C') + utr3.count('G')) / len(utr3))
    else:
        utr3DictGC[txID].append(0)

# Collate data into a single large dataframe.
merged = pd.DataFrame.from_dict(seqDict, orient='index')
merged.columns = kmers
merged.to_csv(analysis_dir + 'transcriptome_kmerFreq.csv', index=True)

seqDfGC = pd.DataFrame.from_dict(seqDictGC, orient='index')
seqDfGC.columns = ['tx_GC']
cdsDfGC = pd.DataFrame.from_dict(cdsDictGC, orient='index')
cdsDfGC.columns = ['CDS_GC']
utr5DfGC = pd.DataFrame.from_dict(utr5DictGC, orient='index')
utr5DfGC.columns = ['UTR5_GC']
utr3DfGC = pd.DataFrame.from_dict(utr3DictGC, orient='index')
utr3DfGC.columns = ['UTR3_GC']
seqDf = pd.DataFrame.from_dict(sequenceDict, orient='index')
seqDf.columns = ['sequence']

merged = seqDfGC.join(cdsDfGC).join(utr5DfGC).join(utr3DfGC)
merged.to_csv(args.o, index=True)
