'''Running BLAST over the Internet'''

# from Bio.Blast import NCBIWWW

# help(NCBIWWW.qblast)

# result_handle = NCBIWWW.qblast("blastn", "nt", "8332116")

# fasta_string = open("sequences/ls_orchid.fasta").read()
# result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string)

# from Bio import SeqIO
# record = SeqIO.read("sequences/ls_orchid.fasta", format="fasta")
# result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)

'''Running BLAST locally'''

from Bio.Blast.Applications import NcbiblastxCommandline
# help(NcbiblastxCommandline)

# blastx_cline = NcbiblastxCommandline(query="sequences/ls_orchid.fasta", db="nr", evalue=0.001, outfmt=5, out="ls_orchid.xml")
# print(blastx_cline)

# Parsing output

# result_handle = open("my_blast.xml")
# from Bio.Blast import NCBIXML
# blast_record = NCBIXML.read(result_handle)

# For multiple results

# from Bio.Blast import NCBIXML
# blast_record = NCBIXML.read(result_handle)
# blast_record = next(blast_records)
#for blast_record in blast_records:
    # Do something with blast_record

'''The BLAST record class'''

# E_VALUE_THRESH = 0.04

# for alignment in blast_record.alignments:
#     for hsp in alignment.hsps:
#         if hsp.expect < E_VALUE_THRESH:
#             print("****Alignment****")
#             print("sequence:", alignment.title)
#             print("length:", alignment.length)
#             print("e value:", hsp.expect)
#             print(hsp.query[0:75] + "...")
#             print(hsp.match[0:75] + "...")
#             print(hsp.sbjct[0:75] + "...")




