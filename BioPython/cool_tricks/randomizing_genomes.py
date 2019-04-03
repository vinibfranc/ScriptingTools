import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def make_shuffle_record(record, new_id):
    nuc_list = list(record.seq)
    random.shuffle(nuc_list)
    return SeqRecord(Seq("".join(nuc_list), record.seq.alphabet),
                     id=new_id, description="Based on %s" % original_rec.id)

original_rec = SeqIO.read("sequences/NC_005816.gb","genbank")
shuffled_recs = (make_shuffle_record(original_rec, "Shuffled%i" % (i+1))
                 for i in range(30))
SeqIO.write(shuffled_recs, "outputs/shuffled.fasta", "fasta")