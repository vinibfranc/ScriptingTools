# The first step to using the mapper is to get the exons from a GenBank or similar file.
# The mapper will accept exons as a sequence of pairs, a SeqRecord with a CDS feature, or a CDS SeqFeature.
# The file used in this example is located in the Tests directory of the Biopython source code.

from Bio.SeqUtils.Mapper import CoordinateMapper
from Bio import SeqIO
 
def get_first_CDS(parser):
    for rec in parser:
        for feat in rec.features:
            if feat.type == "CDS" and len(feat.location.parts) > 1:
                return feat
exons = get_first_CDS(SeqIO.parse("sequences/cor6_6.gb", "genbank"))
cm = CoordinateMapper(exons)

# Once the mapper is constructed, its methods can be used to transform positions located within the given CDS.
# Note that attempting to transcribe and translate a genomic position that is not within CDS will raise an exception.
# Also note that printing the list will show the repr of the positions, which are in Python 0-based coordinates.

sample_g_values = [50, 150, 250, 350, 450, 550, 650]
 
protein_positions = []
for raw_g_pos in sample_g_values:
    # EAFP method
    try:
        # Dialect argument converts g_pos from Genbank to Python coordinates
        p_pos = cm.g2p(raw_g_pos, dialect="genbank")
    except ValueError:
        p_pos = None
    protein_positions.append(p_pos)
print(protein_positions)

# Here's an example function that prints a table of the genomic, CDS, and protein coordinates
# given a coordinate mapper and a list of genomic values.

from Bio.SeqUtils.Mapper import GenomePosition
 
def gcp_table(mapper, g_list):
    """Print a table of genomic, CDS, and protein coordinates"""
    # Print header
    print("%4s | %6s | %2s" % ("g", "CDS", "p"))
    print("-" * 20)
    for g_pos in g_list:
        # Directly convert g_pos from Genbank to Python coordinates
        g_pos = GenomePosition.from_dialect("genbank", g_pos)
        c_pos = mapper.g2c(g_pos)
        # LBYL method:
        if c_pos.pos_type == "exon":
            p_pos = mapper.c2p(c_pos)
        else:
            p_pos = ""
        # Print formatted row
        print("%4s | %6s | %2s" % (g_pos, c_pos, p_pos))
 
gcp_table(cm, sample_g_values)