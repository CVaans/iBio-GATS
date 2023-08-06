import sys

from modeller import *
from modeller.automodel import *
#from modeller import soap_protein_od

env = Environ()

if len(sys.argv) == 4:
    alignment_file_name = sys.argv[1]
    pdb_file_name = sys.argv[2]
    sequence_file_name = sys.argv[3]

    print("Parameter 1:", alignment_file_name)
    print("Parameter 2:", pdb_file_name)
    print("Parameter 3:", sequence_file_name)
else:
    print("Please provide three command-line arguments.")

# alignment_file_name = str(sys.argv[1])
# pdb_file_name = str(sys.argv[2])
# sequence_file_name = str(sys.argv[3])

print("Alignment File name is {0}".format(alignment_file_name))
print("PDF File name is {0}".format(pdb_file_name))
print("Sequence File name is {0}".format(sequence_file_name))

if not alignment_file_name:
    print("Warning", "{} is not a valid alignment file".format(alignment_file_name))
if not pdb_file_name:
    print("Warning", "{} is not a valid pdb file".format(pdb_file_name))
if not sequence_file_name:
    print("Warning", "{} is not a valid sequence file".format(sequence_file_name))

a = automodel(env, alnfile=alignment_file_name,
              knowns=pdb_file_name, sequence=sequence_file_name,
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 5

a.make()
