from pyopenms import *

targets = list()
decoys = list()
FASTAFile().load("C:/Users/Electronica/Downloads/uniprot-yourlist_M202112026320BA52A5CE8FCD097CB85A53697A35338CFBY.fasta", targets) 
decoy_generator = DecoyGenerator()
for entry in targets:
    rev_entry = FASTAEntry(entry) 
    rev_entry.identifier = "DECOY_" + rev_entry.identifier 
    aas = AASequence().fromString(rev_entry.sequence) 
    rev_entry.sequence = decoy_generator.reverseProtein(aas).toString() 
    decoys.append(rev_entry)

target_decoy_database = "search_td.fasta"
FASTAFile().store(target_decoy_database, targets + decoys) 


