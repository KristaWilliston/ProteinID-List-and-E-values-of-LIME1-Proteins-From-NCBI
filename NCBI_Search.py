import os.path
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

Entrez.email = 'kmw305@georgetown.edu'

#search for RefSeq human LIME1 protein from NCBI
#E-utilities search query and get LIME1 IDs
query = "Homo sapiens[Organism] AND LIME1[Gene Name] AND REFSEQ"

handle = Entrez.esearch(db='protein', term=query)
result = Entrez.read(handle)
handle.close()

#Extracting list of IDs
id_list=result["IdList"]
print("ID List:",id_list)
print()

#Blast protein IDs against RefSeq mouse proteins
protein_ids = '\n'.join(id_list)
for protein_id in id_list:
    filename = ("blastp-refseq-LIME1"+protein_id+".xml")
    #check if name exists
    if os.path.exists(filename):
        blast_save = open(filename, "r")
        blast_results = blast_save.read()
        blast_save.close()
    else:
        #BLAST
        result = NCBIWWW.qblast("blastp","refseq_protein",protein_id,
                                entrez_query="Mus musculus[Organism]",expect=1e-3)
        blast_results = result.read()
        result.close()
    
        #save BLAST resulting XML
        blast_save = open(filename, "w")
        blast_save.write(blast_results)
        blast_save.close()

    #read the blast result
    blast_result = NCBIXML.read(open(filename))

    #find the E-value using the alignment and description parameters 
    for description in blast_result.descriptions:
        e_value = description.e

        if e_value <= 1e-5:
            print("Protein ID:", protein_id)
            print("Description:", description.title) 
            print("E-Value:", e_value)
            print()
