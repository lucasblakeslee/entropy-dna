from Bio.Blast import NCBIWWW
from Bio import SeqIO

def main():
    record = SeqIO.read("test_genomes/gamma_test_genomes/amphritea-japonica.fasta", format = "fasta")
    result_handle = blast_search(record)
    blast_result = open("my_blast.xml", "w+")
    blast_result.write(result_handle.read())
    blast_result.close()
    result_handle.close()

def blast_search(SeqIOrecord):
    result_handle = NCBIWWW.qblast("blastn", "nt", SeqIOrecord.seq[1:100])
    return result_handle



if __name__ == "__main__":
    main()

