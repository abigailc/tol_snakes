# #!/usr/bin/python

# created abigailc@Leviathan on december 3 2017 6pm


###########classes#############
#barebones fasta seqio-like object - can read in, populate sequences, be modified, and write itself.
class Fasta_Object:
    def __init__(self, filename="error", name="unspecified"):
        self.sequences = []
        self.populate(filename)
    def populate(self, filename):
        start = True
        index = 0
        self.file = filename
        with open (filename) as old:
            for line in old:
                if ">" in line:
                    if start == True:
                        start = False
                    else:
                        AAseq = AAseq.strip()
                        self.sequences.append(Single_Seq(newid,AAseq))
                        index += 1
                    AAseq = ""
                    # format the new seqID
                    newid = line.strip()
                    newid = newid.strip(">")
                else:
                    AAseq = AAseq + line
            AAseq=AAseq.strip()
            # catch the last AAseq pass
            self.sequences.append(Single_Seq(newid,AAseq))
        print("generated a fasta object containing "+str(len(self.sequences))+ "sequences")
    def write(self, outputfile, data=False):
        if data == False:
            data = self.sequences
        with open (outputfile, "w") as new:
            for item in data:
                new.write(">"+item.id+"\n"+item.seq+"\n")
        print("wrote sequences to file: "+outputfile)

#such that each sequence can have nice attached annotation and be more modular than in seqio
class Single_Seq:
    def __init__(self, ide, seq):
        self.id = ide
        self.seq = seq
        self.taxid = ""
        self.taxonomy = {}

#classes to keep track of files
class Individual_Gene:
    def __init__(self, name="whatever"):
        self.name = name
        #query_sequence
        #query_file
        #blast_file
        self.outgroup_seq = "None"
        self.plain_fasta = ""
        #fasta_ob
        pass