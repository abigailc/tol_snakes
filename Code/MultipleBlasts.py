# #!/usr/bin/python

# created abigailc@Leviathan on december 3 2017 3pm

# run from terminal or import and call run_multiple_blasts(input_fasta, projectname, directory, entrez, number

#Preform Multiple Blasts from NCBI 
#given a list of queries, with the gene name before the first | 
#given a restriction for the entrez-query, and or max hits
#save to file all hits
#convert all hits to a useable format

#support for adding taxonomic information (? maybe seperate module)
#support for being called from something else (subsampler?)


#imports
#used to split initial queries only
#from Bio import SeqIO #why are these butts such a pain i hate seqio
#used to make directories
import os
import Snakes_classes
 


#functions
def get_query_sequences(input_fasta, projectname):
	#this will open your input file, and write one query sequence per file in the folder 
	# /projectname_queries. it will also load each filename to a list to be passed to the blast-preforming fxn
	os.system("mkdir "+projectname+"_queries") #create dir for q files
	list_of_query_genes = []	#initialize list
	query_fasta_ob = Snakes_classes.Fasta_Object(input_fasta, "queries")
	for record in query_fasta_ob.sequences: #parse the input one seq at a time
		gene_name = str(record.id).split("|")[0]	#get the name of the gene
		obj = Snakes_classes.Individual_Gene(gene_name)
		obj.query_sequence = record.seq
		obj.query_file = projectname+"_queries/"+gene_name+"_query.fasta"
		query_fasta_ob.write(projectname+"_queries/"+gene_name+"_query.fasta", [record])	#write the sequence to file
		list_of_query_genes.append(obj)	#store this Gene_Query object
		
	return list_of_query_genes	#ret list of query gene objects

def run_blast_on_each_query(list_of_query_files, projectname, entrez, number):
	os.system("mkdir "+projectname+"_blasts") #create dir for b files
	for obj in list_of_query_files: #for each query
		if is_non_zero_file(obj.plain_fasta) is True:
			print ("blast output for gene "+obj.name+" appears to exist... skipping!")
			continue
		queryfile = obj.query_file
		querygene = obj.name
		outputfile = projectname+"_blasts/"+querygene+"_blast.txt"	#specify an output file 
		obj.blast_file = outputfile #save the blast file location
		os.system("blastp -remote -query " + queryfile + " -db nr -out " + outputfile +" -max_target_seqs "+number+" -entrez_query " +
 entrez + " -evalue 1e-4 -outfmt \"6 sallacc staxids sallids sseq\"")	#run the blast
		print("blast complete for gene:" +querygene)
	print("all blasts are complete")

def format_blast_outputs_as_fasta(list_of_query_genes, projectname):
	os.system("mkdir "+projectname+"_fastas")
	for obj in list_of_query_genes:
		if is_non_zero_file(obj.plain_fasta) is True:
			print ("blast output for gene "+obj.name+" appears to exist... skipping!")
			continue
		newfile = projectname+"_fastas/"+obj.name+".fasta"
		obj.plain_fasta = newfile
		with open (obj.blast_file) as old:
			with open (newfile, "w") as new:
				for line in old:
					if ";" in line:
						#for simplicity we are going to keep the initial accession number, since you 
						#should be easily able to click view identical proteins to
						#see the relevent other entry you might want. 
						#    one entry each for seperate taxids, though.
						line = line.strip()
						bits = line.split("\t")
						if ";" in bits[0]:
							accnums = bits[0].split(";")
							accnum = accnums[0]+"MSH"
						else:
							accnum = bits[0]
						taxids = bits[1].split(";")
						for n in range(len(taxids)):
							newline = ">"+taxids[n]+"|"+accnum+"\n"+bits[2]+"\n"
							new.write(newline)
						#print("PANIC. YOU NEED TO HANDLE MULTISPECIES HITS")
					else:
						line = line.strip()
						bits = line.split("\t")
						newline = ">"+bits[1]+"|"+bits[0]+"\n"+bits[2]+"\n"
						new.write(newline)
		print("Rewrote: "+obj.name)
	print("finished .fasta conversion")


#save time -- don't redo the blast OR the conversion if it's already done!
#currently all or nothing.
def check_if_blast_fasta_done_already(list_of_query_genes, projectname):
	answer = True
	for obj in list_of_query_genes:
		newfile = projectname+"_fastas/"+obj.name+".fasta"
		if is_non_zero_file(newfile) is True:
			obj.plain_fasta = newfile
		else:
			answer = False
	return answer

def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def run_multiple_blasts(input_fasta, projectname, directory, entrez, number):
	os.system("chdir "+directory)
	list_of_query_genes = get_query_sequences(input_fasta, projectname)
	if check_if_blast_fasta_done_already(list_of_query_genes, projectname) is False:
		run_blast_on_each_query(list_of_query_genes, projectname, entrez, number)
		format_blast_outputs_as_fasta(list_of_query_genes, projectname)
	return list_of_query_genes    

#parser for testing
if __name__ == "__main__":
    print("Running in terminal")
    parser = argparse.ArgumentParser(description="All")

    #optional directory
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in eg MakeSpeciesTree")
    parser.add_argument("-p", "--projectname", action = "store", help="type projectname eg Snakes_Species")
    parser.add_argument("-f", "--fasta", action = "store", help="type name of fasta containing queries eg snakes.fasta")
    parser.add_argument("-e", "--entrez", action = "store", default = False, help="type a string to restrict search to eg Serpentes")
    parser.add_argument("-n", "--number", action = "store", default = 100, help="type number of max hits")
   
    args = parser.parse_args()


    print("here we go... preforming blast on each query sequence in "+args.fasta)
    if args.projectname is False:
        print("...specificity? in the projectname please?")
        raise SystemExit
    run_multiple_blasts(args.fasta, args.projectname, args.directory, args.entrez. args.number)