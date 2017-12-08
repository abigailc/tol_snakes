# #!/usr/bin/python

# created abigailc@Leviathan on december 3 2017 5pm






 #eval `ssh-agent -s`
 #ssh-add /mnt/c/Users/acaro/GradSchool/TreeOfLife/Snakes/id_rsa
 #ssh-add -l
 #ssh -l abigailc -i /mnt/c/Users/acaro/GradSchool/TreeOfLife/Snakes/id_rsa eofe5.mit.edu


##########imports####################
import os
import MultipleBlasts
import TaxonomyAndSubsampling
import Snakes_align
import Snakes_classes
import Snakes_concat
import Snakes_raxml
from xml.dom import minidom
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
import time

###########################??
########1 p f except 1 p SF in colubroidea
input_fasta = "FourGenes.fasta"
projectname = "4G_Genus"
directory = "./"
entrez = "Serpentes"
number = "10000"
ranklist = "superfamily family subfamily genus ScientificName"
number_to_keep = 1
rank = "genus"
AA_MODEL = "PROTGAMMALGF"
outgroup_str = "Iguana_iguana"
outgroup_tid = "8517"

# ###########################??
# ########1 p f except 1 p SF in colubroidea
# input_fasta = "FourGenes.fasta"
# projectname = "4G_S"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family subfamily genus ScientificName"
# number_to_keep = 1
# rank = "subfamily"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"

# ############################??
# ########1 p f except 1 p SF in colubroidea
# input_fasta = "FourGenes.fasta"
# projectname = "4G_F"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family subfamily genus ScientificName"
# number_to_keep = 1
# rank = "family"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"


# # ###############set variables############BLUE
# #initial run#
# input_fasta = "snakes_test_two.fasta"
# projectname = "snakes_two"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family genus ScientificName"
# number_to_keep = 2
# rank = "genus"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"

###########################GREEN
# #less sequences
# input_fasta = "snakes_test_two.fasta"
# projectname = "snakes_three"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family genus ScientificName"
# number_to_keep = 1
# rank = "genus"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"

# ############################RED
# #less sequences
# input_fasta = "snakes_test_two.fasta"
# projectname = "SnakesFamily"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family genus ScientificName"
# number_to_keep = 1
# rank = "family"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"


# ############################cyan
# #less sequences
# input_fasta = "SnakeVenoms.fasta"
# projectname = "SnakeVenoms"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family genus ScientificName"
# number_to_keep = 1
# rank = "genus"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"

# ############################??
# ########1 p f except 1 p SF in colubroidea
# input_fasta = "snakes_test_two.fasta"
# projectname = "SnakesSubF"
# directory = "./"
# entrez = "Serpentes"
# number = "10000"
# ranklist = "superfamily family subfamily genus ScientificName"
# number_to_keep = 1
# rank = "subfamily"
# AA_MODEL = "PROTGAMMALGF"
# outgroup_str = "Iguana_iguana"
# outgroup_tid = "8517"

####needed functions
def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0


############main function#####################
#run blasts
gene_objects = MultipleBlasts.run_multiple_blasts(
	input_fasta, projectname, directory, entrez, number)

#TODO
#add_an_outgroup(taxid, gene_objects)
outgroup_objects = MultipleBlasts.run_multiple_blasts(input_fasta, projectname+"OG", directory, outgroup_str, "3")
#auto add it?
for obj in outgroup_objects:
	if is_non_zero_file(obj.plain_fasta) is False:
		obj.outgroup_seq = "None"
		continue
	obj.fas = Snakes_classes.Fasta_Object(obj.plain_fasta)
	name = obj.name
	for record in obj.fas.sequences:
		tid = record.id.split("|")[0]
		if tid == outgroup_tid:
			for item in gene_objects:
				if item.name == name:
					record.taxid = tid
					a = outgroup_str
					for i in range(len(ranklist.split())-1):
						a = "OG|"+a
					record.id = a+"|"+record.id
					item.outgroup_seq = record
#get taxonomic information
TaxonomyAndSubsampling.add_taxonomic_information_to_gene_objects(
	gene_objects, projectname, ranklist)
#subsample by rank and majority
subsampled_folder = TaxonomyAndSubsampling.preform_majority_subsampling_at_rank(number_to_keep, rank, gene_objects, projectname)

#user input if you wanna put a specific outgroup in?
for obj in gene_objects:
	if obj.outgroup_seq == "None":
		continue
	with open (obj.current, "a") as addto:
		addto.write(">"+obj.outgroup_seq.id+"\n"+obj.outgroup_seq.seq+"\n")
print("should have added outgroup sequences to the subsampled fastas")		
#actually add the outgroup seqs
input("please hit any key to continue. this is your last chance to add in any outgroups you might want - make sure you add to all gene fastas in "+subsampled_folder)

#replace seqs with full instead of aligned bit
#check if this has already been done

subsampled_folder = TaxonomyAndSubsampling.remove_duplicates_same_species(gene_objects, projectname)
print("removed dups if they're the same species -- this is toggleable")
subsampled_folder = TaxonomyAndSubsampling.ReplaceWithFullProtein(gene_objects, projectname)
print("replaced aligned blast with full protein sequence -- this is toggleable")
##################################if you tab this out also change like 125########
#align the subsampled data

#fetch names, since muscle on cluster skips that part
for obj in gene_objects:
	exist_already = True
	obj.aligned_file = projectname+"_aligned/"+obj.name+"_Muscle.fasta"
	if is_non_zero_file(obj.aligned_file) is False:
		exist_already = False
		#if they exist, skip the aligning.
print("checked to see if already aligned")
if exist_already is False:
	#align your files
	print("aligning")
	Input_List = os.listdir(subsampled_folder)
	os.chdir(subsampled_folder)
	aligned_list = Snakes_align.muscle_array_on_cluster(Input_List, projectname)
	os.system("mkdir ../"+projectname+"_aligned")
	os.system("mv *_Muscle* ../"+projectname+"_aligned/")
	os.chdir("..")
print("concatenating")
#concatenate your data correctly
new_cc_fasta_name = projectname+"_aligned/Concatenated.fasta"
cc_out = Snakes_concat.Concatenate(gene_objects, new_cc_fasta_name)
print("adding taxonomy info to the concatenated file")
#add the concat aligned file to list of genes
concat_organizer = Snakes_classes.Individual_Gene("concat")
concat_organizer.plain_fasta=cc_out
a = TaxonomyAndSubsampling.add_taxonomic_information_to_concat_object([concat_organizer], projectname, ranklist)
concat_organizer.aligned_file=a
gene_objects.append(concat_organizer)

#re-fetch the list of stuff to run
Input_List = os.listdir(projectname+"_aligned")
os.chdir(projectname+"_aligned")

#make raxml trees
print("running raxml")
Snakes_raxml.raxml_array_on_cluster(Input_List, projectname, AA_MODEL)
os.system("mkdir ../"+projectname+"_raxml_trees")
os.system("mv *RAxML_* ../"+projectname+"_raxml_trees/")
os.chdir("..")

print("reached the end")
#check if everything moved
#concatenate the subsampled data
#???
#make a tree of the subsampled data
#Cluster 



#1. do multiple blasts
#2. add taxonomy info to each if not already there
#3. do subsampling by most-common-taxa per rank
#4. extract and align

