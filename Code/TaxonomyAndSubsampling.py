# #!/usr/bin/python

# created abigailc@Leviathan on december 3 2017 4pm

#using bits of append_taxonomy abigailc@actaeon october 2016

#input: list of objects with fasta files stored as obj.plain_fasta
#       list of taxonomic ranks to append
#       optional subsampling criteria
#output: list of objects with additional obj.taxo_fasta
#added replace with full protein on 12/4/17 5pm


import os
import Snakes_classes
from xml.dom import minidom
import time
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
#input: fasta file
#output: list of the taxonomy for each 
def add_taxonomic_information_to_gene_objects(list_of_gene_objects, projectname, ranklist):
    os.system("mkdir "+projectname+"_taxon") #make dir for new stuff
    os.system("mkdir "+projectname+"_subset_taxon") #make dir for new stuff
    ranks = ranklist.split() #get a list of ranks to append
    done_this_pass = {}
    overall_taxo_file = "all_taxa.txt"
    allranks = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus ScientificName".split()
    load_donethispass_from_overall(overall_taxo_file, done_this_pass)

    for obj in list_of_gene_objects: #one gene at a time
        all_taxo_file = projectname+"_taxon/"+obj.name+"_taxo.txt" #output file name
        sub_taxo_fasta_name = projectname+"_subset_taxon/"+obj.name+"_subset_taxo.fasta" #output file name
        obj.sub_taxo_fasta = sub_taxo_fasta_name #save it
        obj.all_taxo_file = all_taxo_file #save it

        if is_non_zero_file(all_taxo_file) is True:
            obj.fasta_ob = add_taxonomic_information_to_fasta_from_store(obj.plain_fasta, all_taxo_file, obj.name, done_this_pass)
        else:
            print("adding taxonomy to "+obj.name)
            obj.fasta_ob = add_taxonomic_information_to_fasta(obj.plain_fasta, obj.name, done_this_pass)
            write_info(obj.fasta_ob, all_taxo_file, allranks)
            print("should have written "+all_taxo_file)
        
        write_info_project(obj.fasta_ob, overall_taxo_file, allranks, done_this_pass)
        write_specific_ranks(obj.fasta_ob, ranks, sub_taxo_fasta_name) #write only the ranks you specified >class|subfamily|or whatever.
        print("should have written "+sub_taxo_fasta_name)
    print("finished attaching taxonomic dictionaries and writing output files")

#for use only with the species tree
def add_taxonomic_information_to_concat_object(gene_objects, projectname, ranklist = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus"):
    ranks = ranklist.split() #get a list of ranks to append
    done_this_pass = {}
    allranks = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus".split()
    for obj in gene_objects: #one gene at a time
        all_taxo_file = projectname+"_taxon/"+obj.name+"_taxo.txt" #output file name
        sub_taxo_fasta_name = projectname+"_aligned/"+obj.name+"_taxo.fasta" #output file name
        obj.sub_taxo_fasta = sub_taxo_fasta_name #save it
        obj.all_taxo_file = all_taxo_file #save it
        if is_non_zero_file(all_taxo_file) is True:
            obj.fasta_ob = add_taxonomic_information_to_fasta_from_store(obj.plain_fasta, all_taxo_file, obj.name, done_this_pass)
        else:
            print("adding taxonomy to "+obj.name)
            obj.fasta_ob = add_taxonomic_information_to_fasta(obj.plain_fasta, obj.name, done_this_pass)
            write_info(obj.fasta_ob, all_taxo_file, allranks)
            print("should have written "+all_taxo_file)
        write_specific_ranks(obj.fasta_ob, ranks, sub_taxo_fasta_name) #write only the ranks you specified >class|subfamily|or whatever.
        print("should have written "+sub_taxo_fasta_name)
        return sub_taxo_fasta_name
    

def write_specific_ranks(fasta_ob, ranks, output_name):
    to_write = [] #initial list
    #i = 0
    for record in fasta_ob.sequences:
        newid = record.id
        taxo = ""
        for rank in ranks: #go one by one through requested ranks
            #print(rank)
            # print(record.taxonomy[rank])
            try:
                a = record.taxonomy[rank]
            except KeyError:
                print(record.taxonomy)
                print("had a key error for key "+rank)
                print(record.id)
                print(fasta_ob.name)
                record.taxonomy[rank] = "NA"
            taxo = taxo+"|"+record.taxonomy[rank] #add them to the seqid
        taxo = taxo[1:]
        record.id = taxo+"|"+newid
        #if i < 15:
            #i += 1
            #print (record.id) #so i know im doing it write #REMOVE
        #to_write.append(record)
        #raise SystemExit
    fasta_ob.write(output_name)
    print("wrote taxon annotated sequences to "+output_name)

def add_taxonomic_information_to_fasta(input_fasta, name, done_this_pass):
    fasta_ob = Snakes_classes.Fasta_Object(input_fasta, name)
    for record in fasta_ob.sequences: 
        taxid = record.id.split("|")[0] #get the taxid maybe
        #print (taxid) #check
        if taxid in done_this_pass:
            taxdict = done_this_pass[taxid]
            #print(taxid)
            # print(taxdict)
            #print("line 111")
        else:
            taxdict = get_ranks_from_taxid(taxid, done_this_pass)
        record.taxid = taxid
        record.taxonomy = taxdict
    return fasta_ob

def write_info(fasta_ob, all_taxo_file, allranks):
    with open (all_taxo_file, "w") as new:
        for record in fasta_ob.sequences:
            info = ""
            for r in allranks:
                info = record.taxonomy[r]+" "+info
            info = info.strip()+"\n"
            new.write(record.taxid+" "+info)

def write_info_project(fasta_ob, all_taxo_file, allranks, done_this_pass):
    with open (all_taxo_file, "a") as new:
        for record in fasta_ob.sequences:
            if record.taxid in done_this_pass:
                continue
            info = ""
            for r in allranks:
                info = record.taxonomy[r]+" "+info
            info = info.strip()+"\n"
            new.write(record.taxid+" "+info)

def add_taxonomic_information_to_fasta_from_store(input_fasta, all_taxo_file, name, done_this_pass):
    list_of_records = []
    fasta_ob = Snakes_classes.Fasta_Object(input_fasta, name)
    for record in fasta_ob.sequences: 
        taxid = record.id.split("|")[0] #get the taxid maybe
        #print (taxid) #check
        if taxid in done_this_pass:
            taxdict = done_this_pass[taxid]
        else:
            taxdict = get_ranks_from_taxid_from_store(taxid, all_taxo_file ,done_this_pass)
        record.taxid = taxid
        record.taxonomy = taxdict
    return fasta_ob


def get_ranks_from_taxid_from_store(taxid, all_taxo_file, done_this_pass):
    ranklist = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus ScientificName"
    ranklist = ranklist.split()
    ranklist = ranklist[::-1] #reverses it
    taxdict = {} #initialize
    with open(all_taxo_file) as store:
        for line in store:
            bits = line.split()
            if taxid == bits[0]:
                iteration = 0
                for r in ranklist:
                    iteration += 1
                    try:
                        taxdict[r]=bits[iteration]
                    except:
                        print(line)
                        print(iteration)
                        print(r)
                        print(bits[iteration])
                        raise SystemExit
    done_this_pass[taxid] = taxdict
    return taxdict

def load_donethispass_from_overall(all_taxo_file, done_this_pass):
    #this is just ALWAYS loading the last entrant WTF

    if is_non_zero_file(all_taxo_file) is False:
        return ""
    ranklist = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus ScientificName"
    ranklist = ranklist.split()
    ranklist = ranklist[::-1] #reverses it
    
    with open(all_taxo_file) as store:
        for line in store:
            taxdict = {} #initialize
            skip = False
            bits = line.split()
            iteration = 0
            taxid = bits[0]
            for r in ranklist:
                iteration += 1
                try:
                    taxdict[r]=bits[iteration]
                except:
                    print("error in loading from overall taxonomy file")
                    skip = True
                    return
            if skip is False:
                done_this_pass[taxid] = taxdict
    return done_this_pass
          
  
    
#copied from internet
#should be true if file exists and is pop
def is_non_zero_file(fpath):  
    return os.path.isfile(fpath) and os.path.getsize(fpath) > 0

def get_ranks_from_taxid(taxid, done_this_pass):
    ranklist = "superkingdom kingdom phylum class order infraorder superfamily family subfamily genus"
    ranklist = ranklist.split()
    taxdict = {} #initialize
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id='+taxid # define XML location w. current taxid
    try:
        dom = minidom.parse(urlopen(url)) # parse the data -- try bc my wifi is BAD today
        #no seriously im at like 50% packet loss
        #can you believe these conditions
        #its basically torture
    except:
        print("error... waiting 10s in case it's your dumb wifi")
        try:
            time.sleep(10)
            dom = minidom.parse(urlopen(url)) 
        except:
            print("error with url:"+ url+" check your internet connection or protocols")
            print(taxid+" is getting all NA ranks")
            for r in ranklist:
                taxdict[r] = "NA"
                return taxdict
    staffs = dom.getElementsByTagName("Taxon")  #make the xml readable

    first = True
    for staff in staffs:    #just xml parsing stuff
        tid = staff.getElementsByTagName("TaxId")[0]
        tname = staff.getElementsByTagName("ScientificName")[0]
        trank = staff.getElementsByTagName("Rank")[0]
        taxid = tid.firstChild.data
        taxname = tname.firstChild.data
        taxrank = trank.firstChild.data
        if first is True:   #this is the scientific name
            first = False
            taxname = taxname.replace(".", "_") #cleanup weirdshit
            taxname = taxname.replace("-", "_")
            taxname = taxname.replace("=", "_")
            taxname = taxname.replace(",", "_")
            taxname = taxname.replace(" ", "_")
            taxname = taxname.replace("__", "_")
            taxname = taxname.replace("__", "_")
            taxdict["ScientificName"] = taxname #save it
        if taxrank in ranklist:     #cleanup
            taxname = taxname.replace(".", "")
            taxname = taxname.replace(" ", "")
            taxdict[taxrank] = taxname #save it
    for r in ranklist:  #fill in any blanks with NA
        if r in taxdict:
            pass
        else:
            taxdict[r] = "NA"
    done_this_pass[taxid] = taxdict #add it to ongoing taxonomy dictionary for SPEED reasons
    #taxduct is something like {"species":"rerio", "genus":"danio", "class":"fuckifiknow", "ScientificName":"Danio_rerio"}
    return taxdict












#############subsampling##############












#main call
def preform_majority_subsampling_at_rank(number_to_keep, rank, gene_objects, projectname):
    #pick N sequences to keep per unique string in rank. keep the most common ones possible.
    sorted_dict = optimize_choices(rank, gene_objects) #order taxid to keep the ones most common ACROSS gene datasets
    keep_taxids = [] #initialize
    os.system("mkdir "+projectname+"_subsampled") #create dir for q files
    newfolder = projectname+"_subsampled"
    for UNIQUESTR in sorted_dict: #pick top n taxids to keep per rank
        for n in range(number_to_keep):
            try:
                keep_taxids.append(sorted_dict[UNIQUESTR][n])
            except IndexError:
                pass
    for obj in gene_objects: #open each fasta 
        newseqs = []
        for sequence in obj.fasta_ob.sequences: #parse the sequences, keep those with correct taxid
            if sequence.taxid in keep_taxids:
                newseqs.append(sequence)
        outputfile = projectname+"_subsampled/"+obj.name+".fasta"
        obj.current = outputfile
        obj.fasta_ob.write(outputfile, newseqs) #write a fasta containing the correct taxid sequences only
    return newfolder


#input: records including taxonomic dictionary and what rank to extract
#optimize for the taxid(s) within each rank shared most widely
#for each unique NAME, sort the list of TAXIDs based on how many files each is present in
def optimize_choices(rank, gene_objects):
    tablist, overdict = tabulate(rank, gene_objects)    #parse the gene datasets at given rank
    sorted_dict = {}    #initial dict
    for UNIQUESTR in overdict: #for each unique string...
        inner_score = []    #initial list for storing tuples in form TAXID:appearence_value
        for taxid in overdict[UNIQUESTR]: #count in how many genes this taxonID appears
            appearence_value = 0    #set val to 0
            for tab_dict in tablist:    #go through the gene tabulations one by one
                if UNIQUESTR in tab_dict:   #because maybe gene 5 has no "cats"
                    if taxid in tab_dict[UNIQUESTR]:    #see if cat 143 exists
                        appearence_value += 1   #if so bump it up
            inner_score.append((taxid,appearence_value)) #add a tuple form: taxid, times_seen
        appearence_sorted = sorted(inner_score,key=lambda x: x[1], reverse=True) #reverse sort the tuples by second character
        list_sorted = []
        for item in appearence_sorted:
            list_sorted.append(item[0])
        sorted_dict[UNIQUESTR] = list_sorted
    return sorted_dict

def tabulate(rank, gene_objects):
    tablist = [] #this will be a list containing dictionaries in the same order as list_of_records. 
                 #each dict will contain key:RANKSTRING value:[list of taxids]
                 #i can then. hmm. make sure each taxid is only listed once per dict -- in case i forget to fix that later.
    overdict = {}   #this will contain key:RANKSTRING val:taxid(only one copy)
    for obj in gene_objects: #gene by gene
        tabdict_inner = {}
        for record in obj.fasta_ob.sequences: #sequence by sequence
            current_str = record.taxonomy[rank] #get the current string value at RANK
            if current_str in tabdict_inner:
                tabdict_inner[current_str].append(record.taxid) #add it to a tabulation of this gene alone
            else:
                tabdict_inner[current_str] = [record.taxid]
            if current_str in overdict: #add it to a tabulation accross all genes
                if record.taxid not in overdict[current_str]:
                    overdict[current_str].append(record.taxid)
            else:
                overdict[current_str] = [record.taxid]
        tablist.append(tabdict_inner)
    return tablist, overdict



def ReplaceWithFullProtein(gene_objects, projectname):
    subs_rep_folder = projectname+"_rep"
    os.system("mkdir "+subs_rep_folder) #create dir for q files
    for obj in gene_objects:
        outputfile = subs_rep_folder+"/"+obj.name+".fasta"
        if is_non_zero_file(outputfile):
            continue
        ss_fas = Snakes_classes.Fasta_Object(obj.current, "ss_"+obj.name)
        obj.current = outputfile
        for record in ss_fas.sequences:
            recordid = record.id.strip()
            bits = record.id.split("|")
            record.taxid = bits[-2]
            record.accnum = bits[-1]
            if record.accnum[-3:] == "MSH":
                record.accnum = record.accnum[:-3]
            newseq = record.seq #in case we error, keep something at least
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="+record.accnum+"&retmode=xml"
            try:
                dom = minidom.parse(urlopen(url)) # parse the data
            except:
                print("error with url:"+ url)
                print("trying once more, than will skip this sequence.... ")
                try:
                    dom = minidom.parse(urlopen(url))
                except:
                    print("didn't fix, skipping this replace")
                    continue
            staffs = dom.getElementsByTagName("GBSeq_sequence")
            for staff in staffs:
                s=staff.firstChild.data
                s = s.upper()
                #print("should be a sequence: "+s)
               # sequence = staff.getElementsByTagName("GBSeq_sequence")[0]
               # print("seq"+sequence)
            newseq=s
            record.seq = newseq
        
        ss_fas.write(outputfile)
    print("finished sequence replace with full")
    return subs_rep_folder

def remove_duplicates_same_species(gene_objects, projectname):
    nodup_folder = projectname+"_nodup"
    os.system("mkdir "+nodup_folder) #create dir for q files
    for obj in gene_objects:
        outputfile = nodup_folder+"/"+obj.name+".fasta"
        nd_fas = Snakes_classes.Fasta_Object(obj.current, "nd_"+obj.name)
        stored = []
        towrite = []
        for record in nd_fas.sequences:
            recordid = record.id.strip()
            bits = record.id.split("|")
            record.taxid = bits[-2]
            if record.taxid not in stored:
                stored.append(record.taxid) 
                towrite.append(record)
        obj.current = outputfile
        nd_fas.write(outputfile,towrite)
    print("finished keeping one per taxid")
    return nodup_folder