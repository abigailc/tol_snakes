#!/usr/bin/python

# created abigailc@Leviathan on december 3 2017 10pm

#concatenates aligned files by taxonid
#format of seq ids needs to be >blah|blah|TAXID|accessionnumorwhatever
# taxid must be seperated from other things by BAR and must be second to last.
# change this around line 58 if need be.

#input: list of objects each containing a variable called aligned_file pointing to the alignment
#requires Snakes_classes class Fasta_Object

#output: a concatenated fasta filename
#        a generated concat file
import Snakes_classes

def correlate_ids(list_of_lists_of_taxids):
    output_list = [] #init final output
    used = [] #init
    for item in list_of_lists_of_taxids: #for each list
        for taxid in item:  #go through each taxid and see if you've already correlated it
            if taxid in used:
                pass
            else: #if not, find it's index in each sublist
                used.append(taxid) #so we don't re-correlate it
                id_index_list = []
                id_index_list.append(taxid)
                # check the index of that id in each list and append to "id_index_list"
                # if the id is not in that list, should append "NA"
                for eachlist in list_of_lists_of_taxids:
                    try:
                        #for concat tree -- only take the FIRST HIT of any given taxid per gene
                        index = eachlist.index(taxid)
                    except ValueError:
                        index = "NA"
                    id_index_list.append(index)
                # add the result of scanning that id to overall output list.
                output_list.append(id_index_list)
    # output list looks like:
    # outputlist = [ ["CatTID",1,2,3,4] , ["DogTID", 2,1,13,14] ]
    print("Concat_Correlation finished")
    return output_list


#assuming you aligned just fine >.>
#modified from abigailc@Actaeon july '17 MakeSpeciesTree.py
def Concatenate(gene_objects, new_cc_fasta_name):
    list_of_lists_of_taxids = [] #will contain a list of the order taxids appear in each gene to concat
    list_of_aligned_fasta_objects = []
    for obj in gene_objects:
        obj.list_of_taxids = [] #will be list in same order as seqs of taxids.
        ali_fas = Snakes_classes.Fasta_Object(obj.aligned_file)
        for record in ali_fas.sequences:
            recordseq = record.seq.strip()
            recordseq = recordseq.replace("\n", "") #because you will mess that up at some point
            seqlen = len(recordseq)
            ali_fas.number_of_sites = seqlen
            bits = record.id.split("|")
            record.taxid = bits[-2]
            obj.list_of_taxids.append(record.taxid)
        list_of_aligned_fasta_objects.append(ali_fas)
        list_of_lists_of_taxids.append(obj.list_of_taxids) #add this gene list to overall list
    # go though first list, save name (index list1 ,index list2, indexlist3), and add name to "used" list
    # then go through second list, save same thing but skip any that are already in "used" list.
    indexed_ids_list = correlate_ids(list_of_lists_of_taxids)
    # indexed_ids_list in form [ ["Cat", 1,2,3],["dog",2,1,2] ]
    # this part will create a new, concatenated .fasta file
    lenlist_final = []
    with open(new_cc_fasta_name, "w") as new:
        for item in indexed_ids_list:
            lenlist = []
            new.write(">" + item[0].strip()+"\n")
            fas_num = 0
            allseq = ""
            #print(len(list_of_aligned_fasta_objects))
            #print("should be 7")
            #ensure that for each fas in fasta class list, all sequences are of the same length.
            for fas in list_of_aligned_fasta_objects:
                fas_num += 1
                # fas_num keeps track of what number fasta we are on, which
                # correlates to the index of the index in indexed_ids_list
                search_index = item[fas_num]
                # search_index will be something like "22"
                # if search_index is NA, generate a str of "-" that is n
                # characters long, where n is the return from
                # fas.number_of_sites
                if search_index == "NA":
                    ndash = fas.number_of_sites
                    retreived_seq = ""
                    for i in range(int(ndash)):
                        retreived_seq = retreived_seq + ("-")
                else:
                    #print(search_index)
                    #print(len(fas.sequences))
                    retreived_sequence = fas.sequences[search_index]
                    retreived_seq = retreived_sequence.seq.strip()
                    retreived_seq = retreived_seq.replace("\n", "")
                    # retreived_seq wil be something like "the 22nd sequence in
                    # object Fas's sequence list... " or "BLAHSEQUENCEDATA"
                count = 0
                #print(len(retreived_seq))
                allseq = allseq + retreived_seq
            #print(len(allseq))
            lenlist_final.append(len(allseq))
            newseq = ""
            for letter in allseq:
                if count > 79: #wait is fasta supposed to be 60 or 80 character wraps
                    count = 0
                    newseq = newseq + ("\n")
                newseq = newseq + letter
                count += 1
                
            new.write(newseq.strip()+"\n") #write the new sequence with correct line breaks

        for length in lenlist_final:
            if length == lenlist_final[0]:
                pass
            else:
                print("ERROR your concat sequences are not all of the same length something has gone horribly wrong! aborting.")
                raise SystemExit
        print("Should be finished generating new concatenated fasta at: " + new_cc_fasta_name)
    print("done w cc gen!")
    return new_cc_fasta_name