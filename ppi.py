import re
import os
from collections import defaultdict as ddict

################################################################################
#	read in PPI file from PSICQUIC datasets
################################################################################

#reads the lines of a text file into a list of the lines
def readLines(filename):
	return [l for l in open(filename,"r").read().splitlines() if len(l) > 0]

#load in a homologene conversion
def loadHomologene(filename):
	homologene_lines = readLines(filename)
	homologene_dict = dict()
	
	for line in homologene_lines:
		
		split_line = line.split("\t")
		index, species, data = split_line[0], split_line[1], split_line[2:]
		
		if index not in homologene_dict: homologene_dict[index] = dict()
		homologene_dict[index][species] = data

	return homologene_dict

#make a mouse to human conversion dictionary
def mouse2Human(filename):
	homologene_dict = loadHomologene(filename)
	result_dict = dict()

	for index in homologene_dict:
		entry = homologene_dict[index]
		if ("9606" in entry) and ("10090" in entry):
			result_dict[entry["10090"][1]] = entry["9606"][1]

	return result_dict 

#loads in a uniprot style conversion table of gene symbols
def loadUniprotConversions(filename):
	lines = readLines(filename)
	convert_dict = dict()
	for line in lines:
		split_line = line.split("\t")
		if len(split_line) == 3:
			uniprot = split_line[0]
			if uniprot not in convert_dict:
				convert_dict[uniprot] = dict()
			convert_dict[uniprot][split_line[1]] = split_line[2]
	return convert_dict

#loads in a PSICQUIC style PPI file - XMI 2.5 format (I think it's called)
def readPSICQUIC(ppi_file):
	ppi_lines = readLines(ppi_file)
	dict_list = []
	col_labels = ppi_lines[0].split("\t")
	for line in ppi_lines[1:]:
		split_line = line.split("\t")
		interaction = dict()
		for i in range(len(col_labels)):
			#check that the file is properly structured
			assert(len(col_labels) == len(split_line)) 
			interaction[col_labels[i]] = split_line[i]
		dict_list.append(interaction)
	return dict_list

#get the interactors form a MINT dataset line
def mintInteractors(interaction):
	interactorA_field = interaction["#ID(s) interactor A"]
	interactorB_field = interaction["ID(s) interactor B"]
	idA = re.findall("uniprotkb:([^\W_]{6})",interactorA_field)[0]
	idB = re.findall("uniprotkb:([^\W_]{6})",interactorB_field)[0]
	return idA, idB

#get the interactors from an innateDB dataset line
def innateDBInteractors(interaction):
	interactorA_field = interaction["Alias(es) interactor A"]
	interactorB_field = interaction["Alias(es) interactor B"]
	idA = re.findall("uniprotkb:([^\W_]{6})",interactorA_field)[0]
	idB = re.findall("uniprotkb:([^\W_]{6})",interactorB_field)[0]
	return idA, idB

#get the interactors from a BAR dataset line
def barInteractors(interaction):
	interactorA_field = interaction["Alt. ID(s) interactor A"]
	interactorB_field = interaction["Alt. ID(s) interactor B"]
	idA = re.findall("uniprotkb:([^\W_]{6})",interactorA_field)[0]
	idB = re.findall("uniprotkb:([^\W_]{6})",interactorB_field)[0]
	return idA, idB

#get the interactors from a BioGrid dataset line
def bioGridInteractors(interaction):
	interactorA_field = interaction["#ID(s) interactor A"]
	interactorB_field = interaction["ID(s) interactor B"]
	idA = re.findall("entrez gene/locuslink:(\w+)",interactorA_field)[0]
	idB = re.findall("entrez gene/locuslink:(\w+)",interactorB_field)[0]
	return idA, idB

#get the interactors from an IntAct dataset line
def intActInteractors(interaction):
	interactorA_field = interaction["#ID(s) interactor A"]
	interactorB_field = interaction["ID(s) interactor B"]
	idA = re.findall("uniprotkb:([^\W_]{6})",interactorA_field)[0]
	idB = re.findall("uniprotkb:([^\W_]{6})",interactorB_field)[0]
	if "(association)" in interaction['Interaction type(s)']:
		raise NameError("Not a physical association")
	return idA, idB

#load the Entrez gene symbol conversion table for humans
def loadHumanEntrez(filename):
	convert_dict = dict()
	entrez_file = readLines(filename)
	for line in entrez_file:
		try:
			split_line = line.split("\t")
			entrez_id = split_line[1]
			gene_symbol = split_line[0]
			convert_dict[entrez_id] = gene_symbol
		except:
			pass
	return convert_dict

#load the Entrez gene symbol conversion table for mice
def loadMouseEntrez(filename):
	convert_dict = dict()
	mgi_entrez_file = readLines(filename)
	for line in mgi_entrez_file:
		try:
			split_line = line.split("\t")
			entrez_id = split_line[8]
			gene_symbol = split_line[1]
			convert_dict[entrez_id] = gene_symbol
		except:
			pass
	return convert_dict

#load in both Entrez gene symbol conversion tables (humans and mice)
def loadEntrez(mouse_file, human_file):
	mouse_convert = loadMouseEntrez(mouse_file)
	human_convert = loadHumanEntrez(human_file)
	return mouse_convert, human_convert

#convert an Entrez Gene ID to a gene symbol using tables loaded in
def entrezGeneConvert(interactor, species, mouse_convert, 
					  human_convert, mouse_2_human):
	if species == "human":
		return human_convert[interactor]
	elif species == "mouse":
		mouse_ID = mouse_convert[interactor]
		return mouse_2_human[mouse_ID]
	else:
		raise NameError("Unsuported Species")

#get the species for both interactors - only human and mouse supported
def getSpecies(interaction):
	speciesA_field = interaction["Taxid interactor A"]
	speciesB_field = interaction["Taxid interactor B"]

	if "9606" in speciesA_field: speciesA = "human"
	elif "10090" in speciesA_field: speciesA = "mouse"
	else: raise NameError("Unsuported Species")

	if "9606" in speciesB_field: speciesB = "human"
	elif "10090" in speciesB_field: speciesB = "mouse"
	else: raise NameError("Unsuported Species")

	return speciesA,speciesB

def loadUniprot(mouse_file,human_file):
	mouse_convert = loadUniprotConversions(mouse_file)
	human_convert = loadUniprotConversions(human_file)
	return mouse_convert, human_convert

#convert from a uniprot ID to a gene symbol
def uniprotConvert(interactor, species, mouse_convert, 
				   human_convert, mouse_2_human):
	if species == "human":
		return human_convert[interactor]["Gene_Name"]
	elif species == "mouse":
		mouse_ID = mouse_convert[interactor]["Gene_Name"]
		return mouse_2_human[mouse_ID]
	else:
		raise NameError("Unsuported Species")

#additional function - find the type of an interactor
def interactionType(interaction):
	interaction_field = interaction['Interaction type(s)']
	return re.findall("\((.+)\)",interaction_field)[0]

#main processing function - uses a lot of the other functions as arguments
def processPSICQUIC(convertLoadFn, filename, mouse_file, human_file, 
	homologene_file, getInteractors, getSpecies, convertID, 
	additionalInfoFn = None):
	
	mouse_convert, human_convert = convertLoadFn(mouse_file,human_file)
	mouse_2_human = mouse2Human(homologene_file)
	data_file = readPSICQUIC(filename)
	error_list = []
	output_list = []
	print "----processing PSICQUIC file----"
	print "initial data file size : %d" % len(data_file)

	for interaction in data_file:
		
		try:
			interactorA, interactorB = getInteractors(interaction)
			speciesA, speciesB = getSpecies(interaction)
			
			idA = convertID(interactorA, speciesA, mouse_convert,
							human_convert, mouse_2_human)
			idB = convertID(interactorB, speciesB, mouse_convert,
							human_convert, mouse_2_human)
			pmid_field = interaction['Publication Identifier(s)']
			pmid = re.search("pubmed:(\w+)",pmid_field).groups(1)[0]

			if additionalInfoFn != None:
				added_info = additionalInfoFn(interaction)
				output_list.append((idA, idB, pmid, added_info))
			else:
				output_list.append((idA, idB, pmid))

		except NameError:
			pass

		except:
			error_list.append(interaction)
	print "final list size w/ duplicates: %d" % len(output_list)
	output_list = list(set(output_list))
	print "final list size w/o duplicates: %d" % len(output_list)	
	return output_list, error_list

#removes interactions found in high throughput studies
def removeHighThroughput(ppi_list, thresh):
	pmid_dict = ddict(list)
	final_list = []
	for ppi in ppi_list:
		pmid_dict[ppi[2]].append(ppi)
	for pmid in pmid_dict:
		if len(pmid_dict[pmid]) <= thresh:
			final_list.extend(pmid_dict[pmid])
	return final_list

#writes out the final PPI list to a text file
def writePPI(ppi_list, output_filename):
	output = open(output_filename,"w")
	for ppi in ppi_list:
		line_list = [ppi[0],"NA","NA","NA","NA",ppi[1],"NA","NA","NA","NA",
					 "0","Binding",ppi[2]]
		if len(ppi) == 4 and ppi[3] != "physical association":
			line_list[11] = ppi[3].replace(" ","_")
		line = " ".join(line_list)
		output.write(line + "\n")
	output.close()

################################################################################
# pre-canned functions
################################################################################

def runMINT():
	MINT = processPSICQUIC(loadUniprot, 
						   "PPI/MINT.txt",
						   "PPI/MOUSE_10090_idmapping.dat",
							"PPI/HUMAN_9606_idmapping_shortened.dat", 
							"homologene.data", 
							mintInteractors,
							getSpecies,
							uniprotConvert)
	MINT_fixed = removeHighThroughput(MINT[0],10)
	writePPI(MINT_fixed,"MINT.sig")

def runInnateDB():
	innateDB = processPSICQUIC(loadUniprot, 
							   "PPI/innateDB.txt",
							   "PPI/MOUSE_10090_idmapping.dat",
							   "PPI/HUMAN_9606_idmapping_shortened.dat", 
							   "homologene.data", 
							   innateDBInteractors,
							   getSpecies,
							   uniprotConvert)
	innateDB_fixed = removeHighThroughput(innateDB[0],10)
	writePPI(innateDB[0],"innateDB.sig")

def runDIP():
	DIP = processPSICQUIC(loadUniprot, 
						  "PPI/DIP.txt",
						  "PPI/MOUSE_10090_idmapping.dat",
					      "PPI/HUMAN_9606_idmapping_shortened.dat", 
					      "homologene.data", 
					      mintInteractors,
						  getSpecies,
						  uniprotConvert)
	DIP_fixed = removeHighThroughput(DIP[0],10)
	writePPI(DIP_fixed,"DIP.sig")

def runIntAct():
	IntAct = processPSICQUIC(loadUniprot,
							 "PPI/IntAct.txt",
							 "PPI/MOUSE_10090_idmapping.dat",
							 "PPI/HUMAN_9606_idmapping_shortened.dat", 
							 "homologene.data", 
							 intActInteractors,
						     getSpecies,uniprotConvert,
						     additionalInfoFn=interactionType)
	IntAct_fixed = removeHighThroughput(IntAct[0],10)
	writePPI(IntAct_fixed,"IntAct.sig")

def runBioGrid():
	BioGrid = processPSICQUIC(loadEntrez, 
							  "PPI/BioGrid.txt", 
							  "MGI_EntrezGene.rpt", 
						  	  "entrez-gene-id.txt", 
						  	  "homologene.data",
						  	  bioGridInteractors, 
						  	  getSpecies, 
						  	  entrezGeneConvert)
	BioGrid_fixed = removeHighThroughput(BioGrid[0],10)
	writePPI(BioGrid_fixed,"BioGrid.sig")
