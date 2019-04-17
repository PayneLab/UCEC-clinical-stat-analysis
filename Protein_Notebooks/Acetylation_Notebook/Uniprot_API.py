import requests, sys, json, time
import xml.etree.ElementTree as ET
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget,QTableWidgetItem,QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from xml.etree.ElementTree import Element, SubElement, Comment 
import pickle, os, numpy as np, pandas as pd

#Gather statistically significant sites
#from acetyle data
def gatherSignificantSites():
	import scipy.stats
	import CPTAC.Endometrial as en

	acetyl = en.get_acetylproteomics()
	clinical_attribute = "Proteomics_Tumor_Normal"
	acetyle_clinical_df = en.compare_clinical(acetyl, clinical_attribute)

	acetyle_clinical_df.loc[acetyle_clinical_df.Proteomics_Tumor_Normal == 'Adjacent_normal', 'Proteomics_Tumor_Normal'] = 'Non_Tumor'
	acetyle_clinical_df.loc[acetyle_clinical_df.Proteomics_Tumor_Normal == 'Enriched_normal', 'Proteomics_Tumor_Normal'] = 'Non_Tumor'
	acetyle_clinical_df.replace('Myometrium_normal', np.nan, inplace = True)
	tumor = acetyle_clinical_df.loc[acetyle_clinical_df[clinical_attribute] == "Tumor"]
	non_tumor = acetyle_clinical_df.loc[acetyle_clinical_df[clinical_attribute] == "Non_Tumor"]

	threshold = .05 / len(acetyle_clinical_df.columns)
	P_VALUE_INDEX = 1
	print("Threshold:", threshold)
	significantSites = []
	for num in range(1,len(acetyle_clinical_df.columns)):
		site = acetyle_clinical_df.columns[num]
		ttestRes = scipy.stats.ttest_ind(tumor[site], non_tumor[site])
		if (ttestRes[P_VALUE_INDEX] < threshold):
			significantSites.append(site)
	return significantSites

#This function takes the head(gene name)
#And maps it to its sites
#CREBBP = [1, 2, 3] etc.
def splitHeadAndTail(significantSites):
	headToTail = {}
	newListProteins = []
	listOfSites = []
	for dashedProtein in significantSites:
		head, sep, tail = dashedProtein.partition('-K')
	
		if head not in newListProteins:
			newListProteins.append(head)
		if head in headToTail:
			listOfSites.append(tail)
			headToTail[head] = listOfSites
		else:
			listOfSites = []
			listOfSites.append(tail)
			headToTail[head] = listOfSites
	return newListProteins, headToTail

def insertIntoModifiedResidues(geneName, nameToResidues, residuePosition):
	if geneName in nameToResidues:
		tempSiteArray = nameToResidues[geneName]
		tempSiteArray.append(residuePosition)
		nameToResidues[geneName] = tempSiteArray
	else:
		tempSiteArray = []
		tempSiteArray.append(residuePosition)
		nameToResidues[geneName] = tempSiteArray


def insertIntoSites(geneName, nameToSites, bindingPosition):
	if geneName in nameToSites:
		tempBindingSiteArray = nameToSites[geneName]
		tempBindingSiteArray.append(bindingPosition)
		nameToSites[geneName] = tempBindingSiteArray
	else:
		tempBindingSiteArray = []
		tempBindingSiteArray.append(bindingPosition)
		nameToSites[geneName] = tempBindingSiteArray
		
def insertIntoAltName(geneName, nameToAlternateName, altName):
	if geneName in nameToAlternateName:
		tempAlternateNameArray = nameToAlternateName[geneName]
		tempAlternateNameArray.append(altName)
		nameToAlternateName[geneName] = tempAlternateNameArray
	else:
		tempAlternateNameArray = [] 
		tempAlternateNameArray.append(altName)
		nameToAlternateName[geneName] = tempAlternateNameArray
		
def insertIntoFunction(geneName, nameToFunction, functionText):
	if geneName in nameToFunction:
		tempNameToFunction = nameToFunction[geneName]
		tempNameToFunction.append(functionText)
		nameToFunction[geneName] = tempNameToFunction
	else:
		tempNameToFunction = []
		tempNameToFunction.append(functionText)
		nameToFunction[geneName] = tempNameToFunction
		
def gatherModifiedResidues(geneName, nameToResidues, feature_elements, positions):
	for residue_elements in Element.iter(feature_elements):
		if len(positions) > 0:
			if residue_elements.get('position') != None and residue_elements.get('position') in positions:
				insertIntoModifiedResidues(geneName, nameToResidues, residue_elements.get('position'))
		else:
			if residue_elements.get('position') != None:
				insertIntoModifiedResidues(geneName, nameToResidues, residue_elements.get('position'))
				
def gatherBindingSites(geneName, nameToSites, positions, feature_elements):
	for binding_elements in Element.iter(feature_elements):
		if len(positions) > 0:
			if binding_elements.get('position') != None and binding_elements.get('position') in positions:
				insertIntoSites(geneName, nameToSites, binding_elements.get('position'))
				
		else:
			if binding_elements.get('position') != None:
				insertIntoSites(geneName, nameToSites, binding_elements.get('position'))
				
				

#Will set a T/F value for each column, that way the user can specify what they want in their df
def geneToDataFrame(geneName, modifiedResidueDescription="all", siteDescription="all",goDescription="all", 
fullName = True, alternateName = True, function=True, bindingSite=True, modifiedResidues=True,
positions=[]):

	nameToResidues = {}
	nameToFullName = {}
	nameToAlternateName = {}
	nameToFunction = {}
	nameToSites = {}
	
	tempSiteArray = []
	tempAlternateNameArray = []
	tempNameToFunction = []

	start_time_loop = time.time()

	try:
		#Call Uniprot API to get the XML response body
		requestURL = "https://www.uniprot.org/uniprot/?query="+geneName+" gene:"+geneName+" AND reviewed:yes AND organism:\"Homo sapiens (Human) [9606]\"&sort=score&format=xml"
		responseBody = requests.get(requestURL)

		#Bad request
		if not responseBody.ok:
		  responseBody.raise_for_status()
		  sys.exit()

		#The root contains all the elements inside the XML
		root = ET.fromstring(responseBody.content)

		for root_elements in Element.iter(root):
			if 'entry' in root_elements.tag:
				#Loop through elements inside entry tag
				for entry_elements in Element.iter(root_elements):
					if 'protein' in entry_elements.tag:
						for protein_elements in Element.iter(entry_elements):
							if 'recommendedName' in protein_elements.tag:
								for recName_element in Element.iter(protein_elements):
										nameToFullName[geneName] = recName_element.text.strip()
							elif 'alternativeName' in protein_elements.tag:
								for altName_elements in Element.iter(protein_elements):
									insertIntoAltName(geneName, nameToAlternateName, altName_elements.text.strip())
									
					if 'dbReference' in entry_elements.tag:
						for db_elements in Element.iter(entry_elements):
							if 'GO' == db_elements.get('type'):
								for go_elements in Element.iter(db_elements):
									if 'term' == go_elements.get('type'):
										#print(go_elements.get('value'))
										if goDescription == "all":
											go_elements.get('value')
										elif goDescription in go_elements.get('value'):
											print(geneName, go_elements.get('value'))
									
					#TO DO: Grab all the interacting proteins
									
							
					if 'feature' in entry_elements.tag:
						for feature_elements in Element.iter(entry_elements):
							if 'modified residue' == feature_elements.get('type'):
								if modifiedResidueDescription == "all":
									gatherModifiedResidues(geneName, nameToResidues, feature_elements, positions)
									
								elif modifiedResidueDescription in feature_elements.get('description'):
									gatherModifiedResidues(geneName, nameToResidues, feature_elements, positions)
										
							if 'binding site' == feature_elements.get('type'):
								if siteDescription == "all":
									gatherBindingSites(geneName, nameToSites, positions, feature_elements)
									
								elif siteDescription in feature_elements.get('description'):
									gatherBindingSites(geneName, nameToSites, positions, feature_elements)
							
								
					if 'comment' in entry_elements.tag:
						for comment_elements in Element.iter(entry_elements):
							if 'function' == comment_elements.get('type'):
								for function_elements in Element.iter(comment_elements):
									insertIntoFunction(geneName, nameToFunction, function_elements.text.strip())
									

	except:
		print("Error:", geneName)
		
	end_time_loop = time.time()
	print("Loop time:", (end_time_loop-start_time_loop))
	
	#Map Dictionaries to the DF
	df = pd.DataFrame(nameToFullName.items(), columns=['Gene', 'Full Name'])
	df['Alternate Name(s)'] = df['Gene'].map(nameToAlternateName)
	df['Modified Residues'] = df['Gene'].map(nameToResidues)
	df['Binding Sites'] = df['Gene'].map(nameToSites)
	
	return df


#Check: EP300, CREBBP
#significantSites = gatherSignificantSites()
#newListProteins, headToTail = splitHeadAndTail(significantSites)
geneToDataFrame("CREBBP", goDescription="p53 binding")
#for geneName in newListProteins:
#	if geneName == "CREBBP" or geneName == "EP300":
#		testDataFrame = geneToDataFrame(geneName, modifiedResidueDescription="acetyl", 
#		positions = headToTail[geneName], goDescription="p53 binding")
#		testDataFrame.to_csv(geneName+'.csv')
