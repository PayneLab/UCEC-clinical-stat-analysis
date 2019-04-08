import requests, sys, json, time
import xml.etree.ElementTree as ET
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget,QTableWidgetItem,QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from xml.etree.ElementTree import Element, SubElement, Comment 
import pickle, os, numpy as np, pandas as pd
from tabulate import tabulate

######################################################################################
#Here we are gathering all the significant sites for tumor vs non tumor
######################################################################################
'''import pandas as pd
import numpy as np
import scipy.stats
import statsmodels.stats.multitest
import matplotlib.pyplot as plt
import seaborn as sns
import CPTAC.Endometrial as en
import math

start_time_sig = time.time()
acetyl = en.get_acetylproteomics()
clinical_attribute = "Proteomics_Tumor_Normal"
acetyle_clinical_df = en.compare_clinical(acetyl, clinical_attribute)

acetyle_clinical_df.loc[acetyle_clinical_df.Proteomics_Tumor_Normal == 'Adjacent_normal', 'Proteomics_Tumor_Normal'] = 'Non_Tumor'
acetyle_clinical_df.loc[acetyle_clinical_df.Proteomics_Tumor_Normal == 'Enriched_normal', 'Proteomics_Tumor_Normal'] = 'Non_Tumor'
acetyle_clinical_df.replace('Myometrium_normal', np.nan, inplace = True) # Set 'Myometrium_normal' to NaN
tumor = acetyle_clinical_df.loc[acetyle_clinical_df[clinical_attribute] == "Tumor"]
non_tumor = acetyle_clinical_df.loc[acetyle_clinical_df[clinical_attribute] == "Non_Tumor"]

threshold = .05 / len(acetyle_clinical_df.columns) #Here we find the significant p-value that we want to work with
P_VALUE_INDEX = 1
print("Threshold:", threshold)
significantSites = []
for num in range(1,len(acetyle_clinical_df.columns)):
    site = acetyle_clinical_df.columns[num]
    oneSite = acetyle_clinical_df[[clinical_attribute, site]]
    ttestRes = scipy.stats.ttest_ind(tumor[site], non_tumor[site])
    if (ttestRes[P_VALUE_INDEX] < threshold): #Check if there is a significant enough difference between data points
        significantSites.append(site)

print("Number of significant sites: ", len(significantSites))
end_time_sig = time.time()
print("total time:", (end_time_sig-start_time_sig))'''
######################################################################################
#We now should be able to pass anything.. Into a function and call uniprot's API to return a data frame
######################################################################################

def inputToDataFrame(name, residues="all", sites="all", function=True, positions=[]):

	#Goal: Have the gene name be the key and link all items to this key, such as function,
	#Amino acid modifications, binding sites
	nameToResidues = {}
	nameToFullName = {}
	nameToAlternateName = {}
	nameToFunction = {}
	nameToSites = {}
	
	tempSiteArray = []
	tempAlternateNameArray = []
	tempNameToFunction = []

	#for name in newListProteins:

	#intTemp += 1
	start_time_loop = time.time()

	try:
		#Call Uniprot API to get the XML response body
		requestURL = "https://www.uniprot.org/uniprot/?query="+name+" gene:"+name+" AND reviewed:yes AND organism:\"Homo sapiens (Human) [9606]\"&sort=score&format=xml"
		responseBody = requests.get(requestURL)

		#Bad request
		if not responseBody.ok:
		  responseBody.raise_for_status()
		  sys.exit()

		#The root contains all the elements inside the XML
		root = ET.fromstring(responseBody.content)

		#Parse through XML to gather the Accession ID
		for root_elements in Element.iter(root):
			if 'entry' in root_elements.tag:
				#Loop through elements inside entry tag
				for entry_elements in Element.iter(root_elements):
					if 'protein' in entry_elements.tag:
						for protein_elements in Element.iter(entry_elements):
							if 'recommendedName' in protein_elements.tag:
								for recName_element in Element.iter(protein_elements):
										nameToFullName[name] = recName_element.text.strip()
							elif 'alternativeName' in protein_elements.tag:
								for altName_elements in Element.iter(protein_elements):
									if name in nameToAlternateName:
										tempAlternateNameArray = nameToAlternateName[name]
										tempAlternateNameArray.append(altName_elements.text.strip())
										nameToAlternateName[name] = tempAlternateNameArray
									else:
										tempAlternateNameArray = []
										tempAlternateNameArray.append(altName_elements.text.strip())
										nameToAlternateName[name] = tempAlternateNameArray
								
					if 'feature' in entry_elements.tag:
						for feature_elements in Element.iter(entry_elements):
							if 'modified residue' == feature_elements.get('type'):
								if residues == "all":
									for residue_elements in Element.iter(feature_elements):
										if residue_elements.get('position') != None:
											print(name, residue_elements.get('position'))
											if name in nameToResidues:
												tempSiteArray = nameToResidues[name]
												tempSiteArray.append(residue_elements.get('position'))
												nameToResidues[name] = tempSiteArray
											else:
												tempSiteArray = []
												tempSiteArray.append(residue_elements.get('position'))
												nameToResidues[name] = tempSiteArray
								elif residues in feature_elements.get('description'):
									for residue_elements in Element.iter(feature_elements):
										if residue_elements.get('position') != None:
											print(name, residue_elements.get('position'))
											if name in nameToResidues:
												tempSiteArray = nameToResidues[name]
												tempSiteArray.append(residue_elements.get('position'))
												nameToResidues[name] = tempSiteArray
											else:
												tempSiteArray = []
												tempSiteArray.append(residue_elements.get('position'))
												nameToResidues[name] = tempSiteArray
								
									
					if 'comment' in entry_elements.tag:
						for comment_elements in Element.iter(entry_elements):
							if 'function' == comment_elements.get('type'):
								for function_elements in Element.iter(comment_elements):
									if name in nameToFunction:
										tempNameToFunction = nameToFunction[name]
										tempNameToFunction.append(function_elements.text.strip())
										nameToFunction[name] = tempNameToFunction
									else:
										tempNameToFunction = []
										tempNameToFunction.append(function_elements.text.strip())
										nameToFunction[name] = tempNameToFunction

	except:
		print("Error:", name)
		
	end_time_loop = time.time()
	print("Loop time:", (end_time_loop-start_time_loop))
	
	#Map Dictionaries to the DF
	df = pd.DataFrame(nameToFullName.items(), columns=['Gene', 'Full Name'])
	df['Alternate Name(s)'] = df['Gene'].map(nameToAlternateName)
	df['Modified Residues'] = df['Gene'].map(nameToResidues)
	df['Binding Sites'] = df['Gene'].map(nameToSites)
	
	return df
		

#Removed dashes from proteins (head) and K from numbers (tail)
'''headToTail = {}
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
        headToTail[head] = listOfSites'''

#Questions: Should we allow an array of gene names? And they can just pass in a single gene in 
#an array if thats all they want? modifiedResidue description, binding site description, specific sites,
#

#GeneName, fullName, alternateName, modifiedResidue Description, binding site description, specific sites, 
testDataFrame = inputToDataFrame("ACIN1", "all", )

#testDataFrame.to_html('Uniprot.html')
testDataFrame.to_csv('Uniprot.csv')

#print(tabulate(testDataFrame, headers='keys', tablefmt='psql'))
