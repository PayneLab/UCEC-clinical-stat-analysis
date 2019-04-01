import requests, sys, json, time
import xml.etree.ElementTree as ET
from PyQt5.QtWidgets import QMainWindow, QApplication, QWidget, QAction, QTableWidget,QTableWidgetItem,QVBoxLayout
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from xml.etree.ElementTree import Element, SubElement, Comment 
import pickle, os, numpy as np, pandas as pd

######################################################################################
#Here we are gathering all the significant sites for tumor vs non tumor
######################################################################################
import pandas as pd
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
print("total time:", (end_time_sig-start_time_sig))
######################################################################################
#We now should be able to pass anything.. Into a function and call uniprot's API to return a data frame
######################################################################################

def inputToDataFrame(name):

	nameToSites = {}
	nameToFullName = {}
	nameToAlternateName = {}
	nameToFunction = {}
	
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
							if 'modified residue' == feature_elements.get('type') and 'N6-acetyllysine' in  feature_elements.get('description'):
								for residue_elements in Element.iter(feature_elements):
									if residue_elements.get('position') in headToTail[name]:
										print(name, residue_elements.get('position'))
										if name in nameToSites:
											tempSiteArray = nameToSites[name]
											tempSiteArray.append(residue_elements.get('position'))
											nameToSites[name] = tempSiteArray
										else:
											tempSiteArray = []
											tempSiteArray.append(residue_elements.get('position'))
											nameToSites[name] = tempSiteArray
										headToTail[name].remove(residue_elements.get('position'))
									
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
	df = pd.DataFrame(nameToSites.items(), columns=['Gene', 'Modified Residues'])
	#Here we can map other dictionaries to this dataframe...
	#df['exampleColumnName'] = df['key'].map(exampleDictionary)
	return df
		

#Removed dashes from proteins (head) and K from numbers (tail)
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

testDataFrame = inputToDataFrame(newListProteins[0])
print(testDataFrame)

'''#Display the data
class App(QWidget):
 
    def __init__(self):
        super().__init__()
        self.title = 'Proteins and their functions'
        self.left = 0
        self.top = 0
        self.width = 1200
        self.height = 800
        self.initUI()
 
    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
 
        self.createTable()
 
        # Add box layout, add table to box layout and add box layout to widget
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.tableWidget) 
        self.setLayout(self.layout) 
 
        # Show widget
        self.show()
 
    def createTable(self):
        w, h = 8, len(nameToSites);
        Matrix = [[0 for x in range(w)] for y in range(h)]

        #Create table
        self.tableWidget = QTableWidget()
        self.tableWidget.setRowCount(len(Matrix))
        self.tableWidget.setColumnCount(len(Matrix[0]))

        self.tableWidget.setItem(0,0, QTableWidgetItem("Name"))
        self.tableWidget.setItem(0,1, QTableWidgetItem("Full Name"))
        self.tableWidget.setItem(0,2, QTableWidgetItem("Alternate Name"))
        self.tableWidget.setItem(0,3, QTableWidgetItem("Lysine Modified Residues"))
        self.tableWidget.setItem(0,4, QTableWidgetItem("Key Words"))
        self.tableWidget.setItem(0,5, QTableWidgetItem("Linked to..."))
        self.tableWidget.setItem(0,6, QTableWidgetItem("Function"))
        self.tableWidget.setItem(0,7, QTableWidgetItem("Pathway"))
        
        #for i in range(1,len(nameToSites)):
        tempInt = 0
        for key, values in nameToSites.items():
        	tempInt += 1
        	self.tableWidget.setItem(tempInt,0, QTableWidgetItem(key))
        	if key in nameToFullName:
        		self.tableWidget.setItem(tempInt,1, QTableWidgetItem(nameToFullName[key]))
        	if key in nameToAlternateName:
        		self.tableWidget.setItem(tempInt,2, QTableWidgetItem(repr(nameToAlternateName[key])))
        	self.tableWidget.setItem(tempInt,3, QTableWidgetItem(repr(values)))
        	if key in nameToFunction:
        		self.tableWidget.setItem(tempInt,6, QTableWidgetItem(repr(nameToFunction[key])))
        		

        #Have a dictionary with an ID mapping to an array

        self.tableWidget.move(0,0)
 
        # table selection change
        self.tableWidget.doubleClicked.connect(self.on_click)
 
    @pyqtSlot()
    def on_click(self):
        print("\n")
        for currentQTableWidgetItem in self.tableWidget.selectedItems():
            print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = App()
    sys.exit(app.exec_())

print(namesToAccession)'''
