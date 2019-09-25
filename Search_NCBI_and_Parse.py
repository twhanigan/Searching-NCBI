import numpy as np
import pandas as pd
import math 
import sys
import numpy as d
import pandas as pd
import math 
from pandas import ExcelWriter
from pandas import ExcelFile
import scipy
from scipy import stats
from scipy.stats import pearsonr
import warnings
from scipy.stats import ttest_ind
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.metrics import jaccard_similarity_score
from more_itertools import unique_everseen
import numbers
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pyplot as pyplot
import matplotlib.animation as anim
from matplotlib.animation import FFMpegWriter
from matplotlib import style
import statistics
from Bio import Entrez, SeqIO
from Bio import Entrez
import xml.etree.ElementTree as ET
from xml.etree import ElementTree
import re
import xml.etree.ElementTree as ET
import lxml.etree as et

style.use('seaborn-muted')

###Import dataset with column of gene id's
Final_Table = pd.read_excel('Final_Table_2.xlsx', dtype={'Entrez':str})

Entrez.email = "twhanigan@gmail.com"

def test_apply(x):
    try:
        return int(x)
    except ValueError:
        return None

def search_genes(id_list,search_field):
    term = " OR ".join(map(lambda x:x+"["+search_field+"]",id_list))
    esearch_result = Entrez.esearch(db="gene",term=term,retmod="xml")
    parsed_result = Entrez.read(esearch_result)
    return parsed_result['IdList']

def retrieve_annotation(id_list):
	request = Entrez.epost("gene", id = ",".join(id_list))
	try:
		result = Entrez.read(request)
	except RuntimeError as e:
		print ("An error occured while retrieving the annotations.")
		print ("The error returned was %s" % e)
		sys.exit(-1)

	webEnv = result["WebEnv"]
	queryKey = result["QueryKey"]
	data = Entrez.esummary(db = 'gene', webenv=webEnv, query_key = queryKey, retmod='xml')
	#annotations = Entrez.read(data)
	#print ("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))
	return data

#set accession from dataset to id_list
id_list = Final_Table['Entrez'].dropna().copy().tolist()
#Retrieve anotations
annotations = retrieve_annotation(id_list)
#flatten the tuple
tree = ET.parse(annotations)
root = tree.getroot()
Options = [elem.tag for elem in root.iter()]
# prints the whole document: print(ET.tostring(root, encoding = 'utf8').decode('utf8'))
#Store the elements from your data of interest(ex. 'Chromosome', or 'Name') from the options given in options output
Chrome = []
for element in root.iter('Chromosome'):
	Chrome.append(element.text)
Name = []
for name in root.iter('Name'):
	Name.append(name.text)

#Write you output from the data of interest to a new dataframe
Name_Chrome = pd.DataFrame(columns =['Name', 'Chrome'])
Name_Chrome['Name'] = Name
Name_Chrome['Chromosome'] = Chrome
#merge your data of interest output form NCBI to the original dataframe
Final_Table_Appended = Final_Table.merge(Name_Chrome, left_on= 'Site', right_on='Name',how='left')
#Remove Duplicate Values
Final_Table_Dropped = Final_Table_Appended.drop_duplicates(['Entrez']).copy()
#Replace failed string assignments with 'MT' for mitochondrial
Final_Table_Dropped['Chromosome'] = Final_Table_Dropped['Chromosome'].fillna('MT')
#change strings to numeric for sorting
Final_Table_Dropped['Chromosome'] = pd.to_numeric(Final_Table_Dropped['Chromosome'], errors = 'ignore',downcast='signed')
#write Final Table to excel
Final_Table_Dropped.to_excel("Final_Table_appended.xlsx")
#Sort By Chromosome
Final_Table_Sorted = Final_Table_Dropped.sort_values(by='Chromosome',ascending = True)

#Plot data of interest on per 'option' basis, here the option is by chromosome
fig, axes = plt.subplots()
axes.set_ylabel("Pearson", fontname="Arial", fontsize=16)
axes.set_title("Correlation: Protein Expression and Sensitivity", fontname='Arial', fontsize=18)
Final_Table_Sorted.boxplot(column=['Pearson'], by = 'Chromosome',fontsize = 16)
sns.violinplot('Chromosome', 'Pearson', data = Final_Table_Sorted, ax = axes, scale ='width', cut = 20)
plt.show()
