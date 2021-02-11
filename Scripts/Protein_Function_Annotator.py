"""
Author: Nicky de Vrij
Date last update: 24/11/2020
This script annotates epitope prediction results with the protein function after preliminary data processing in R
Use reticulate to source this script in R
This script firstly reads in the proteome fasta files to create a dictionary with the protein ID as a key, and the protein function as a value.
Then, the proteinID's from the epitope prediction results are read in and compared with the dictionary, to append the protein function to the table.
"""
from Bio import SeqIO
import pandas as pd
import os

dir = os.getcwd() + '/Data/'

ldonovani_fasta_dict = SeqIO.to_dict(SeqIO.parse(dir + "LDonovaniProteome.fa", "fasta"))


def fasta_identifier_corrector():
    lbrazil_fasta_dict = {}
    unc_lbrazil_fasta_dict = SeqIO.to_dict(SeqIO.parse(dir + "LBraziliensisProteome.fa", "fasta"))
    for keys, values in unc_lbrazil_fasta_dict.items():
        lbrazil_fasta_dict[keys.replace(".", "_")[0:15]] = values
    lmajor_fasta_dict = {}
    unc_lmajor_fasta_dict = SeqIO.to_dict(SeqIO.parse(dir + "LMajorProteome.fa", "fasta"))
    for keys, values in unc_lmajor_fasta_dict.items():
        lmajor_fasta_dict[keys.replace(".", "_").replace(":", "_")[0:15]] = values
    lmexicana_fasta_dict = {}
    unc_lmexicana_fasta_dict = SeqIO.to_dict(SeqIO.parse(dir + "LMexicanaProteome.fa", "fasta"))
    for keys, values in unc_lmexicana_fasta_dict.items():
        lmexicana_fasta_dict[keys.replace(".", "_")[0:15]] = values
    linfantum_fasta_dict = {}
    unc_linfantum_fasta_dict = SeqIO.to_dict(SeqIO.parse(dir + "LInfantumProteome.fa", "fasta"))
    for keys, values in unc_linfantum_fasta_dict.items():
        linfantum_fasta_dict[keys.replace(".", "_")[0:15]] = values

    return lbrazil_fasta_dict, lmajor_fasta_dict, lmexicana_fasta_dict, linfantum_fasta_dict


lbrazil_fasta_dict, lmajor_fasta_dict, lmexicana_fasta_dict, linfantum_fasta_dict = fasta_identifier_corrector()

def protective_annotator(protdf):
    function_list = []
    for identity in protdf['Identity']:
        if identity[0:4] == 'LmxM':
            function_list.append(lmexicana_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:5] == 'LdBPK':
            function_list.append(ldonovani_fasta_dict[identity.replace(".", "_") + '.1'].description.split('.')[1].strip("1 "))
        elif identity[0:4] == 'LINF':
            function_list.append(linfantum_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:4] == 'LbrM':
            function_list.append(lbrazil_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:4] == 'LmjF':
            function_list.append(lmajor_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])

    protdf['Protein'] = function_list
    return protdf


def risk_annotator(riskdf):
    function_list = []
    for identity in riskdf['Identity']:
        if identity[0:4] == 'LmxM':
            function_list.append(lmexicana_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:5] == 'LdBPK':
            function_list.append(ldonovani_fasta_dict[identity.replace(".", "_") + '.1'].description.split('.')[1].strip("1 "))
        elif identity[0:4] == 'LINF':
            function_list.append(linfantum_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:4] == 'LbrM':
            function_list.append(lbrazil_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
        elif identity[0:4] == 'LmjF':
            function_list.append(lmajor_fasta_dict[identity.replace(".", "_")].description.split(" | ")[4][13:])
            
    riskdf['Protein'] = function_list
    return riskdf
