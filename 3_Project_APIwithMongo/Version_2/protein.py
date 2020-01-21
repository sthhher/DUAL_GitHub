import requests
import urllib.request
from urllib.request import urlopen
from xml.dom import minidom
from Bio.SeqUtils import seq1
from Bio import SeqIO
from itertools import groupby
import string
import os.path as path
import string

import pymongo
from pymongo import MongoClient
client = MongoClient('localhost', 27017)
db = client["project"]
collection = db["proteins"]

import extract_atoms_information
from atom import Atom
from model import Model
from chain import Chain
from aminoacid import Aminoacid

from Bio.PDB import PDBParser

class Protein:

    def __init__(self, protein_name):
        self.protein_name = protein_name
        """
        for model in db.get_models(self.protein_name):
            db.proteins.aggregate([{$project:{"_id":1, Model.model_identifier:0}}])
            self.model_list.append(Chain(protein_name, chain_identifier))
        """
        extract_atoms_information.load_protein_pdb(self.protein_name) #It creates automatically so, we dont have to do it. #Protein_pdb
        extract_atoms_information.load_protein_fasta(self.protein_name) #Protein_fasta
        
        try: #Check it if exists in mongodb to insert
            extract_atoms_information.general_dictionary(self.protein_name)
        except: 
            db.protein.find({"_id":self.protein_name}) == True
        
    def model_list(self):
        pass

    def get_sequence_aminoacids (self): #Returns a list with chain and the sequence
        """
        :param: self
        :return: sequences
        :rtype: chain dictionary with aminoacids (one letter)
        """
        sequences = {}
        open_fasta = list(SeqIO.parse(self.protein_name + "_fasta", 'fasta')) #Show all sequences with information in the file
        chain = (len(open_fasta)-1) #Counts how many chains there are, but we have to keep in mind the 0. We subtract 1 
        chains = self.get_chain_list()
        while chain >=0:
            sequence = open_fasta[chain].seq
            sequences[chains[chain]] = str(sequence)
            chain -=1
        return sequences

    def get_chain_list(self): #Return the differents chains there are
        """
        :param: self
        :return: list_differents_chains
        :rtype: a list with all the chains that may have
        """
        list_differents_chains = []
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    Chain.chain_identifier = line[21]
                    chain_identifier = Chain.chain_identifier.strip()
                    if chain_identifier not in list_differents_chains:
                        list_differents_chains.append(chain_identifier)
            return list_differents_chains

    def get_aminoacid_list(self): #Return the differents aminoacids with the differents sequence number classified in chains
        """
        :param: self
        :return: protein
        :rtype: a dictionary with all the chains, inside of that another dictionary with all the aminoacids that exist and inside of aminoacids a list with all the sequence number that have
        """
        protein = {}
        with open(self.protein_name, 'r+') as text_file:
            for line in text_file:
                if line.startswith("ATOM"):
                    Chain.chain_identifier = line[21]
                    chain_identifier = Chain.chain_identifier.strip()

                    Aminoacid.residue_name = line[17:20]
                    residue_name = Aminoacid.residue_name.strip()

                    Aminoacid.residue_sequence_number = line[22:26]
                    residue_sequence_number = Aminoacid.residue_sequence_number.strip()
                    
                    chain_letter = protein.setdefault(chain_identifier, {})
                    residue_letter = chain_letter.setdefault(residue_name, [])
                    if residue_sequence_number not in residue_letter:
                        residue_letter.append(residue_sequence_number)
        return protein

    def get_similar_protein(self): 
        """
        :param: self
        :return: dicc_hit_def
        :rtype: a dictionary with all the chains as a key, and the value is the most similar protein
        """
        #SALE LA MISMA PROTEINA QUE INTRODUZCO EN ALGUNAS (EJ. 2F40)
        dicc_hit_def = {}
        open_fasta = list(SeqIO.parse(self.protein_name + "_fasta", 'fasta')) #Show all sequences with information in the file
        chain = (len(open_fasta)-1) #Counts how many chains there are, but we have to keep in mind the 0. We subtract 1 
        while chain >=0:
            fasta_sequence = open_fasta[chain].seq #We want the sequence
            url = "https://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=" + str(fasta_sequence) + "&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML" #TENGO QUE HACER QUE ME DIGA EL NOMBRE DE LA PROTEINA
            request = requests.get(url) #Here is where im getting the status_code
            file_protein_pdb = "similar_protein_" + self.protein_name #The name that the file will have

            if request.status_code == 200: #Check out if a file exists
                Download = urlopen(url) #Download pdb
                file = open(file_protein_pdb, "wb") #I want to write in bytes in file_protein_pdb
                file.write(Download.read()) 
                file.close()
                doc = minidom.parse(file_protein_pdb) #This is for read the xml
                hits = doc.getElementsByTagName("Hit") #Keep what find in the tag hit in hits.
                if hits != []: #There are wrong chains so they are empty because of there isnt "hit"
                    hit_def = hits[0].getElementsByTagName("Hit_def")[0] #Define Hit_def
                    hit_def = hit_def.firstChild.data #I want values of hit_def
                    hit_def = hit_def[0:4] #I only want 4 first values of hit_def
                    dicc_hit_def[open_fasta[chain].id[5]]=hit_def #Add to the dictionary
                    chain = chain-1
                else: 
                    dicc_hit_def[open_fasta[chain].id[5]] = None
                    chain = chain-1
        return dicc_hit_def

Protein("2ki5")