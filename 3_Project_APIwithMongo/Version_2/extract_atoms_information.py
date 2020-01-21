import urllib.request
from urllib.request import urlopen
import requests
import pymongo
import os.path

import pymongo
from pymongo import MongoClient
client = MongoClient('localhost', 27017)
db = client["project"]
collection = db["proteins"]

from atom import Atom
from model import Model
from chain import Chain
from aminoacid import Aminoacid

def load_protein_pdb(name_protein):
        """
        .. note:: You will download the pdb of the protein you choose.
        """
        fileDownload = "https://files.rcsb.org/download/" + name_protein + ".pdb" #Get the link
        request = requests.get(fileDownload) #Here is where im getting the status_code
        if request.status_code == 200: #Check out if a file exists
                Download = urlopen(fileDownload) #Open url
                file = open(name_protein, "wb") #Open the file, write binary
                file.write(Download.read()) #Read Download and write in file

def load_protein_fasta(name_protein):
        """
        .. note:: You will download the fasta file of the protein you choose. Fasta file will be for the function get_similar_protein
        """
        url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList=" + name_protein + "&compressionType=uncompressed"
        request = requests.get(url)
    
        if os.path.isfile(name_protein):
                if request.status_code == 200: #Check out if a file exists
                        Download = urlopen(url) #Open url
                        file = open(name_protein + '_fasta', "wb") #Open the file, write binary
                        file.write(Download.read())

def general_dictionary(name_protein): #To insert in mongodb. Its automatically
        """
        db.proteins.find({"MODEL 1.A.46.1.atom_number":{$eq: "1"}})

        db.proteins.aggregate([
                {$project:
                        {"_id":0, "MODEL 1.A.46.residue_name":0
                         }
                }
        ])
        
        :param: self
        :return: there is no return. 
        :function: insert in mongodb all the information. Protein - Models - Chains - Residue name - Residue Sequence Number - Element symbol - Atom name
        """
        model_count = 0

        protein = {}
        model_dictionary = {}

        with open(name_protein, 'r+') as text_file:
            for line in text_file:
                if line.startswith("MODEL"):
                    #MODEL
                    model_count +=1

                    Model.model_identifier = line[11:17]
                    model_identifier = Model.model_identifier.strip() #Model identifier has blanks, I want this value without them
                    model_dictionary.setdefault("MODEL " + model_identifier, {})
                    #setdefault equals the first value to the second if it doesnt exist

                if not model_count: #When we have only one model we dont have any model count neither model_identifier
                    model_identifier = "1"
                    model_dictionary["MODEL 1"]={}

                if line[:4] == 'ATOM':
                    #CHAIN
                    Chain.chain_identifier = line[21]
                    chain_identifier = Chain.chain_identifier.strip()

                    #AMINOACID
                    Aminoacid.residue_sequence_number = line[22:26]
                    residue_sequence_number = Aminoacid.residue_sequence_number.strip()

                    Aminoacid.residue_name = line[17:20]
                    residue_name = Aminoacid.residue_name.strip()

                    #ATOM
                    Atom.atom_number = line[6:11]
                    atom_number = Atom.atom_number.strip()

                    Atom.atom_name = line[12:16]
                    atom_name = Atom.atom_name.strip()

                    Atom.x_coordinate = line[30:38]
                    x_coordinate = float(Atom.x_coordinate.strip())

                    Atom.y_coordinate = line[38:46]
                    y_coordinate = float(Atom.y_coordinate.strip())

                    Atom.z_coordinate = line[46:54]
                    z_coordinate = float(Atom.z_coordinate.strip())

                    Atom.occupancy = line[54:60]
                    occupancy = float(Atom.occupancy.strip())

                    Atom.temperature_factor = line[60:66]
                    temperature_factor = float(Atom.temperature_factor.strip())

                    Atom.element_symbol = line[76:78]
                    element_symbol = Atom.element_symbol.strip()

                    model_id = protein.setdefault("MODEL " + model_identifier, {"id":model_identifier})
                    chain_letter = model_id.setdefault(chain_identifier, {"id":chain_identifier})
                    residue_number = chain_letter.setdefault(residue_sequence_number, {"id":residue_sequence_number, "residue_name":residue_name})

                    residue_number[atom_number] = {
                        "atom_number":atom_number,
                        "atom_name":atom_name,
                        "x_coordinate":x_coordinate, 
                        "y_coordinate":y_coordinate, 
                        "z_coordinate":z_coordinate, 
                        "occupancy":occupancy, 
                        "temperature_factor":temperature_factor,
                        "element_symbol":element_symbol
                    }
            collection.insert_one(protein) #Insert in mongodb
