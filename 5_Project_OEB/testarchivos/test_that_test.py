#!/usr/bin/env python3

import os , shutil
import owlready2

CachePath = '/tmp/onto'

shutil.rmtree(CachePath,ignore_errors=True)
os.makedirs(CachePath,exist_ok=True)
owlready2.onto_path.append(CachePath)

import time

w = owlready2.World(filename='/tmp/onto/firstDirect.sqlite3', exclusive=False)
print("FirstDirect")
onto1 = w.get_ontology('https://raw.githubusercontent.com/inab/OEB-ontologies/test/oebDatasets.owl').load()
w.save()
onto1.save()
res1 = onto1.search(label = "participant")
import pprint
pprint.pprint(res1)

theRes = res1[0]
print("{0} => {1}".format(theRes.label,theRes.iri))

#w = owlready2.World(filename='/tmp/onto/firstIndirect1.sqlite3', exclusive=False)
#print("FirstIndirect1")
#onto1 = w.get_ontology('https://w3id.org/oebDatasets/').load()
#w.save()
#onto1.save()
#res1 = onto1.search(label = "participant")
#import pprint
#pprint.pprint(res1)
#
#theRes = res1[0]
#print("{0} => {1}".format(theRes.label,theRes.iri))
#
#w = owlready2.World(filename='/tmp/onto/firstIndirect2.sqlite3', exclusive=False)
#print("FirstIndirect2")
#onto1 = w.get_ontology('https://w3id.org/oebDatasets').load()
#w.save()
#onto1.save()
#res1 = onto1.search(label = "participant")
#import pprint
#pprint.pprint(res1)
#
#theRes = res1[0]
#print("{0} => {1}".format(theRes.label,theRes.iri))
#
#########res1b = onto1.search(iri = "*EFO_0003042")
#########pprint.pprint(res1b)
#########res1c = onto1.search(iri = "*EFO_000304")
#########pprint.pprint(res1c)
#########
#########print("Second")
#########onto2 = w.get_ontology('http://purl.obolibrary.org/obo/cl/releases/2018-07-07/cl.owl').load()
#########w.save()
#########onto2.save()
#########
#########print("Third")
#########onto3 = w.get_ontology('http://purl.obolibrary.org/obo/uberon/releases/2018-07-30/uberon.owl').load()
#########w.save()
#########onto3.save()
#########
#########print("Fourth")
#########onto4 = w.get_ontology('http://purl.obolibrary.org/obo/obi/2018-08-27/obi.owl').load()
#########w.save()
#########onto4.save()
#########res4 = onto4.search(label = "ChIP-seq assay")
#########pprint.pprint(res4)
#########
#########res4b = res4[0].ancestors()
#########pprint.pprint(res4b)

# http://www.ebi.ac.uk/efo/EFO_0003042
