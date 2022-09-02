"""
condense_alleles_18.py
Coded by: Nicolas Lopez-Rozo

This script reads the information of alleles of a gene and the 
expression matrix for all of the alleles and consolidates the 
expression matrix for the representative genes according to the 
allele information. Some prerequisites for making the approach 
feasible are:

- The graph of genes should correspond to a directed acyclic graph 
(DAG), where the root nodes are the representative genes
- No gene should be reachable to more than 1 representative gene

The file should be comma-separated (CS). If an allele has more than 
one gene, it should be separated with the vertical bar ('|')
"""

import os, sys, time, random
from tqdm import tqdm
import numpy as np
import pandas as pd
import networkx as nx

def test_feasibility(file):
	"""Verify the two preconditions for the algorithm to work:
	1. Graph is a DAG, whose root nodes are representative genes
	2. No gene should be reachable to more than 1 representative gene

	The function may sends a warning if any of the preconditions 
	is not met. Otherwise, it returns a triple with: 
	(i) a bool specifying whether the graph is a DAG
	(ii) a dict with representative genes as keys and subjugated 
	genes as value (CS-string) and
	(iii) a dict of representative genes for any given gene
	
	The function continues its execution despite the warning, using 
	the gene appearing first as representative.
	"""
	subjugate, repr_gene = dict(), dict()
	flag = False
	G = nx.DiGraph()
	with open(file, "r") as f:
		lines = f.readlines()
	for line in lines:
		line = line.replace("|",",").strip()
		# trick for "cleaning" up to 16 succesive commas (empty alleles)
		for i in range(4): line = line.replace(",,",",")
		name, data = line.split(",", 1)
		subjugate[name] = data
		for gene in data.split(","):
			if gene == "": continue
			if gene in subjugate:
				flag = True
				print(f"WARNING: Gene {gene} appeared as representative and as secondary\n{gene},{subjugate[gene]}\n{line.strip()}")
			elif gene in repr_gene:
				flag = True
				print(f"WARNING: Gene {gene} appeared twice\n{repr_gene[gene]},{subjugate[repr_gene[gene]]}\n{line.strip()}")				
			else:
			    G.add_edge(name,gene)
			    repr_gene[gene]=name
	if flag: print("Using the first gene as representative...")
	for gene in subjugate:
		repr_gene[gene] = gene
	return nx.is_directed_acyclic_graph(G), subjugate, repr_gene

if __name__ == '__main__':
	file_expression = "expression_data/sugar_cane_expression_2018.csv"
	file_alleles = "expression_data/alleles_2018_all.csv"
	file_results = "expression_data/Sspon18_gene_exp.csv"
	# first, check the alleles information
	isDAG, subjugate, repr_gene = test_feasibility(file_alleles)
	if not isDAG:
		raise Exception("Allele information is inconsistent (Not a DAG)")

	# second, read the expression data
	lines = None
	with open(file_expression, "r") as f:
		headers = f.readline()
		lines = f.readlines()

	# third, consolidate expression profiles
	exp_mat = [None for x in subjugate]
	names = sorted(subjugate.keys())
	pos = {v:u for u,v in enumerate(names)}
	flag = False
	for line in lines:
		name, data = line.strip().split(",", 1)
		if name not in repr_gene:
			#if not flag: print(f"Warning: Following genes not included in {file_alleles}. Included as-is:")
			#if flag: print(",", end="")
			#print(f"{name}", end="")
			flag = True
			pos[name]=len(names)
			names.append(name)
			exp_mat.append(None)
			rep = name
			#continue
		else:
			rep = repr_gene[name]
		if name == 'Sspon.ctg0212320': print(name, rep, pos[rep])
		exp = np.array([float(x) for x in data.split(",")])
		if exp_mat[pos[rep]] is None:
			exp_mat[pos[rep]] = exp
		else:
			exp_mat[pos[rep]] += exp
	if flag: print()
	print(len(subjugate), len(names))

	# finally, writing expression profiles on file_results
	with open(file_results, "w") as f:
		f.write(headers)
		for i, gene in enumerate(names):
			if exp_mat[i] is None: continue
			f.write(gene)
			for x in exp_mat[i]:
				f.write(f",{round(float(x),8)}")
			f.write("\n")

