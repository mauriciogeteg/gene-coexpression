import os, sys, time, random
from tqdm import tqdm
import numpy as np
import pandas as pd

if __name__ == '__main__':
	file_expression = "expression_data/sugar_cane_allele_2019.csv"
	file_results = "expression_data/sugar_cane_expression_2019.csv"
	with open(file_expression, "r") as f:
		headers = f.readline().strip()
		lines = f.readlines()
	
	last_name, last_exp = "", []
	with open(file_results, "w") as f:
		f.write(headers)
		for i, line in tqdm(enumerate(lines), desc="writing file"):
			name, data = line.strip().split(",", 1)
			stem = name.strip().split("-")[0]
			if stem == name: raise Exception(f"CHECK GENE DESCRIPTORS! {name} on line {i}")
			exp = np.array([float(x) for x in data.split(",")])
			if stem == last_name:
				last_exp += exp
			else:
				# write current data on the file
				f.write(last_name)
				for x in last_exp:
					f.write(f",{float(x):0.6f}")
				f.write("\n")
				# replace last_name and last_exp with current data
				last_name, last_exp = stem, exp
		f.write(last_name)
		for x in last_exp:
			f.write(f",{round(float(x),8)}")
		f.write("\n")
		
