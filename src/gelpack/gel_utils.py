import pandas as pd
import numpy as np
import labkey

def lab_to_df(sql_query, dr):
	"""generate an pandas dataframe from labkey sql query

	Args:
		sql_query (str): an sql query as string.
		dr (str): GEL datarelease version
	"""
	import pandas as pd
	import warnings
	import labkey

	ver = labkey.__version__

	if ver == '1.2.0' or ver == '1.4.0':
		server_context = labkey.utils.create_server_context(
			domain= "labkey-embassy.gel.zone",
			container_path = dr,
			context_path = "labkey",
			use_ssl = True
		)
		
		results =  labkey.query.execute_sql(
			server_context,
			schema_name="lists",
			sql=sql_query,
			max_rows=50000000
		)

	if ver == '2.4.0':
		from labkey.api_wrapper import APIWrapper

		labkey_server = "labkey-embassy.gel.zone"
		project_name = dr  # Project folder name
		contextPath = "labkey"
		schema = 'lists'
		api = APIWrapper(
			labkey_server, 
			project_name, 
			contextPath, 
			use_ssl=True)
		
		results = api.query.execute_sql(sql=sql_query, schema_name=schema)
	
	return(pd.DataFrame(results['rows']))



def assign_groups(dataframe, vars, type='and'):
	"""assigns groups/labels for survival analysis based on a list of variables.

	Args:
		dataframe (pd.DataFrame): A pandas dataframe with samples / participants 
			their survival time, censoring and variables to stratify on.
		vars (list): list of variable names in string format, corresponding to 
			columns of boolean series (if a var is True or False in a sample)
			in dataframe. 
		type (str, optional): determines how the groups are assigned:
			'and': binary groups (1, 2) where var1 AND var2 AND..var* are True.
			'or': binary groups where var1 OR var2 OR var* are True.
			'full': groups on each possible iteration of vars. 
			Defaults to 'and'.
	Returns:
		pd.series of grouping assigned as group_1, group_2, group_n
	"""
	# check if each var is indeed in dataframe
	dataframe=dataframe.copy()
	for var in vars:
		if not var in dataframe:
			raise ValueError(
				f'{var} should only include columns in dataframe.'
				)
		# if dataframe[var].isna().any():
			# raise ValueError(
			# 	f'{var} contains NA values, only True/False allowed'
			# 	)
	# dealing with np.nan values:
	# this depends on the type
	# 'and' = any of the vars group = unknown
	# 'or' = all of the vars group = unknown
	# 'full' = no need to exclude the nan values.
	if type == 'and':
		nan_frame = pd.DataFrame()
		nan_frame = dataframe.loc[dataframe[vars].isna().any(axis=1)]
		nan_frame['group'] = 'unknown'
		dataframe.dropna(subset=vars, how='any', axis=0, inplace=True)
	if type == 'or':
		nan_frame = pd.DataFrame()
		nan_frame = dataframe.loc[dataframe[vars].isna().all(axis=1)]
		nan_frame['group'] = 'unknown'
		dataframe.dropna(subset=vars, how='all', axis=0,inplace=True)
	
	# setting up the grouping
	if type=='and':
		dataframe['group'] = (dataframe[vars]
			.all(axis=1)
			# here we could assign group names instead. taken from Strata.
			# although, it would be neater if the group names can be assigned
			# through the mapping file or a dictionary?
			.replace({True: 'group_1', False: 'group_2'})  
		)
		dataframe = pd.concat([dataframe,nan_frame])
		# return dataframe('group')
	elif type=='or':
		dataframe['group'] = (dataframe[vars]
			.any(axis=1)
			.replace({True: 'group_1', False: 'group_2'})
		)
		# return dataframe['group']
		dataframe = pd.concat([dataframe,nan_frame])
	elif type=='full':
		dataframe['combination'] = dataframe[vars].agg(tuple, axis=1)
		dataframe['group'] = dataframe['combination'].factorize()[0]
		# dataframe.groupby(['snv1','snv2'], sort=False).ngroup()
		mapping = dataframe[['combination','group']].drop_duplicates()
		return dataframe, mapping
	else:
		raise ValueError(f'Unrecognized type: {type}.')

	return dataframe 




# TODO test allowing of NAN (will they be included in the groups?)
# TODO move this to gel_utils (and adjust script to search for it there.)
# string joining based on the True / False combinations.
# all False should be WT, handled seperately as the zip will be empty.
def create_name_map(
	mapping,
	genes):
	"""Create a dicitonary to rename groups generated by assign_groups
		
	Args:
		mapping (pd.dataframe): a 2 column dataframe with a (boolean) tuple
			matching the order of the variable genes.
		genes (list): list of genes being queried, order must match the 
			order of the mapping, as well as the order of the columns
			in the dataframe on which assign_groups was applied.
	Returns:
		dictionary with names (item) associated with group (key).
	"""
	name = []
	for map in mapping.iloc[:,0]:
		if not any(map):
			name.append('WT')
		else:
			name.append(f'''{
				'_'.join(
					[x for x,y in zip(tuple(genes), map) if y]
					)
				}'''
			)
	return dict(zip(mapping.iloc[:,1], name))


def force_list(iter):
	# a test to make sure the gene is in a list.
	# so that a single gene is not parsed like 'K', 'R', 'A', 'S'
	try:
		from collections.abc import Iterable
	except ImportError:
		from collections import Iterable

	if not isinstance(iter, Iterable) or isinstance(iter, str):
		iter = [iter]
	return iter




def translateicd(icd_vec, lookup=None):
	"""As we are importing ICD-10 codes from different RWE sources their 
	formatting has to be unified. this function cleans the codes and returns
	a matching disease type. 
	

	Args:
		icd_vec (list): list of str, ICD-10 codes from different sources.
		icd_dictionary (DataFrame): reference dictionary of disease types
			and regex-form ICD-10 codes
	"""
	import re
	if lookup is None:
		lookup = [
			('ADULT_GLIOMA',r"C71[0-9]{0,1}|D43[0-4]{0,1}|D32[0-4]{0,1}|D33[0-4]{0,1}"),
			('BLADDER',r"C67[0-9]{0,1}|D090|D414"),
			('BREAST',r"C50[0-9]{0,1}|D05[0-9]{0,1}"),
			('CARCINOMA_OF_UNKNOWN_PRIMARY',r"C80[0-9]{0,1}|D489"),
			# ('OTHER',r"C80[0-9]{0,1}|D489"),  # putting unknown in OTHER, like childhood.
			('COLORECTAL',r"C18[0-9]{0,1}|C20|C19|D010|D012|D011|C17[0-9]{0,1}"
				r"|C21[0-9]{0,1}|C78$|C784|C785|C788"),  # changed from 
			('ENDOCRINE',r"C73|C74[0-9]{0,1}|C75[0-9]{0,1}|D093|D44[0-9]{0,1}"),
			('ENDOMETRIAL_CARCINOMA',r"C53[0-9]{0,1}|C54[0-9]{0,1}"
				r"|C55[0-9]{0,1}|D070|D06[0,1,7,9]{0,1}"),
			('HAEMONC',r"C81[0-9]{0,1}|C82[0-9]{0,1}|C83[0-9]{0,1}|C84[0-9]{0,1}"
				r"|C85[0-9]{0,1}|C86[0-9]{0,1}|C88[0-9]{0,1}|C90[0-9]{0,1}|C91[0-9]{0,1}"
				r"|C92[0-9]{0,1}|C93[0-9]{0,1}|C94[0-9]{0,1}|C95[0-9]{0,1}|C96[0-9]{0,1}|"
				r"D47[0-9]{0,1}"),
			('HEPATOPANCREATOBILIARY',r"C22[0-9]{0,1}|C23|C24[0-9]{0,1}"
				r"|C25[0-9]{0,1}|C787|D015|D376"),
			('LUNG',r"C34[0-9]{0,1}|D022|C450|C459"),
			('MALIGNANT_MELANOMA',r"C43[0-9]{0,1}|D03[0-9]{0,1}"),
			('NASOPHARYNGEAL',r"C11[0-9]{0,1}"),
			('ORAL_OROPHARYNGEAL',r"C01[0-9]{0,1}|C02[0-9]{0,1}|C03[0-9]{0,1}|C04[0-9]"
				r"{0,1}|C05[0-9]{0,1}|C06[0-9]{0,1}|C07[0-9]{0,1}|C08[0-9]{0,1}"
				r"|C09[0-9]{0,1}|C10[0-9]{0,1}|C00[0-9]{0,1}|C14[0-9]{0,1}|D000|D370"),
			('OVARIAN',r"C56|C796|D391|C57[0-9]{0,1}"),
			('PROSTATE',r"C61|D075"),
			('RENAL',r"C64|C65|C66|C68[0-9]{0,1}|C790|D410|D411|D412"),
			('SARCOMA',r"C41[0-9]{0,1}|C45\b|C45[1-8]{1}|C46[0-9]{0,1}|C48[0-9]{0,1}"
				r"|C49[0-9]{0,1}|C40[0-9]{0,1}|C786|D481|D483|D484"),
			('SINONASAL',r"C12|C30[0-9]{0,1}|C31[0-9]{0,1}"),
			('TESTICULAR_GERM_CELL_TUMOURS',r"C62[0-9]{0,1}"),
			('UPPER_GASTROINTESTINAL',r"C15[0-9]{0,1}|D001|C16[0-9]{0,1}"
				r"|D002|C26[0-9]{0,1}|D017|D019|D379")
		]

	def transicd(icd, lookups=lookup):
		for value, pattern in lookups:
			if re.search(pattern, icd):
				return value
		return 'OTHER'
	# remove trailing 'X' and any '.'
	icd_vec_clean = [re.sub('X$|\[.\]|[.]' , '', str(x)) for x in icd_vec]
	return list(map(transicd, icd_vec_clean))
