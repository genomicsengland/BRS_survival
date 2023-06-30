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
		print(mapping)
		return dataframe, mapping
	else:
		raise ValueError(f'Unrecognized type: {type}.')

	return dataframe


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