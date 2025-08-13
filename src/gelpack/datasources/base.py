from abc import ABC, abstractmethod
import pandas as pd
import re


class BaseLoader(ABC):
	"""one door API for users and programmes. Users shouls call .grab_data(...).
	Backends implement ._execute(spec) to turn a backend-native spec into a pd.DataFrame
	the inclusion of IN batching or LIKE chains are optional and invisible to users.

	Args:
		ABC (class): abstract base classes
	"""


	# Parameter names that, if present in a template and passed as sequences,
	# will be auto-rendered as SQL IN (...) lists or batched.
	_IN_KEYS = {"participants", "platekeys", "disease_types", "abbrs"}

	def __init__(self, version, queries, default_batch_size=5000):
		self.version = version
		self.queries = queries or {}
		self.default_batch_size=default_batch_size
	

	# backend hook
	@abstractmethod
	def _execute(self, spec, **kwargs):
		"""backend-native execution, works wether we use SQL strings,
		s3 URIs or REST paths... etc)
		"""
		raise NotImplementedError
	
	# public / user facing functions
	# no need to overwhelm the user with batching and such.
	# TODO: test if batching is faster than running the query in full too.
	def grab_data(
			self,
			key_or_spec,
			*,
			batch_param=None,
			batch_size=None,
			likes=None,
			**params
			):
		"""smart fetch 
		- if key_or_spec matches a template in self.queries ( the pre-fitted SQL queries accompanying the package) render
		placeholders using params and optional LIKE chains from 'likes', then execute the query/grab.
		- otherwise treat key_or_spec as a raw backend spec and execute it.

		LIKE chains:
		to prevent SQL injection risks we've removed f-string replacements of participants and disease terms from the
		SQL queries. Provide `likes` as {placeholder:(field_expr, terms), ...} the template/table should include those placeholders
		(e.g. {diag_like}) no default field names are assumed. 
		
		Batching:
		If the template SQL includes an IN-like placeholder (e.g. {participants}) with many values 
		(default batch_size more than 5000) the result will be concatenated accross batches.


		Args:
			key_or_spec (sql string, S3 URI or a key): what should the loader grab?
			batch_param (str, optional): what should be batched. Defaults to None.
			batch_size (int, optional): how many should be included in 1 batch. Defaults to None.
			likes (dict, optional): to build WHERE/ IN chains on specific fields. Defaults to None.
		"""

		# pull from the existing queries.
		template = self.queries.get(key_or_spec)

		# if there are no templates (e.g. a raw SQL string):
		if template is None:
			return self._execute(key_or_spec)
		
		placeholders = set(re.findall(r"{(\w+)}"), template)
		fmt = dict(params)

		# inject LIKE chains only when its requested by the template.
		# this essentially performs some of the wrangling we used to do
		# with tables in Cohort and SurvDat.
		if likes:
			for ph, (field_expr, terms) in likes.items():
				if ph in placeholders:
					fmt[ph] = self._or_like(field_expr,terms)
		
		# determine batching target:
		to_batch = batch_param or next(
			(
				k for k in self._IN_KEYS
				if k in placeholders and isinstance(fmt.get(k),(list,tuple,pd.Series))
			),
			None,
		)
		if to_batch and fmt.get(to_batch) is not None:
			ids = list(fmt[to_batch])
			bs = batch_size or self.default_batch_size

			static_fmt = {k: v for k,v in fmt.items() if k != to_batch}
			# this is where the actual batching occurs:
			return self._run_batched(
				base_sql = template,
				param_name = to_batch,
				ids=ids,
				batch_size=bs,
				extra_format=static_fmt
			)
		

		# explode small lists
		for k in (self._IN_KEYS & placeholders):
			v = fmt.get(k)
			if isinstance(v, (list, tuple, pd.Series)):
				fmt[k] = self._format_in_clause(v)

		spec = template.format(**fmt)
		return self._execute(spec)