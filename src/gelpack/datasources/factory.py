# pattern match and plug the right loader in.
import os
import importlib.util
from .labkey_loader import LabKeyLoader


def get_loader(
		api, 
		version, 
		queries_override=None, 
		allow_env_override=True,
        queries_env_var="GELPACK_LABKEY_QUERIES",
		**kwargs):
	
	
	if api=="labkey":
		
		### select Cohort and Survdat queries ###
		## this needs to be updated on new releases.
		from gelpack.queries.labkey_queries_v19 import QUERIES as DEFAULTS

		# as we dont want to have to update the package every time a minor change
		# occurs in the back-end structures, we want to be able to pass a path to
		# either a completely new query file, or a query file that contains updated
		# queries:
		env_q = {}
		if allow_env_override:
			path = os.getenv(queries_env_var)
			if path:
				env_q = _load_queries_from_path(path)
        
		# explicit override passed by caller
		user_q = {}
		if isinstance(queries_override, (list, tuple)):
			for d in queries_override:
				if d:
					if not isinstance(d, dict):
						raise TypeError("Each item in queries_override must be a dict")
					user_q.update(d)
		elif isinstance(queries_override, dict):
			user_q = queries_override
		elif queries_override is not None:
			raise TypeError("queries_override must be a dict or list/tuple of dicts")

		queries = _merge_query_dicts(DEFAULTS, env_q, user_q)

		### pass through the Labkey loader ###
		return LabKeyLoader(version=version, queries=queries)

	raise ValueError(f"Unknown API: {api}")


def _load_queries_from_path(path):
    spec = importlib.util.spec_from_file_location("user_queries", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    q = getattr(mod, "QUERIES", None)
    if not isinstance(q, dict):
        raise ValueError(f"{path} does not define a dict named QUERIES")
    return q


def _merge_query_dicts(*dicts):
    merged = {}
    for d in dicts:
        if not d:
            continue
        if not isinstance(d, dict):
            raise TypeError("queries_override items must be dicts")
        merged.update(d)  # later dicts override earlier ones
    return merged