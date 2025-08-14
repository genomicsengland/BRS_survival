# pattern match and plug the right loader in.

from .labkey_loader import LabKeyLoader

def get_loader(
		api, 
		version, 
		queries_override=None, 
		allow_env_override=True,
		**kwargs):
	
	if api=="labkey":
		
		### select Cohort and Survdat queries ###
		## this needs to be updated on new releases.
		from gelpack.queries.labkey_queries_v19 import QUERIES as DEFAULTS

		# Build a list so updated entries can override establishe defaults
		qsrc = [DEFAULTS]
		if queries_override is not None:
			# Accept list/tuple or single source
			if isinstance(queries_override, (list, tuple)):
				qsrc.extend(queries_override)
			else:
				qsrc.append(queries_override)

		
		### pass through the Labkey loader ###

		return LabKeyLoader(
			version=version,
			queries=qsrc,
			allow_env_override=allow_env_override,
			**kwargs
		)
	raise ValueError(f"Unknown API: {api}")
