## build on the base loader to implement the labkey API

from .base import BaseLoader


class LabKeyLoader(BaseLoader):
    """
    LabKey-backed loader.

    Uses existing gelpack.gel_utils.lab_to_df under the hood.
    """

    def _execute(self, spec, **kwargs):
        from gelpack.gel_utils import lab_to_df
        return lab_to_df(sql_query=spec, dr=self.version)
    # for now just grabbing the lab_to_df function for backwards compatibility
    # but we should implement the lab_to_df function within this _execute later.
