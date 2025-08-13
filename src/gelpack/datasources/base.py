from abc import ABC, abstractmethod
import pandas as pd


class BaseLoader(ABC):
    """one door API for users and programmes. Users shouls call .grab_data(...).
    Backends implement ._execute(spec) to turn a backend-native spec into a pd.DataFrame
    the inclusion of IN batching or LIKE chains are optional and invisible to users.

    Args:
        ABC (class): abstract base classes
    """

    def __init__(self, version, queries, default_batch_size=5000):
        self.version = version
        self.queries = queries or {}
        self.default_batch_size=5000
    

    # backend hook
    @abstractmethod
    def _execute(self, spec, **kwargs):
        """backend-native execution, works wether we use SQL strings,
        s3 URIs or REST paths... etc)
        """
        raise NotImplementedError
    
    # public / user facing functiono
    # no need to overwhelm the user with batching and such.
    # TODO: test if batching is faster than running the query in full too.