#### Christian J. Bouwens
#### Matthieu Vizuete-Forster
#### Miruna Carmen Barbu
#### Chris Odhams
#### BRSC team
#### generate a cohort from HPO terms, ICD-10 codes, cancer study abbreviations
# 	 and/or well as disease types
#### can be applied for cancer as well as rare disease cases.
#### last update: 2025.09.15


import warnings
from functools import reduce

import pandas as pd
import numpy as np


from gelpack.gel_utils import force_list
from gelpack.survival import Survdat


# grab a default LabKey loader if the user doesn't inject one
try:
    from gelpack.datasources.factory import get_loader as _get_loader
except Exception:
    _get_loader = None



class Cohort(object):
    """build and manage a cohort

    Inputs can be:
        - featdict with keys in {'icd10','hpo','terms','cancer_terms','cancer_abbr'}
        - explicit participant IDs
        - sample platekeys

    State (populated lazily by methods):
        - self.pids: Dict[str, pd.Series]         # source -> participant_ids
        - self.feature_tables: Dict[str, Dict]    # e.g. {'age': {...}, 'sex': {...}, ...}
        - self.sample_tables: Dict[str, Dict]     # e.g. {'cancer_samples': {...}, ...}
        - self.all_data: pd.DataFrame             # merged long view (participant-keyed merge; sample cols can duplicate PIDs)
        - self.all_cancer_samples: pd.DataFrame   # sample-level cancer rows (if fetched)
        - self.all_rd_samples: pd.DataFrame       # sample-level RD rows (if fetched)
        - self.participant_view: pd.DataFrame     # normalized 1-row-per-participant view (first-match aggregation)
        - self.sample_view: pd.DataFrame          # normalized 1-row-per-sample view with unified `sample_id`
        - self.provenance: List[dict]             # decisions log (e.g., which platekeys were dropped and why)
    """
    def __init__(
        self,
        version=None,
        loader=None,
        featdict=None,
        participants=None,
        platekeys=None,
        name=None,
    ):
        # Loader setup
        if loader is None:
            if version is None:
                raise ValueError("Provide a loader or a LabKey data-release `version`.")
            if _get_loader is None:
                raise RuntimeError(
                    "Could not import gelpack.datasources.factory.get_loader; "
                    "pass an explicit loader=... or ensure the package structure is available."
                )
            loader = _get_loader(api="labkey", version=version)

        self.loader = loader
        self.version = version or getattr(loader, "version", None)

        # Optional: initial feature dictionary
        self.featdict = featdict or {}
        if self.featdict:
            # validate + normalise keys and cache convenience lists
            valid_keys = {"icd10", "hpo", "terms", "cancer_terms", "cancer_abbr"}
            if not any(k in valid_keys for k in self.featdict.keys()):
                raise KeyError(
                    "featdict must contain at least one of: icd10, hpo, terms, cancer_terms, cancer_abbr"
                )
            if "icd10" in self.featdict:
                self.icd10s = force_list(self.featdict["icd10"])
            if "hpo" in self.featdict:
                self.hpo = force_list(self.featdict["hpo"])
            if "terms" in self.featdict:
                self.dterms = force_list(self.featdict["terms"])
            if "cancer_terms" in self.featdict:
                self.cterms = force_list(self.featdict["cancer_terms"])
            if "cancer_abbr" in self.featdict:
                self.cabbr = force_list(self.featdict["cancer_abbr"])

        # Core state
        self.pids = {}                 # source_key -> Series[participant_id]
        self.feature_tables = {}       # feature_name -> {source_key: DataFrame}
        self.sample_tables = {}        # sample_name  -> {source_key: DataFrame}
        self.name = name
        self.platekeys = None          # union/derived platekeys
        self.input_platekeys = None    # exactly what the user passed (if any)
        self.plate_map = None          # DataFrame[plate_key, participant_id]
        self.groups = {}
        self.survdat = None
        self.provenance = []           # decisions log

        # Require at least one seed
        if participants is None and platekeys is None and not self.featdict:
            raise ValueError("Provide featdict and/or participants and/or platekeys to create a cohort.")

        # Participants intake (explicit list or "all_cancer" convenience)
        if participants is not None:
            if isinstance(participants, str) and participants == "all_cancer":
                ca_pids = self.loader.grab_data("cancer_parts")
                if ca_pids is not None and "participant_id" in ca_pids.columns:
                    self.add_participants(ca_pids["participant_id"])
                else:
                    warnings.warn("No participants found for 'all_cancer'.")
            else:
                self.add_participants(participants)

        # Platekeys intake -> resolve to participants
        if platekeys is not None:
            pk_series = pd.Series(force_list(platekeys), dtype=str).dropna().drop_duplicates()
            plate_pids = self.loader.grab_data("pids_by_platekeys", platekeys=pk_series)
            if plate_pids is None or plate_pids.empty:
                warnings.warn("No participants found for provided platekeys.")
            else:
                # Store mapping + platekeys + participants
                self.plate_map = plate_pids[["plate_key", "participant_id"]].drop_duplicates()
                self.platekeys = pk_series
                self.input_platekeys = pk_series  # remember the exact input set
                self.add_participants(plate_pids["participant_id"])


    # small helper function to ensure participant_id are strings.
    @staticmethod
    def _ensure_pid_str(df):
        """Ensure participant_id column uses pandas 'string' dtype."""
        if isinstance(df, pd.DataFrame) and "participant_id" in df.columns:
            df = df.copy()
            df["participant_id"] = df["participant_id"].astype("string")
        return df

    @staticmethod
    def _series_as_pid(s):
        """Coerce any list/Series to de-duplicated 'string' Series named participant_id."""
        if isinstance(s, pd.Series):
            out = s.dropna().astype("string").rename("participant_id")
        else:
            out = pd.Series(force_list(s), dtype="string", name="participant_id").dropna()
        return out.drop_duplicates().reset_index(drop=True)
    ## ---------------------------------------------------------------##
    ########## adding or removing participants to cohort ###############
    ## ---------------------------------------------------------------##

    # enable featdicts to be added post initialisation of the cohort.
    def set_featdict(self, featdict, merge=True):
        """
        (Optionally) merge a new featdict into this cohort and normalize cached lists.
        If merge=False the featdict is replaced.
        """
        if not isinstance(featdict, dict) or not featdict:
            raise ValueError("set_featdict: provide a non-empty dict.")

        valid_keys = {"icd10", "hpo", "terms", "cancer_terms", "cancer_abbr"}
        unknown = [k for k in featdict.keys() if k not in valid_keys]
        if unknown:
            raise KeyError(f"Unknown featdict keys: {unknown}. Allowed: {sorted(valid_keys)}")

        if merge and isinstance(self.featdict, dict):
            # merge lists with de-dup
            merged = dict(self.featdict)
            for k, v in featdict.items():
                prev = merged.get(k, [])
                merged[k] = list(pd.unique(force_list(prev) + force_list(v)))
            self.featdict = merged
        else:
            self.featdict = {k: force_list(v) for k, v in featdict.items()}

        self._normalize_featdict_keys()

        self.provenance.append({
            "action": "set_featdict",
            "merge": bool(merge),
            "keys": sorted(self.featdict.keys())
        })

    def _normalize_featdict_keys(self):
        """Cache axis lists from self.featdict onto attributes for quick access."""
        d = self.featdict if isinstance(self.featdict, dict) else {}
        # Clear then (re)populate to reflect current featdict
        for attr in ("icd10s", "hpo", "dterms", "cterms", "cabbr"):
            if hasattr(self, attr):
                delattr(self, attr)

        if "icd10" in d:
            self.icd10s = force_list(d["icd10"])
        if "hpo" in d:
            self.hpo = force_list(d["hpo"])
        if "terms" in d:
            self.dterms = force_list(d["terms"])
        if "cancer_terms" in d:
            self.cterms = force_list(d["cancer_terms"])
        if "cancer_abbr" in d:
            self.cabbr = force_list(d["cancer_abbr"])

    def _ensure_axis(self, feat_key, attr_name):
        """
        Ensure e.g. self.dterms / self.icd10s exists; if missing, try to hydrate
        from self.featdict. Raise a clear error if still unavailable.
        """
        vals = getattr(self, attr_name, None)
        if vals:
            return  # already cached

        if isinstance(self.featdict, dict) and feat_key in self.featdict:
            setattr(self, attr_name, force_list(self.featdict[feat_key]))
            return

        raise ValueError(
            f"No '{feat_key}' provided. Pass featdict={{'{feat_key}': [...]}} "
            f"at construction or call set_featdict(...) before this method."
        )

    ## ------------------------------------------------
    # split old costum_pids into several functions:
    # _read_participants_from_file
    # add_participants
    # remove_participants
    # filter_participants
    ## ------------------------------------------------
    def _read_participants_from_file(self, path):
        import csv

        with open(path, "r") as f:
            try:
                csv.Sniffer().has_header(f.read(1024))
            except csv.Error:
                raise ValueError("Participants file has no header.")
            dialect = csv.Sniffer().sniff(f.read(1024), [",", "\t"])
            if dialect.delimiter not in [",", "\t"]:
                raise ValueError("Unsupported delimiter (use ',' or '\\t').")
            f.seek(0)
            reader = csv.DictReader(f, dialect=dialect)
            lines = [line for line in reader]
            if not lines or "participant_id" not in lines[0].keys():
                raise ValueError("Header must include 'participant_id'.")
        part_table = pd.read_csv(path, sep=dialect.delimiter)
        
        return part_table["participant_id"]
    
    def add_participants(self, ids):
        """add participant_id to a cohort membership (self.pids['custom']).
        This does not affect cached tables - call concat/ features again as needed.

        Args:
            ids (str, list, pd.Series): participant_ids to add.
        """

        if isinstance(ids, str) and ids.endswith((".tsv", ".csv", ".txt")):
            ids = self._read_participants_from_file(ids)

        ids = self._series_as_pid(ids)  # <-- canonicalize here

        prev = self.pids.get("custom")
        if prev is not None:
            prev = self._series_as_pid(prev)
            combined = (pd.concat(
                [prev, ids], 
                ignore_index=True)
            .drop_duplicates()
            .reset_index(drop=True))
            added_n = int(len(combined) - len(prev))
            self.pids["custom"] = combined
        else:
            self.pids["custom"] = ids
            added_n = int(len(ids))

        self.provenance.append({
            "action": "add_participants",
            "n_added": added_n,
            "total_custom": int(len(self.pids["custom"]))
        })


    def remove_participants(self, ids):
        """remove participant_id from cohort membership, also purging their 
        data from the cached feature tables. equivalent to old costum_participant(..., action='exclude')
        """
        if isinstance(ids, str) and ids.endswith((".tsv", ".csv")):
            ids = self._read_participants_from_file(ids)
        elif isinstance(ids, list):
            ids = pd.Series(ids, name="participant_id")
        elif isinstance(ids, pd.Series):
            ids = ids.rename("participant_id")
        else:
            warnings.warn("remove_participants: unsupported type; ignoring.")
            return

        id_set = set(self._series_as_pid(ids).astype(str))

        # membership
        for key, series in list(self.pids.items()):
            self.pids[key] = pd.Series([x for x in series if str(x) not in id_set], name="participant_id")

        # features
        for name, feature in list(self.feature_tables.items()):
            if isinstance(feature, dict):
                for k, tbl in list(feature.items()):
                    mask = ~tbl["participant_id"].astype(str).isin(id_set)
                    self.feature_tables[name][k] = tbl.loc[mask].copy()

        # samples
        for name, smp in list(self.sample_tables.items()):
            for k, tbl in list(smp.items()):
                mask = ~tbl["participant_id"].astype(str).isin(id_set)
                self.sample_tables[name][k] = tbl.loc[mask].copy()

        # all_data & normalized views
        if hasattr(self, "all_data") and isinstance(self.all_data, pd.DataFrame):
            self.all_data = self.all_data.loc[~self.all_data["participant_id"].astype(str).isin(id_set)].copy()
        if hasattr(self, "participant_view"):
            self.participant_view = self.participant_view.loc[
                ~self.participant_view["participant_id"].astype(str).isin(id_set)
            ].copy()
        if hasattr(self, "sample_view"):
            self.sample_view = self.sample_view.loc[~self.sample_view["participant_id"].astype(str).isin(id_set)].copy()

        self.provenance.append({
            "action": "remove_participants", 
            "n": len(id_set)
            })
        
    def filter_participants(self, query_string, mode='exclude',commit=True):
        """Run a pandas query on self.all_data to select participants, then 
        either remove the matched participant_ids, or keep only the matched participant_id 
        (remove the inverse).
        If commit = False: returns, but does not mutate the cohort.

        Args:
            query_string (str): pd.query like string
            mode (str, optional): exclude = remove matched participant_ids. 
                include = keep matched participant_ids. Defaults to 'exclude'.
            commit (bool, optional): if false, return (matched_df, matched_pids) 
                and does not mutate the cohort. Defaults to True.
        """
        if not hasattr(self, "all_data"):
            self.concat_all()
        df = self.all_data.query(query_string, engine="python")

        matched_pids = pd.Series(df["participant_id"].dropna().astype(str).unique(), name="participant_id")
        if not commit:
            return df, matched_pids

        if mode == "exclude":
            self.remove_participants(matched_pids)
            self.provenance.append({"action": "filter_participants", "mode": "exclude", "query": query_string})
        elif mode == "include":
            all_pids = pd.Series(self.all_data["participant_id"].astype(str).unique())
            to_drop = all_pids[~all_pids.isin(matched_pids)]
            self.remove_participants(to_drop)
            self.provenance.append({"action": "filter_participants", "mode": "include", "query": query_string})
        else:
            raise ValueError("mode must be 'exclude' or 'include'.")

        # refresh merged views after membership change
        self.concat_all()
        self.build_normalized_views()

    ## ----------------------------------------------------------
    # Groups
    ## ----------------------------------------------------------
    def add_to_group(self, id, group):
        if group not in self.groups:
            self.groups[group] = []
        # confirm the id is in the cohort.
        if any(id in series.values for series in self.pids.values()):
            self.groups[group].append(id)

    ## ---------------------------------------------------------------##
    # link to survival analysis
    ## ---------------------------------------------------------------##
    def initialize_survdat(self, impute=False,df=None):
        """
        Lazily initialize the Survdat object when all necessary data is available.
        """
        self.concat_all()

        if not hasattr(self, "all_pids") or len(self.all_pids) == 0:
            raise ValueError("PIDs must be set before initializing Survdat.")
        
        self.survdat = Survdat(
            df=df,
            pids=self.all_pids,
            version=self.version,
            impute=impute)

    ## ---------------------------------------------------------------##
    ## Query helpers (participant_id as key)
    ## ---------------------------------------------------------------##
    def get_cancer_parts(self):
        return self.loader.grab_data("cancer_parts")

    def get_pids_per_platekey(self, platekeys):
        return self.loader.grab_data("pids_by_platekeys", platekeys=platekeys)


    ## ---------------------------------------------------------------##
    ## PID  (participant_id as key)
    ## ---------------------------------------------------------------##
    def get_term_pids(self):
        """get participant_ids associated with the given disease_terms from
        various labkey tables: 
            - rare_diseases_participant_disease
        This function does not look for cancer disease terms.

        Returns:
            pd.DataFrame: a dataframe with participant_ids and normalised
            disease terms.
        """
        self._ensure_axis("terms", "dterms")
        df = self.loader.grab_data("rd_participant_disease_by_terms", terms=self.dterms)
        df = self._ensure_pid_str(df)
        self.dterm_table = df
        self.pids["dterm"] = df["participant_id"]
    

    def get_hpo_pids(self):
        """get participant_ids, normalised_hpo_id and terms associated with the 
        given hpo terms from the labkey table:
            - rare_diseases_participant_phenotype

        Returns:
            pd.DataFrame: a dataframe with participant_ids and normalised
            disease terms.
        """
        self._ensure_axis("hpo", "hpo") 
        df = self.loader.grab_data("hpo_pids_by_ids", hpo_ids=self.hpo)
        df = self._ensure_pid_str(df)
        self.hpo_table = df
        self.pids["hpo"] = df["participant_id"]


    def get_cterm_pids(self):
        """get participant_ids, cancer_disease_types associated with the 
        given cancer_terms from the labkey table:
            - cancer_participant_disease
        Note: this is not an extensive search / confirmation of the various
        tumours our participants may have gotten. But rather returns
        participants with disease_types reported on enrolment to the 100K genomes
        programme.

        Returns:
            pd.DataFrame: a dataframe with participant_ids and cancer
            disease types.
        """
        self._ensure_axis("cancer_terms", "cterms")
        cpd = self.loader.grab_data("cancer_participant_disease_by_types", cancer_types=self.cterms)
        ca = self.loader.grab_data("cancer_analysis_by_disease_type", cancer_types=self.cterms)
        cpd = self._ensure_pid_str(cpd)
        ca  = self._ensure_pid_str(ca)
        disease_table = pd.concat([cpd, ca], axis=0, ignore_index=True).drop_duplicates()
        self.cterm_table = disease_table
        self.pids["cterm"] = disease_table["participant_id"]


    def get_cabbr_pids(self):
        """Query cancer analysis for participants with a given study abbreviation
        This is currently limited to the GEL 100K cancer cohort.
        """
        self._ensure_axis("cancer_abbr", "cabbr")
        df = self.loader.grab_data("cancer_analysis_by_abbr", abbrs=self.cabbr)
        df = self._ensure_pid_str(df)
        self.cabbr_table = df
        self.pids["cabbr"] = df["participant_id"]


    def _extract_icd10_code(self, df, text_col):
        """helper function to extract and clean ICD10 codes, which often come
        with additional '.' or are appended together in diag_all for HES data."""
        if df is None or df.empty:
            return pd.DataFrame(columns=["participant_id", "code"])
        tmp = df.copy()
        tmp[text_col] = (
            tmp[text_col]
            .astype(str)
            .str.replace(
                ".", 
                "", 
                regex=False
            ))
        tmp["code"] = tmp[text_col].str.extract(r"([A-Z][0-9]+)")
        out = tmp.loc[
            ~tmp["code"].isna(), 
            ["participant_id", "code"]
            ].drop_duplicates()
        
        return out.reset_index(drop=True)
    
    def _and_participants_clause(self):
        """Helper clause to adjust sql queries so they may be limited to participants
        already in the cohort or keep them participant agnostic. 
        
        Return 'AND participant_id IN (...)' or '' if no cohort PIDs yet."""
        try:
            pid_df = self.concat_cohort(self.pids)  # safe for dicts of Series/DFs
        except Exception:
            pid_df = pd.DataFrame(columns=["participant_id"])
        pids = (
            pid_df["participant_id"].astype(str).dropna().drop_duplicates()
            if "participant_id" in pid_df.columns else pd.Series([], dtype=str)
        )
        if pids.empty:
            return ""
        return "AND participant_id IN " + self.loader._format_in_clause(pids)
    

    def get_icd10_pids(self, limit=['all']):
        """retrieve participant_ids and cleaned icd-10 codes associated with the 
        given icd-10 codes from the tables defined in the query module.
        currently these
        hes:
            - hes-apc
            - hes-ae
            - hes-op
        mort:
            - mortality
        mental_health:
            - mhsds_medical_history_previous_diagnosis
            - mhsds_provisional_diagnosis
            - mhsds_primary_diagnosis
            - mhsds_secondary_diagnosis
            - mhsds_care_activity
            - mhsds_indirect_activity
            - mhmd_v4_event
            - mhldds_event
        cancer:
            - cancer_invest_sample_pathology
            - cancer_participant_tumour
            - cancer_registry
            - rtds
            - sact
            - av_tumour
        
        the sources queried can be limited using the limit argument.
        'all', 'hes', 'mort', 'mental_health', 'cancer'.

        Returns:
            pd.DataFrame: a dataframe with participant_ids and clean
            icd-10 codes.
        """

        valid = {"all", "hes", "mort", "mental_health", "cancer"}
        req = set(limit)
        if not req.issubset(valid):
            raise ValueError("limit must be a subset of {'all','hes','mort','mental_health','cancer'}")
        cats = ["hes", "mort", "mental_health", "cancer"] if "all" in req else list(req)

        # ensure codes present from featdict
        self._ensure_axis("icd10", "icd10s")

        from gelpack.queries.labkey_queries_v19 import ICD10_SOURCES

        frames = []
        for cat in cats:
            for spec in ICD10_SOURCES.get(cat, []):
                key = spec["key"]
                extra = {}
                if "table" in spec:
                    extra["table"] = spec["table"]
                if "code_col" in spec:
                    extra["code_col"] = spec["code_col"]

                # this loading is participant agnostic.
                df = self.loader.grab_data(
                    key,
                    likes={"like_diag": ("diag", self.icd10s)},
                    and_participants="",   
                    **extra,
                )
                df = self._ensure_pid_str(df)
                if df is None or df.empty:
                    continue

                # keep only canonical columns if present
                cols = ["participant_id", "diag"]
                if "event_date" in df.columns:
                    cols.append("event_date")
                out = df.loc[:, [c for c in cols if c in df.columns]].copy()

                # attach source metadata without schema coupling
                out["source"] = cat
                out["source_table"] = spec.get("table", key)
                frames.append(out)

        if not frames:
            warnings.warn("No participants found for these ICD-10 codes in selected sources.")
            return

        out = (
            pd.concat(frames, ignore_index=True)
            .assign(code=lambda d: d["diag"].astype(str).str.extract(r"([A-Z][0-9]+)")[0])
            .dropna(subset=["code"])
            .drop_duplicates(subset=["participant_id", "code", "source", "source_table"])
            .reset_index(drop=True)
        )

        self.icd10_table = out
        self.pids["icd10"] = (
            out["participant_id"].dropna().astype(str).drop_duplicates().reset_index(drop=True)
        )

        self.provenance.append({
            "action": "get_icd10_pids",
            "n_rows": int(len(out)),
            "n_unique_participants": int(out["participant_id"].nunique()),
            "categories": cats,
            "global_search": True,
        })


    ## ---------------------------------------------------------------##
    ##################### adding cohort features #######################
    ## ---------------------------------------------------------------##
    def age(self):
        """Get the age at consent and current age for each participant in the 
        self. This function does go over each source (cancer, icd10, hpo...)
        individually, which may lead to a few duplicate entries being calculated.
        But allows for comparison of the source cohorts downstream.

        Returns:
            age_table dictionary: The keys of which correspond to the keys of the sources,
            for each key a pd.dataframe with 'date_of_consent', 'participant_id',
            'year_of_birth', 'current_age' and 'age_at_consent'. In addition a
            QC flag has been added - this flags instances where the year_of_birth
            was set at a default value, the date of consent occured before 
            the start of the 100,000 Genomes Programme, and if participants
            were over 100 years old at the date of consent.
        """
        #### Age ####
        self.age_table = {}
        prov_rows = []

        for key, pid in self.pids.items():
            df = self.loader.grab_data("participant_age", participants=pid)
            df = self._ensure_pid_str(df)

            if df is None:
                df = pd.DataFrame(columns=[
                    "participant_id", "year_of_birth", "date_of_consent",
                    "current_age", "age_at_consent", "age_qc"
                ])

            df = df.copy()
            df["date_of_consent"] = pd.to_datetime(df.get("date_of_consent"), errors="coerce")
            df["year_of_birth"] = pd.to_datetime(df.get("year_of_birth"), format="%Y", errors="coerce")

            def _age_years(born, ref=None):
                if pd.isna(born):
                    return np.nan
                ref = ref or pd.Timestamp.today()
                return ref.year - born.year - ((ref.month, ref.day) < (born.month, born.day))

            df["current_age"] = df["year_of_birth"].apply(lambda b: _age_years(b))
            df["age_at_consent"] = df.apply(
                lambda x: _age_years(x["year_of_birth"], x["date_of_consent"]), axis=1
            )

            df["age_qc"] = ~(
                (df["year_of_birth"] == pd.Timestamp("1900-01-01"))
                | (df["date_of_consent"] < pd.Timestamp("2000-01-01"))
                | (df["age_at_consent"] > 100)
            )

            self.age_table[key] = df

            n_requested = int(len(pid)) if hasattr(pid, "__len__") else 0
            n_rows = int(len(df))
            n_parts = (
                int(df["participant_id"].nunique())
                if "participant_id" in df.columns else pd.NA
            )

            prov_rows.append({
                "group": key,
                "query_key": "participant_age",
                "n_requested": n_requested,
                "n_rows_returned": n_rows,
                "n_participants_returned": n_parts,
                "columns_returned": list(df.columns),
            })

        self.feature_tables["age"] = self.age_table
        self.feature_tables["age_provenance"] = pd.DataFrame(prov_rows)

        self.provenance.append({
            "action": "age",
            "groups": list(self.pids.keys()),
            "query_key": "participant_age",
        })


    def ancestry(self, threshold=0.8):
        """Get the predicted Ancestry for each participant in the 
        self. like age(), this function sets a threshold of 0.8, above which a
         caclucated ancestry will be assigned to the participant id.

        Returns:
            ancestry_table (dictionary): The keys of which correspond to the keys of 
            the sources, for each key a pd.dataframe with 'participant_id' 
            and 'predicted_ancestry'. Ancestries are set by the highest scoring
            predicted ancestry in aggregate_gvcf_sample_stats. If no single 
            score was higher than 0.8 the ancestry is set to unassigned (UNA). 
            If the participant is not part of this table the ancestry is set 
            to unknown (UNO).
        """
        ## ancestry  ##
        self.ancestry_table = {}
        self.ancestry_best = {}
        prov_rows = []

        CANON = {"AFR", "SAS", "EAS", "EUR", "AMR"}

        for key, pid in self.pids.items():
            df = self.loader.grab_data("aggregate_gvcf_ancestry", participants=pid)
            df = self._ensure_pid_str(df)

            if df is None or df.empty:
                out = pd.DataFrame({"participant_id": list(pid), "predicted_ancestry": "UNO"})
                self.ancestry_table[key] = out
                prov_rows.append({
                    "group": key, "strategy": "no_data", "threshold": threshold,
                    "score_cols": "", "n_total": int(len(pid)),
                    "n_assigned": 0, "n_UNA": 0, "n_UNO": int(len(pid)),
                })
                continue

            df = df.copy()

            # Prefer canonical alias columns if present; else, any numeric score columns.
            alias_cols = [c for c in df.columns if c in CANON]
            strategy = "canonical" if alias_cols else "numeric_infer"
            score_cols = alias_cols or [
                c for c in df.columns
                if c != "participant_id" and pd.api.types.is_numeric_dtype(df[c])
            ]


            if not score_cols:
                out = pd.DataFrame({"participant_id": list(pid), "predicted_ancestry": "UNO"})
                self.ancestry_table[key] = out
                prov_rows.append({
                    "group": key, "strategy": "no_numeric_cols", "threshold": threshold,
                    "score_cols": "", "n_total": int(len(pid)),
                    "n_assigned": 0, "n_UNA": 0, "n_UNO": int(len(pid)),
                })
                continue

            long = df.melt(
                id_vars="participant_id",
                value_vars=score_cols,
                var_name="axis",
                value_name="score",
            )
            long["score"] = pd.to_numeric(long["score"], errors="coerce")

            best = (
                long.sort_values(["participant_id", "score"], 
                    ascending=[True, False])
                .drop_duplicates("participant_id")
                .rename(columns={"axis": "best_axis", "score": "best_score"})
            )

            best["predicted_ancestry"] = best.apply(
                lambda r: r["best_axis"] if pd.notna(r["best_score"]) and float(r["best_score"]) >= 0.8 else "UNA",
                axis=1,
            )
            assigned = best.copy()
            assigned["predicted_ancestry"] = assigned.apply(
                lambda r: r["best_axis"] if pd.notna(r["best_score"]) and float(r["best_score"]) >= threshold else "UNA",
                axis=1,
            )
            out = assigned[["participant_id", "predicted_ancestry"]]

            # Participants with missing ancestry data => UNO
            have = set(out["participant_id"])
            miss = [p for p in list(pid) if p not in have]
            if miss:
                out = pd.concat(
                    [out, pd.DataFrame({"participant_id": miss, "predicted_ancestry": "UNO"})],
                    ignore_index=True,
                )

            # Save artifacts
            self.ancestry_table[key] = out
            self.ancestry_best[key] = assigned[["participant_id", "best_axis", "best_score"]]

            n_assigned = int((assigned["best_score"] >= threshold).sum())
            n_UNA = int((assigned["best_score"] < threshold).sum())
            n_UNO = int(len(miss))
            prov_rows.append({
                "group": key,
                "strategy": strategy,
                "threshold": threshold,
                "score_cols": ",".join(score_cols),
                "n_total": int(len(pid)),
                "n_assigned": n_assigned,
                "n_UNA": n_UNA,
                "n_UNO": n_UNO,
            })

        # Expose in feature_tables for easy access
        self.feature_tables["ancestry"] = self.ancestry_table
        self.feature_tables["ancestry_best"] = self.ancestry_best
        self.feature_tables["ancestry_provenance"] = pd.DataFrame(prov_rows)

        # High-level breadcrumb
        self.provenance.append({
            "action": "ancestry",
            "groups": list(self.pids.keys()),
            "threshold": threshold,
            "n_groups": len(self.pids),
        })


    def sex(self):
        """query the labkey tables for phenotypic sex per participant_id.
        This data is the participant's stated sex by the clinician at the GMC.
        """
        ## Sex ##
        self.sex_table = {}
        prov_rows = []

        for key, pid in self.pids.items():
            df = self.loader.grab_data("participant_sex", participants=pid)
            df = self._ensure_pid_str(df) 

            if df is None:
                df = pd.DataFrame(columns=["participant_id","sex"])
            self.sex_table[key] = df

            n_requested = int(len(pid)) if hasattr(pid, "__len__") else 0
            n_rows = int(len(df)) if df is not None else 0
            n_parts = (
                int(df["participant_id"].nunique())
                if df is not None and "participant_id" in df.columns
                else pd.NA
            )
            cols = list(df.columns) if df is not None else []

            prov_rows.append({
                "group": key,
                "query_key": "participant_sex",
                "n_requested": n_requested,
                "n_rows_returned": n_rows,
                "n_participants_returned": n_parts,
                "columns_returned": cols,
            })

        self.feature_tables["sex"] = self.sex_table
        self.feature_tables["sex_provenance"] = pd.DataFrame(prov_rows)

        self.provenance.append({
            "action": "sex",
            "groups": list(self.pids.keys()),
            "query_key": "participant_sex",
        })


    def mortality(self, pick='earliest'):
        """Extract the death date from all registered mortality sources
            for a set of participant_id.
        Expects queries named 'mortality_source_*' that each return:
            - participant_id
            - dod  (date-of-death; canonical alias in SQL)

        Args:
            pids (list): list of participant ids to include
            version (str): Data release version

        Returns:
            mortality_table (pd.DataFrame): date of death per participant_id
        """
        if not self.pids:
            raise ValueError("mortality() requires participants.")

        # Fixed list of sources; SQLs must alias date to 'dod'
        sources = [
            "death_details_by_pids",
            "mortality_by_pids",
            "rare_diseases_pedigree_member_dod",
        ]

        self.mortality_table = {}
        events_prov = []

        for group, pid in self.pids.items():
            frames = []
            for key in sources:
                df = self.loader.grab_data(key, participants=pid)
                df = self._ensure_pid_str(df)
                if df is not None and not df.empty:
                    df = df.loc[:, ["participant_id", "dod"]].copy()
                    frames.append(df)

            if not frames:
                # Still return a well-formed frame for the group
                out = pd.DataFrame({"participant_id": pd.Series(pid).astype(str).drop_duplicates()})
                out["date_of_death"] = pd.NaT
                out["status"] = "Alive"
                self.mortality_table[group] = out
                events_prov.append({"group": group, "n_events": 0})
                continue

            events = (pd.concat(frames, ignore_index=True)
                        .dropna(subset=["dod"])
                        .assign(dod=lambda d: pd.to_datetime(d["dod"], errors="coerce"))
                        .dropna(subset=["dod"])
                        .drop_duplicates())

            # pick per participant
            if pick == "latest":
                picked = (events.sort_values(["participant_id", "dod"])
                                .groupby("participant_id", as_index=False)
                                .tail(1))
            else:
                picked = (events.sort_values(["participant_id", "dod"])
                                .groupby("participant_id", as_index=False)
                                .head(1))
            picked = picked.rename(columns={"dod": "date_of_death"})[["participant_id", "date_of_death"]]

            # ensure all cohort pids appear
            out = (pd.DataFrame({"participant_id": pd.Series(pid, dtype="string").drop_duplicates()})
                    .merge(picked, on="participant_id", how="left"))
            out["status"] = np.where(out["date_of_death"].isna(), "Alive", "Deceased")

            self.mortality_table[group] = out
            events_prov.append({"group": group, "n_events": int(len(events))})

        # register feature
        self.feature_tables["mortality"] = self.mortality_table
        self.feature_tables["mortality_provenance"] = pd.DataFrame(events_prov)
        self.provenance.append({"action": "mortality", "pick": pick, "groups": list(self.pids.keys())})


    ## ---------------------------------------------------------------##
    ############## annotate cohort with ICD10 codes ####################
    ## ---------------------------------------------------------------##

    def tag_icd10(self, codes, *, tag_name, limit=['all']):
        """
        Annotate the existing cohort members with ICD-10 presence flags,
        counts, and earliest event date (where available). also captures a source
        for the data/inclusion.

        Does NOT modify self.pids. Stores:
        - self.feature_tables['icd10_flags'][label]
        - self.icd10_tags[tag_name] 
        """
        if not self.pids:
            raise ValueError("tag_icd10() requires an existing cohort (participants).")

        valid = {"all", "hes", "mort", "mental_health", "cancer"}
        req = set(limit)
        if not req.issubset(valid):
            raise ValueError("limit must be a subset of {'all','hes','mort','mental_health','cancer'}")
        cats = ["hes", "mort", "mental_health", "cancer"] if "all" in req else list(req)

        # normalise codes
        terms = list(dict.fromkeys(str(x) for x in (codes or [])))
        if not terms:
            raise ValueError("Provide at least one ICD-10 code or prefix.")
        if not tag_name:
            raise ValueError("Provide tag_name=... for the boolean flag column.")

        # cohort PIDs (as strings)
        base_pids = self.concat_cohort(self.pids)
        if base_pids.empty or "participant_id" not in base_pids.columns:
            raise ValueError("No participant_ids available in the cohort.")
        pid_list = base_pids["participant_id"].astype(str).dropna().drop_duplicates()

        from gelpack.queries.labkey_queries_v19 import ICD10_SOURCES

        frames = []
        for cat in cats:
            for spec in ICD10_SOURCES.get(cat, []):
                key = spec["key"]
                extra = {}
                if "table" in spec:
                    extra["table"] = spec["table"]
                if "code_col" in spec:
                    extra["code_col"] = spec["code_col"]

                df = self.loader.grab_data(
                    key,
                    participants=pid_list,
                    likes={"like_diag": ("diag", terms)},
                    and_participants=" AND participant_id IN {participants}",  # <— scoped to cohort
                    **extra,
                )
                df = self._ensure_pid_str(df)
                if df is None or df.empty:
                    continue

                cols = ["participant_id", "diag"]
                if "event_date" in df.columns:
                    cols.append("event_date")
                out = df.loc[:, [c for c in cols if c in df.columns]].copy()

                out["code"] = out["diag"].astype(str).str.extract(r"([A-Z][0-9]+)")[0]
                out = out.dropna(subset=["code"])
                out["source"] = cat
                out["source_table"] = spec.get("table", key)
                frames.append(out)

        # Long raw matches (may be empty if no hits)
        tag_df = (pd.concat(frames, ignore_index=True) if frames else
                pd.DataFrame(columns=["participant_id", "code", "source", "source_table", "event_date"]))

        # Store raw
        if not hasattr(self, "icd10_tags"):
            self.icd10_tags = {}
        self.icd10_tags[tag_name] = tag_df

        # Build Boolean per-participant feature over the cohort
        pid_base = pd.DataFrame({"participant_id": pid_list})
        hits = tag_df[["participant_id"]].drop_duplicates() if not tag_df.empty else pd.DataFrame(columns=["participant_id"])
        hits["__hit__"] = True
        flags = pid_base.merge(hits, on="participant_id", how="left")
        colname = f"icd10_{tag_name}"
        flags[colname] = flags["__hit__"].fillna(False).astype(bool)
        flags = flags[["participant_id", colname]]

        # Expose via feature_tables so concat_all() picks it up
        ft = self.feature_tables.setdefault("icd10_flags", {})
        ft[tag_name] = flags

        self.provenance.append({
            "action": "tag_icd10",
            "tag_name": tag_name,
            "n_rows": int(len(tag_df)),
            "n_participants_tagged": int(flags[colname].sum()),
            "categories": cats,
            "cohort_scoped": True,
        })

    ## ---------------------------------------------------------------##
    ##################### sample level data ############################
    ## ---------------------------------------------------------------##
    def rd_sample_data(self):
        self.rd_samples = {}
        prov_rows = []

        seeded = self.input_platekeys is not None and len(self.input_platekeys) > 0

        if seeded:
            # Return only the explicitly provided platekeys
            df = self.loader.grab_data("rd_samples_by_platekeys", platekeys=self.input_platekeys)
            df = self._ensure_pid_str(df)
            for key, _pid in self.pids.items():
                self.rd_samples[key] = df
                prov_rows.append({
                    "action": "rd_sample_data",
                    "group": key,
                    "seed": "platekeys",
                    "queries_used": ["rd_samples_by_platekeys"],
                    "n_rows": 0 if df is None else int(len(df)),
                    "n_participants_in": None,
                    "n_platekeys_in": int(self.input_platekeys.nunique()),
                    "columns_out": [] if df is None else list(df.columns),
                })
        else:
            # No input platekeys: gather all RD samples for the cohort participants
            for key, pid in self.pids.items():
                df = self.loader.grab_data("rd_samples_by_participants", participants=pid)
                df = self._ensure_pid_str(df)
                self.rd_samples[key] = df
                prov_rows.append({
                    "action": "rd_sample_data",
                    "group": key,
                    "seed": "participants",
                    "queries_used": ["rd_samples_by_participants"],
                    "n_rows": 0 if df is None else int(len(df)),
                    "n_participants_in": int(len(pid)),
                    "n_platekeys_in": None,
                    "columns_out": [] if df is None else list(df.columns),
                })

            # Only infer platekeys from results if we were NOT seeded by platekeys
            try:
                all_rd = self.concat_cohort(self.rd_samples)
                if "sample_platekey" in all_rd.columns:
                    self.platekeys = (
                        all_rd["sample_platekey"].dropna().astype(str).drop_duplicates().reset_index(drop=True)
                    )
            except Exception:
                pass

        self.sample_tables["rare_disease_samples"] = self.rd_samples
        self.feature_tables["rd_samples_provenance"] = pd.DataFrame(prov_rows)
        self.provenance.append({
            "action": "rd_sample_data",
            "seed": "platekeys" if seeded else "participants",
            "groups": list(self.pids.keys()),
        })


    def ca_sample_data(self):
        """query labkey for cancer sample information. at default, data is retrieved from 
        cancer analysis (for nsv4 sample stats) and supplemented with dragen 3.2
        realigned samples. But this can be changed by updating the QUERIES file.
        We make a distinction between cohorts generated from platekeys and
        cohorts generated from participant_ids / disease terms. As we want to limit
        the output to specific platekeys in the former, but may want to grab all 
        samples per participant_id for the latter.
        (there are multiple tumour samples available)

        """
        self.cancer_samples = {}
        prov_rows = []

        seeded = self.input_platekeys is not None and len(self.input_platekeys) > 0

        if seeded:
            # Return only the explicitly provided platekeys
            df = self.loader.grab_data("cancer_samples_by_platekeys", platekeys=self.input_platekeys)
            df = self._ensure_pid_str(df)
            for key, _pid in self.pids.items():
                self.cancer_samples[key] = df
                prov_rows.append({
                    "action": "ca_sample_data",
                    "group": key,
                    "seed": "platekeys",
                    "queries_used": ["cancer_samples_by_platekeys"],
                    "n_rows": 0 if df is None else int(len(df)),
                    "n_participants_in": None,
                    "n_platekeys_in": int(self.input_platekeys.nunique()),
                    "columns_out": [] if df is None else list(df.columns),
                })
        else:
            # No input platekeys: gather all CA samples for the cohort participants
            for key, pid in self.pids.items():
                df = self.loader.grab_data("cancer_samples_by_participants", participants=pid)
                df = self._ensure_pid_str(df)
                self.cancer_samples[key] = df
                prov_rows.append({
                    "action": "ca_sample_data",
                    "group": key,
                    "seed": "participants",
                    "queries_used": ["cancer_samples_by_participants"],
                    "n_rows": 0 if df is None else int(len(df)),
                    "n_participants_in": int(len(pid)),
                    "n_platekeys_in": None,
                    "columns_out": [] if df is None else list(df.columns),
                })

            # Infer platekeys ONLY from tumour_sample_platekey (PRIMARY)
            all_ca = self.concat_cohort(self.cancer_samples)
            if isinstance(all_ca, pd.DataFrame) and not all_ca.empty and "tumour_sample_platekey" in all_ca.columns:
                inferred = (all_ca["tumour_sample_platekey"]
                            .dropna().astype(str).drop_duplicates().reset_index(drop=True))
                if not inferred.empty:
                    self.platekeys = inferred  # PRIMARY = tumour
            # (No fallback to germline here by design.)

        self.sample_tables["cancer_samples"] = self.cancer_samples
        self.feature_tables["cancer_samples_provenance"] = pd.DataFrame(prov_rows)
        self.provenance.append({
            "action": "ca_sample_data",
            "seed": "platekeys" if seeded else "participants",
            "groups": list(self.pids.keys()),
        })

    def load_gmc_registration(self):
        """get each pids site of registration through a simple lookup of
        the participant table."""

        self.gmc_registration = {}
        prov_rows = []

        for key, pid in self.pids.items():
            df = self.loader.grab_data("participant_gmc_registration", participants=pid)
            df = self._ensure_pid_str(df)
            self.gmc_registration[key] = df

            prov_rows.append({
                "action": "gmc_registration",
                "group": key,
                "seed": "participants",
                "queries_used": ["participant_gmc_registration"],
                "n_rows": 0 if df is None else int(len(df)),
                "n_participants_in": int(len(pid)),
                "columns_out": [] if df is None else list(df.columns),
            })

        # expose in feature tables
        self.feature_tables["gmc_registration"] = self.gmc_registration

        self.feature_tables["gmc_registration_provenance"] = pd.DataFrame(prov_rows)
        self.provenance.append({
            "action": "gmc_registration",
            "groups": list(self.pids.keys()),
        })


    def omics_sample_metadata(self):
        """Extract sample metadata from labkey. Specifically looking at omics
        sample availability and location. Please keep in mind the sample counts
        are subject to change.
        
        Returns: 
            omics_sample_data (dictionary): The keys of which correspond to the
            different cohorts in self.pids. This table contains info on each sample
            collected.
            omics_sample_counts: the number of samples per sample_type.
        """
        self.omics_sample_data = {}
        self.omics_sample_counts = {}
        self.omics_sample_location = {}
        prov_rows = []

        for key, pid in self.pids.items():
            df = self.loader.grab_data("omics_metadata_by_participants", participants=pid)
            df = self._ensure_pid_str(df)
            self.omics_sample_data[key] = df

            # Counts: histogram of aliquots per sample_type
            if df is None or df.empty:
                counts = pd.DataFrame(columns=["sample_type", "aliquots", "count"])
                location = pd.DataFrame(
                    columns=["sample_type", "laboratory_sample_gmc_trust", "laboratory_sample_gmc_ods_code", "count"]
                )
            else:
                counts = (
                    df.groupby(["sample_type", "aliquots"])
                    .size()
                    .reset_index(name="count")
                )
                location = (
                    df.groupby(["sample_type", "laboratory_sample_gmc_trust", "laboratory_sample_gmc_ods_code"])
                    .size()
                    .reset_index(name="count")
                )

            self.omics_sample_counts[key] = counts
            self.omics_sample_location[key] = location

            prov_rows.append({
                "action": "omics_sample_metadata",
                "group": key,
                "seed": "participants",
                "queries_used": ["omics_metadata_by_participants"],
                "n_rows": 0 if df is None else int(len(df)),
                "n_participants_in": int(len(pid)),
                "columns_out": [] if df is None else list(df.columns),
            })

        # expose main sample-level table in sample_tables
        self.sample_tables["omics_sample_data"] = self.omics_sample_data
        # optional: expose summaries as features for convenience
        self.feature_tables["omics_sample_counts"] = self.omics_sample_counts
        self.feature_tables["omics_sample_location"] = self.omics_sample_location

        # per-call provenance table + overall trail
        self.feature_tables["omics_provenance"] = pd.DataFrame(prov_rows)
        self.provenance.append({
            "action": "omics_sample_metadata",
            "groups": list(self.pids.keys()),
        })

        
    ## ---------------------------------------------------------------##
    ############ aggregating and summarysing the cohort ################
    ## ---------------------------------------------------------------##
    def concat_cohort(self, dc):
        """concatenate the dataframes in a dictionary, note; the dicitonarys
        should have the same structure / columns. 

        Args:
            dc (dictionary): Dictionary value's are pd.DataFrames, one per Key.

        Returns:
            pd.DataFrame: A concatenated pandas dataframe with unique values in
            the dictionary.
        """
        if not isinstance(dc, dict) or not dc:
            return pd.DataFrame()

        vals = list(dc.values())

        # all DataFrames
        if all(isinstance(v, pd.DataFrame) for v in vals):
            frames = []
            for v in vals:
                df = v.copy()
                if "participant_id" in df.columns:
                    df["participant_id"] = df["participant_id"].astype("string")
                frames.append(df)
            return (
                pd.concat(frames, ignore_index=True)
                .drop_duplicates()
                .reset_index(drop=True)
            )

        # all sequences: build a DF with participant_id
        if all(isinstance(v, (pd.Series, list, tuple)) for v in vals):
            frames = []
            for v in vals:
                s = v if isinstance(v, pd.Series) else pd.Series(list(v))
                s = s.dropna().astype("string").rename("participant_id")
                frames.append(pd.DataFrame({"participant_id": s}))
            return (
                pd.concat(frames, ignore_index=True)
                .drop_duplicates()
                .reset_index(drop=True)
            )

        # mixed types: coerce sequences to DF, pass DFs through (making sure pid = str)
        frames = []
        for v in vals:
            if isinstance(v, pd.DataFrame):
                df = v.copy()
                if "participant_id" in df.columns:
                    df["participant_id"] = df["participant_id"].astype("string")
                frames.append(df)
            elif isinstance(v, (pd.Series, list, tuple)):
                s = v if isinstance(v, pd.Series) else pd.Series(list(v))
                s = s.dropna().astype("string").rename("participant_id")
                frames.append(pd.DataFrame({"participant_id": s}))
            else:
                raise TypeError("concat_cohort: unsupported dict value type.")

        return (
            pd.concat(frames, ignore_index=True)
            .drop_duplicates()
            .reset_index(drop=True)
        )


    def concat_all(self):
        """concatenates all the features and sources in the cohort to a long
        type table, matching each patient/diagnosis combination with its features
        this is usefull for filtering, especially when you want to exclude
        participants of a certain feature (e.g. age) and diagnosis - 
        but don't want to blanket remove based on features. (e.g. 
        filter participants age 30 and under for diagnosis A, but include all ages
        for diagnosis B: a participant who is 25 yo and was diagnosed for both
        A and B will be included in the cohort)

        Raises:
            RuntimeError: When no valid cohort is present
        """


        diag_keys = ["icd10", "dterm", "cterm", "cabbr", "hpo", "custom"]
        
        # only return a runtime error if there are no participant ids:
        if all(k not in self.pids for k in diag_keys):
            raise RuntimeError("No valid cohorts to concatenate found")

        # Union of current cohort membership (normalize only the key column)
        pid_union = self.concat_cohort(self.pids)
        if pid_union.empty or "participant_id" not in pid_union.columns:
            raise RuntimeError("No participants available to concatenate.")
        pid_union = pid_union.copy()
        pid_union["participant_id"] = pid_union["participant_id"].astype("string")
        self.all_pids = pid_union["participant_id"]

        # Build diagnosis table
        conc_tables = []
        for key in diag_keys:
            if key not in self.pids:
                continue
            if key == "dterm" and hasattr(self, "dterm_table"):
                t = self.dterm_table.rename(columns={"normalised_specific_disease": "diag"})
                t = t.loc[:, ["participant_id", "diag"]].copy()
                t["participant_id"] = t["participant_id"].astype("string")
                conc_tables.append(t)
            elif key == "icd10" and hasattr(self, "icd10_table"):
                t = self.icd10_table.rename(columns={"code": "diag"})
                t = t.loc[:, ["participant_id", "diag"]].copy()
                t["participant_id"] = t["participant_id"].astype("string")
                conc_tables.append(t)
            elif key == "cterm" and hasattr(self, "cterm_table"):
                t = self.cterm_table.rename(columns={"disease_type": "diag"})
                t = t.loc[:, ["participant_id", "diag"]].copy()
                t["participant_id"] = t["participant_id"].astype("string")
                conc_tables.append(t)
            elif key == "cabbr" and hasattr(self, "cabbr_table"):
                t = self.cabbr_table.rename(columns={"study_abbreviation": "diag"})
                t = t.loc[:, ["participant_id", "diag"]].copy()
                t["participant_id"] = t["participant_id"].astype("string")
                conc_tables.append(t)
            elif key == "hpo" and hasattr(self, "hpo_table"):
                t = self.hpo_table.rename(columns={"normalised_hpo_id": "diag"})
                t = t.loc[:, ["participant_id", "diag"]].copy()
                t["participant_id"] = t["participant_id"].astype("string")
                conc_tables.append(t)
            elif key == "custom":
                t = pd.DataFrame({"participant_id": self.pids["custom"]})
                t["participant_id"] = t["participant_id"].astype("string")
                t["diag"] = "custom"
                conc_tables.append(t.loc[:, ["participant_id", "diag"]])

        if not conc_tables:
            raise RuntimeError("No diagnosis/ontology tables to anchor concat_all().")

        diag_table = (
            pd.concat(conc_tables, ignore_index=True)
            .dropna(subset=["participant_id", "diag"])
            .drop_duplicates()
            .reset_index(drop=True)
        )

        # merge all features & samples on participant_id (keep original dtypes elsewhere)
        tmp_merge = [pid_union.loc[:, ["participant_id"]].copy()]

        # participant-level features
        for _, feature in self.feature_tables.items():
            if isinstance(feature, dict) and feature:
                merged = self.concat_cohort(feature)
                if not merged.empty and "participant_id" in merged.columns:
                    m = merged.copy()
                    m["participant_id"] = m["participant_id"].astype("string")
                    tmp_merge.append(m)

        # sample-level tables 
        self.all_cancer_samples = None
        self.all_rd_samples = None
        for name, samples in self.sample_tables.items():
            merged = self.concat_cohort(samples)
            if name == "cancer_samples":
                self.all_cancer_samples = merged
            elif name == "rare_disease_samples":
                self.all_rd_samples = merged
            if not merged.empty and "participant_id" in merged.columns:
                m = merged.copy()
                m["participant_id"] = m["participant_id"].astype("string")
                tmp_merge.append(m)

        # safe reduce (works even if there’s only the seed frame)
        if len(tmp_merge) == 1:
            df_merged = tmp_merge[0]
        else:
            df_merged = reduce(
                lambda left, right: pd.merge(left, right, on=["participant_id"], how="left"),
                tmp_merge,
            )

        diag_table = diag_table.copy()
        diag_table["participant_id"] = diag_table["participant_id"].astype("string")
        self.all_data = pd.merge(diag_table, df_merged, on="participant_id", how="left")



    def summarize(self, pid=None, anc=None, mort=None, age=None, sex=None, reg=None, omics=None):
        """gather up self data and return summary stats. including:
        - number of unique participants
        - mean age at consent
        - standard deviation of mean age at consent
        - mean current age
        - fraction female
        - fraction diseased
        - fraction of each ancestry.

        Args:
            pid (np.Array): array of participant ids.
            anc (pd.DataFrame): 2 column table with participant_id and 
                predicted ancestry.
            mort (pd.DatFrame): 2 column table with participant_id and 
                date_of_death
            age (pd.DataFrame): Table with age_at_conset and current_age 
                per participant id.
            sex (pd.DataFrame): 2 column table with a sex per participant_id.

        Returns:
            pd.DataFrame: single row DataFrame with summary stats.
        """
        # not all features have to be calculated/present.
        # only include those that are.
        # instead of building a dataframe in one go we'll concat each
        # column that we actually have.
        series = []
        if pid is not None:
            series.append(pd.Series(len(pid), name="n_unique_participants"))
        if age is not None and not age.empty:
            series += [
                pd.Series(np.mean(age["age_at_consent"]), name="mean_consent_age"),
                pd.Series(np.std(age["age_at_consent"]),  name="std_consent_age"),
                pd.Series(np.mean(age["current_age"]),     name="mean_current_age"),
                pd.Series(np.std(age["current_age"]),      name="std_current_age"),
            ]
        if sex is not None and not sex.empty:
            vc = sex["sex"].value_counts(normalize=True)
            series.append(pd.Series(vc.get("Female", np.nan), name="fraction_female"))
        if mort is not None and not mort.empty:
            vc = mort["status"].value_counts(normalize=True)
            series.append(pd.Series(vc.get("Deceased", np.nan), name="fraction_deceased"))

        # omics coverage (by participant_id presence)
        if omics is not None and not omics.empty and pid is not None:
            n_with_omics = int(omics["participant_id"].nunique())
            series += [
                pd.Series(n_with_omics, name="n_participants_with_omics"),
                pd.Series(n_with_omics / len(pid) if len(pid) else np.nan, name="fraction_with_omics"),
                pd.Series(int(len(omics)), name="n_total_omics_rows"),
            ]

        # GMC registration coverage
        if reg is not None and not reg.empty and pid is not None:
            n_with_reg = int(reg["participant_id"].nunique())
            series += [
                pd.Series(n_with_reg, name="n_with_gmc_registration"),
                pd.Series(n_with_reg / len(pid) if len(pid) else np.nan, name="fraction_with_gmc_registration"),
            ]

        summary = pd.concat(series, axis=1) if series else pd.DataFrame(index=[0])

        if anc is not None and not anc.empty:
            anc_dist = (
                anc.predicted_ancestry.value_counts(normalize=True).rename_axis("ancestry").reset_index(name="prop")
            )
            anc_wide = anc_dist.set_index("ancestry").T
            summary = pd.concat([summary.reset_index(drop=True), anc_wide.reset_index(drop=True)], axis=1)

        return summary
    

    def summary_stats(self):
        """calculate feature stats of the cohort, this should ideally be all
        features (age, sex, ancestry, mortality). depends on the summarize() 
        function. First collapse all the participant ids in the featuers (as they
        are generated per diag/ontology, e.g. self.age has got a table for
        icd10, hpo part of the cohort.)

        Returns:
            self.summary: a summary table of the entire cohort. not per ontology.
        """

        # setting tables to none to allow summary() to only use those which we
        # have calculated.
        all_ancestry = all_mortality = all_age = all_sex = all_omics = all_reg = None
        # collapse all the features we have identified
        all_pids = self.concat_cohort(self.pids)["participant_id"]

        all_age       = self.concat_cohort(self.feature_tables.get("age", {}))
        all_ancestry  = self.concat_cohort(self.feature_tables.get("ancestry", {}))
        all_sex       = self.concat_cohort(self.feature_tables.get("sex", {}))
        all_mortality = self.concat_cohort(self.feature_tables.get("mortality", {}))
        all_reg       = self.concat_cohort(self.feature_tables.get("gmc_registration", {}))

        # omics lives in sample_tables:
        all_omics     = self.concat_cohort(self.sample_tables.get("omics_sample_data", {}))

        self.summary = self.summarize(
            pid=all_pids,
            anc=all_ancestry,
            mort=all_mortality,
            age=all_age,
            sex=all_sex,
            reg=all_reg,
            omics=all_omics,
        )


    def build_normalized_views(self, participant_strategy="first"):
        """
        Create normalised dual views of participant and samples:

        - participant_view: one row per participant_id (reduced by `first`)
        - sample_view: one row per sample with:
            - sample_id      platekey (tumour/germline or RD plate_key)
            - sample_origin: {"cancer","rare_disease"}
            - sample_role:  {"tumour","germline"} where applicable
        """

        # participant-normalised view 
        if not hasattr(self, "all_data"):
            try:
                self.concat_all()
            except Exception:
                self.all_data = pd.DataFrame(columns=["participant_id"])

        if not self.all_data.empty and "participant_id" in self.all_data.columns:
            if participant_strategy != "first":
                pass  # only 'first' supported for now
            self.participant_view = (
                self.all_data.sort_values(["participant_id"])
                            .groupby("participant_id", as_index=False)
                            .first()
            )
        else:
            self.participant_view = pd.DataFrame(columns=["participant_id"])

        # sample-normalised view (CA tumour+germline, RD germline)
        frames, prov = [], []

        # cancer samples (expect tumour_sample_platekey & germline_sample_platekey)
        ca_dict = self.sample_tables.get("cancer_samples")
        if isinstance(ca_dict, dict):
            ca_all = self.concat_cohort(ca_dict)
            if ca_all is not None and not ca_all.empty:
                if "tumour_sample_platekey" in ca_all.columns:
                    t = ca_all.dropna(subset=["tumour_sample_platekey"]).copy()
                    if not t.empty:
                        t["sample_id"] = t["tumour_sample_platekey"]
                        t["sample_origin"] = "cancer"
                        t["sample_role"] = "tumour"
                        frames.append(t)
                        prov.append({"source": "cancer", "role": "tumour", "n_rows": int(len(t))})

                if "germline_sample_platekey" in ca_all.columns:
                    g = ca_all.dropna(subset=["germline_sample_platekey"]).copy()
                    if not g.empty:
                        g["sample_id"] = g["germline_sample_platekey"]
                        g["sample_origin"] = "cancer"
                        g["sample_role"] = "germline"
                        frames.append(g)
                        prov.append({"source": "cancer", "role": "germline", "n_rows": int(len(g))})

        # rare disease samples (expect plate_key; treat as germline)
        rd_dict = self.sample_tables.get("rare_disease_samples")
        if isinstance(rd_dict, dict):
            rd_all = self.concat_cohort(rd_dict)
            if rd_all is not None and not rd_all.empty and "sample_platekey" in rd_all.columns:
                rd = rd_all.dropna(subset=["sample_platekey"]).copy()
                if not rd.empty:
                    rd["sample_id"] = rd["sample_platekey"]
                    rd["sample_origin"] = "rare_disease"
                    rd["sample_role"] = "germline"
                    frames.append(rd)
                    prov.append({"source": "rare_disease", "role": "germline", "n_rows": int(len(rd))})


        if frames:
            all_cols = sorted(set(c for f in frames for c in f.columns))
            base = ["participant_id", "sample_id", "sample_origin", "sample_role"]
            ordered = base + [c for c in all_cols if c not in base]
            self.sample_view = pd.concat(frames, ignore_index=True)[ordered]
        else:
            self.sample_view = pd.DataFrame(columns=["participant_id", "sample_id", "sample_origin", "sample_role"])

        # provenance
        prov_df = pd.DataFrame(prov)
        self.feature_tables["normalized_views_provenance"] = prov_df
        self.provenance.append({
            "action": "build_normalized_views",
            "participant_strategy": participant_strategy,
            "n_participants_view": int(len(self.participant_view)),
            "n_samples_view": int(len(self.sample_view)),
            "sources": prov_df.to_dict("records") if not prov_df.empty else [],
        })

    def build_normalised_views(self, participant_strategy="first"):
        return self.build_normalized_views(participant_strategy=participant_strategy)

    def load_icd10_data(self):
        """gelpack contains a translation file for icd-10 codes, this function
        loads in the relevant file - used for generating summary statitics per
        icd-10 code.
        """
        from importlib import resources
        # Prefer the modern API when available; fall back for older Python.
        try:
            with resources.as_file(resources.files("gelpack") / "coding19.tsv") as p:
                df = pd.read_csv(p, sep="\t", usecols=["coding", "meaning"])
        except Exception:
            # Back-compat with resources.path if files/as_file is unavailable
            with resources.path("gelpack", "coding19.tsv") as p:
                df = pd.read_csv(p, sep="\t", usecols=["coding", "meaning"])

        self.icd10_lookup = df.rename(columns={"coding": "code"})

    # Provenance
        self.provenance.append({
            "action": "load_icd10_data",
            "source": "pkg://gelpack/coding19.tsv",
            "n_rows": int(len(self.icd10_lookup))
        })


    def ontology_stats(self):
        """calculate cohort summary per ontology/source (hpo, icd10, cterm,
        dterm, custom).

        Raises:
            RuntimeError: Returns an error when there are no ontologies to
            summarize.

        Returns:
            self.ont_vcount and self.icd10_overlap_matrix
        """
        # Build per-cohort (key) summary frames via existing self.summarize(...)
        summ_dict = {}
        for key, pid in self.pids.items():
            # Defensive fetches (these should exist if their respective methods were run)
            anc = getattr(self, "ancestry_table", {}).get(key, pd.DataFrame())
            mort = getattr(self, "mortality_table", {}).get(key, pd.DataFrame())
            age  = getattr(self, "age_table", {}).get(key, pd.DataFrame())
            sex  = getattr(self, "sex_table", {}).get(key, pd.DataFrame())

            s = self.summarize(
                pid=pid.drop_duplicates(),
                anc=anc.drop_duplicates() if not anc.empty else anc,
                mort=mort.drop_duplicates() if not mort.empty else mort,
                age=age.drop_duplicates() if not age.empty else age,
                sex=sex.drop_duplicates() if not sex.empty else sex,
            )
            s["source"] = key
            summ_dict[key] = s

        self.ont_summary = pd.concat(summ_dict.values(), ignore_index=True) if summ_dict else pd.DataFrame()

        # Participant counts per searched term / ontology
        keys_present = set(self.pids.keys())
        valid_keys = {"dterm", "hpo", "cterm", "cabbr", "icd10", "custom"}
        if not (keys_present & valid_keys):
            raise RuntimeError("No valid cohorts to summarize found")

        self.ont_vcount = {}

        # custom
        if "custom" in keys_present:
            tmp = pd.DataFrame({"participant_id": self.pids["custom"]})
            tmp["diag"] = "custom"
            self.ont_vcount["custom"] = (
                tmp.drop_duplicates(["participant_id", "diag"]).diag.value_counts()
            )

        # dterm (rare disease terms)
        if "dterm" in keys_present and hasattr(self, "dterm_table"):
            self.ont_vcount["dterm"] = (
                self.dterm_table
                .drop_duplicates(["participant_id", "normalised_specific_disease"])
                .normalised_specific_disease
                .value_counts()
            )

        # hpo
        if "hpo" in keys_present and hasattr(self, "hpo_table"):
            hpo_uniq = (
                self.hpo_table
                .drop_duplicates(["participant_id", "normalised_hpo_id"])
                .copy()
            )
            hpo_uniq["term"] = hpo_uniq["normalised_hpo_id"] + " (" + hpo_uniq["normalised_hpo_term"] + ")"
            self.ont_vcount["hpo"] = hpo_uniq["term"].value_counts()

        # cancer terms (disease_type)
        if "cterm" in keys_present and hasattr(self, "cterm_table"):
            self.ont_vcount["cterm"] = (
                self.cterm_table
                .drop_duplicates(["participant_id", "disease_type"])
                .disease_type
                .value_counts()
            )

        # cancer abbreviations (study_abbreviation)
        if "cabbr" in keys_present and hasattr(self, "cabbr_table"):
            self.ont_vcount["cabbr"] = (
                self.cabbr_table
                .drop_duplicates(["participant_id", "study_abbreviation"])
                .study_abbreviation
                .value_counts()
            )

        # icd10
        if "icd10" in keys_present and hasattr(self, "icd10_table"):
            # ensure lookup loaded
            if not hasattr(self, "icd10_lookup"):
                self.load_icd10_data()

            un_icd10 = self.icd10_table.drop_duplicates(["participant_id", "code"])

            # meanings
            full_icd10 = un_icd10.merge(self.icd10_lookup, on="code", how="left")
            self.ont_vcount["icd10_full"] = full_icd10["meaning"].value_counts(dropna=False)

            # simple codes (e.g., C34, D05)
            simple = un_icd10.copy()
            simple["simple_code"] = simple["code"].astype(str).str.extract(r"([A-Z][0-9]{2})")
            self.ont_vcount["icd10_simple"] = (
                simple["simple_code"].value_counts(dropna=False).reset_index()
            )

            # overlap matrix: participants with >1 simple code
            counts = (
                simple.dropna(subset=["simple_code"])
                    .groupby("participant_id", as_index=False)
                    .agg(n_codes=("simple_code", "nunique"))
            )
            multparts = counts.loc[counts["n_codes"] > 1, "participant_id"]
            self.icd10_overlap_matrix = (
                simple[simple["participant_id"].isin(multparts)]
                .groupby(["participant_id", "simple_code"])
                .size()
                .unstack(fill_value=0)
            )

        # Provenance
        self.provenance.append({
            "action": "ontology_stats",
            "sources_considered": sorted(list(keys_present & valid_keys)),
            "ont_summary_rows": int(len(self.ont_summary)),
            "vcount_keys": sorted(list(self.ont_vcount.keys())),
            "has_icd10_overlap_matrix": bool(hasattr(self, "icd10_overlap_matrix"))
        })


    def select_single_ca_sample(self, limit_to_featdict=False):
        """This function attempts to select a single cancer sample from
        participants with multiple samples. In those cases it first removes
        samples that are not in cancer_analysis. If there are still remaining
        samples it will remove non-Primary tumour samples. Finally, if there 
        are still multiple primary tumour samples it will select one with the 
        highest tumour_purity, and highest coverage_homogeneity in case of ties.
        
        With the limit_to_featdict flag samples are first filtered to limit them
        to the associated cancer study abbreviation (cancer_abbr) or disease type
        (cancer_terms).
        """
        self.concat_all()

        if self.all_cancer_samples is None or self.all_cancer_samples.empty:
            warnings.warn("No cancer samples available; call ca_sample_data() first.")
            return

        subcounts = self.all_cancer_samples["participant_id"].value_counts()
        dups = subcounts[subcounts > 1].index

        to_drop = []
        for pid in dups:
            p_samps = self.all_cancer_samples[self.all_cancer_samples["participant_id"] == pid].copy()
            dropped_for_pid = []

            # optional: limit to featdict
            if limit_to_featdict:
                if hasattr(self, "cterms"):
                    miss = p_samps.loc[~p_samps["disease_type"].isin(self.cterms), "tumour_sample_platekey"].tolist()
                    for pk in miss:
                        dropped_for_pid.append({"participant_id": pid, "platekey": pk, "reason": "not_matching_disease_type"})
                if hasattr(self, "cabbr"):
                    miss = p_samps.loc[
                        ~p_samps["study_abbreviation"].isin(self.cabbr), "tumour_sample_platekey"
                    ].tolist()
                    for pk in miss:
                        dropped_for_pid.append({"participant_id": pid, "platekey": pk, "reason": "not_matching_abbreviation"})
                if dropped_for_pid:
                    to_drop += [d["platekey"] for d in dropped_for_pid]
                    p_samps = p_samps[~p_samps["tumour_sample_platekey"].isin(to_drop)]

            # 1) present in cancer_analysis (nsv4 path present preferred)
            no_interp = p_samps.loc[
                p_samps["nsv4_somatic_small_variants_annotation_vcf"].isna(), "tumour_sample_platekey"
            ].tolist()
            for pk in no_interp:
                dropped_for_pid.append({"participant_id": pid, "platekey": pk, "reason": "not_in_cancer_analysis"})
            to_drop += no_interp
            p_samps = p_samps[~p_samps["tumour_sample_platekey"].isin(no_interp)]

            # 2) PRIMARY
            if len(p_samps) > 1:
                non_primary = p_samps.loc[p_samps["tumour_type"] != "PRIMARY", "tumour_sample_platekey"].tolist()
                for pk in non_primary:
                    dropped_for_pid.append({"participant_id": pid, "platekey": pk, "reason": "non_primary"})
                to_drop += non_primary
                p_samps = p_samps[~p_samps["tumour_sample_platekey"].isin(non_primary)]

            # 3) purity > coverage_homogeneity tie-break
            if len(p_samps) > 1:
                keep_idx = (
                    p_samps.sort_values(["tumour_purity", "coverage_homogeneity"], ascending=[False, False]).index[0]
                )
                keep_pk = p_samps.loc[keep_idx, "tumour_sample_platekey"]
                drop_rest = p_samps.loc[p_samps.index != keep_idx, "tumour_sample_platekey"].tolist()
                for pk in drop_rest:
                    dropped_for_pid.append(
                        {"participant_id": pid, "platekey": pk, "reason": "lower_purity_or_coverage"}
                    )
                to_drop += drop_rest

            if dropped_for_pid:
                self.provenance.extend(dropped_for_pid)

        # update platekeys + tables
        if self.platekeys is not None:
            self.platekeys = self.platekeys[~self.platekeys.isin(to_drop)]
        for key, table in list(self.cancer_samples.items()):
            self.cancer_samples[key] = table.loc[~table["tumour_sample_platekey"].isin(to_drop)].copy()
        self.sample_tables["cancer_samples"] = self.cancer_samples

        # rebuild merged + normalised views
        self.concat_all()
        self.build_normalized_views()


    def apply_survdat(self):
        from gelpack.survival import Survdat
        self.concat_all()
        self.initialize_survdat(
            impute=False, 
            df=None)
        