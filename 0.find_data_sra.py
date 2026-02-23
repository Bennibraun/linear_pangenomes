import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import pandas as pd
import io
from pathlib import Path

def make_session():
    session = requests.Session()
    retry = Retry(total=5, backoff_factor=2, status_forcelist=[500, 502, 503, 504])
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    return session

def fetch_ena_portal(platform: str, cache_path: Path, taxon_id: int = 33208) -> pd.DataFrame:
    if cache_path.exists():
        print(f"Loading cached {cache_path}")
        return pd.read_parquet(cache_path)

    url = "https://www.ebi.ac.uk/ena/portal/api/search"
    params = {
        "result": "read_run",
        "query": f'tax_tree({taxon_id}) AND instrument_platform="{platform}" AND library_strategy="WGS"',
        "fields": "run_accession,sample_accession,study_accession,scientific_name,instrument_platform,library_strategy",
        "format": "tsv",
        "limit": 0,
    }
    print(f"Fetching {platform} from ENA Portal API...")
    session = make_session()
    r = session.get(url, params=params, timeout=120)
    r.raise_for_status()
    df = pd.read_csv(io.StringIO(r.text), sep="\t")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_path, index=False)
    print(f"Saved {len(df)} rows to {cache_path}")
    return df

cache_dir = Path("sra_cache")

df_ont      = fetch_ena_portal("OXFORD_NANOPORE", cache_dir / "metazoa_ont_wgs.parquet")
df_pb       = fetch_ena_portal("PACBIO_SMRT",     cache_dir / "metazoa_pacbio_wgs.parquet")
df_illumina = fetch_ena_portal("ILLUMINA",        cache_dir / "metazoa_illumina_wgs.parquet")

# df_lr = pd.concat([df_ont, df_pb])

# lr_counts       = df_lr.groupby("scientific_name")["sample_accession"].nunique().rename("lr_samples")
# illumina_counts = df_illumina.groupby("scientific_name")["sample_accession"].nunique().rename("illumina_samples")

# combined = pd.concat([lr_counts, illumina_counts], axis=1).dropna()
# hits = combined[(combined["lr_samples"] >= 10) & (combined["illumina_samples"] >= 50)]
# print(hits.sort_values("lr_samples", ascending=False).to_string())



df_lr = pd.concat([df_ont, df_pb])

# Find studies present in both LR and Illumina for the same species
lr_studies = df_lr[["study_accession", "sample_accession", "scientific_name"]].drop_duplicates()
il_studies = df_illumina[["study_accession", "sample_accession", "scientific_name"]].drop_duplicates()

# Join on study_accession
shared_studies = lr_studies.merge(il_studies, on=["study_accession", "scientific_name"], suffixes=("_lr", "_il"))

# Count unique samples per study per species
study_counts = (
    shared_studies
    .groupby(["scientific_name", "study_accession"])
    .agg(
        lr_samples=("sample_accession_lr", "nunique"),
        illumina_samples=("sample_accession_il", "nunique"),
    )
    .reset_index()
)

# Filter to your thresholds
hits = study_counts[
    (study_counts["lr_samples"] >= 10) &
    (study_counts["illumina_samples"] >= 50)
].sort_values("lr_samples", ascending=False)

print(hits.to_string())