# Source inspection: `GSE125367_GPL21145`

- Generated at: `2026-08-05T15:27:48.162668+00:00`
- GSE: `GSE125367`
- GPL: `GPL21145`

## Layout

```json
{
  "beta_id_column": "GSM_ID",
  "betas_bytes": 435108762,
  "betas_exists": true,
  "betas_path": "/data/projects/methyl-burden-score/data/raw/cpgcorpus/GSE125367/GPL21145/betas/QCDPB.arrow",
  "exists": true,
  "metadata_bytes": 21018,
  "metadata_exists": true,
  "metadata_path": "/data/projects/methyl-burden-score/data/raw/cpgcorpus/GSE125367/GPL21145/metadata/metadata.arrow",
  "orientation": "samples_as_rows_probes_as_columns",
  "processing_level": "QCDPB",
  "raw_gpl_root": "/data/projects/methyl-burden-score/data/raw/cpgcorpus/GSE125367/GPL21145"
}
```

## Sample alignment

```json
{
  "beta_only_count": 0,
  "intersection_size": 44,
  "metadata_only_count": 0,
  "n_beta_samples": 44,
  "n_duplicate_beta_ids": 0,
  "n_duplicate_metadata_ids": 0,
  "n_metadata_samples": 44,
  "n_unique_beta_ids": 44,
  "n_unique_metadata_ids": 44,
  "perfect_alignment": true
}
```

## Beta value QC

```json
{
  "duplicate_probe_names": 0,
  "fraction_above_1": 0.0,
  "fraction_below_0": 0.0,
  "max": 0.9949491849606493,
  "mean": 0.6458271211423645,
  "min": 0.005998052032217561,
  "missing_fraction": 0.13360818233156177,
  "n_probes": 865919,
  "n_samples": 44,
  "n_values": 38100436,
  "sample_mean_beta": {
    "max": 0.66234772950111,
    "median": 0.6443766018224477,
    "min": 0.6335836552318038
  },
  "sample_missing_fraction": {
    "max": 0.14243941985335812,
    "median": 0.13387972778054297,
    "min": 0.12581430826670856
  },
  "std": 0.37473592028527825,
  "zero_variance_loci": 467
}
```

## Metadata counts (sanitized fields)

```json
{
  "age_numeric": {
    "field": "age (years):ch1",
    "max": 16.0,
    "median": 8.0,
    "min": 0.5,
    "n_missing": 0,
    "n_numeric": 44,
    "n_total": 44,
    "n_unparsed": 0
  },
  "column_names": [
    "GSM_ID",
    "title",
    "geo_accession",
    "status",
    "submission_date",
    "last_update_date",
    "type",
    "channel_count",
    "source_name_ch1",
    "organism_ch1",
    "characteristics_ch1",
    "characteristics_ch1.1",
    "molecule_ch1",
    "extract_protocol_ch1",
    "label_ch1",
    "label_protocol_ch1",
    "taxid_ch1",
    "hyb_protocol",
    "scan_protocol",
    "description",
    "data_processing",
    "platform_id",
    "contact_name",
    "contact_email",
    "contact_laboratory",
    "contact_institute",
    "contact_address",
    "contact_city",
    "contact_state",
    "contact_zip/postal_code",
    "contact_country",
    "supplementary_file",
    "supplementary_file.1",
    "data_row_count",
    "age (years):ch1",
    "Sex:ch1"
  ],
  "fields": {
    "Sex:ch1": {
      "F": 21,
      "M": 23
    },
    "age (years):ch1": {
      "0.5": 1,
      "1": 4,
      "1.08": 1,
      "10": 3,
      "11": 4,
      "11.33": 1,
      "11.92": 1,
      "12": 2,
      "13": 1,
      "14": 3,
      "15": 1,
      "16": 2,
      "2": 1,
      "3": 1,
      "3.5": 1,
      "3.83": 1,
      "4": 7,
      "5": 1,
      "6": 1,
      "7": 2,
      "8": 3,
      "9": 2
    },
    "molecule_ch1": {
      "genomic DNA": 44
    },
    "organism_ch1": {
      "Homo sapiens": 44
    },
    "platform_id": {
      "GPL21145": 44
    },
    "type": {
      "genomic": 44
    }
  },
  "n_columns": 36
}
```

## Warnings

_None._
