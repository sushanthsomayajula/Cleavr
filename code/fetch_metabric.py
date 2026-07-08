"""
Fetch METABRIC (brca_metabric) clinical + expression data from the public
cBioPortal REST API, for cross-cohort validation of the TCGA-derived screen
results (results/screen_results.json).

Why the API instead of a bulk tarball: the cBioPortal datahub bulk-download
mirrors (S3 / download.cbioportal.org) were unreachable/returned AccessDenied
when tried directly. The REST API on www.cbioportal.org is the production
backend and responded fine, so this pulls only what's needed (clinical
attributes + a fixed gene list) through it instead.

METABRIC is Illumina microarray, not RNA-seq V2 RSEM like TCGA -- different
platform, different normalization. That's actually the point of a
cross-cohort check: if a gene's TNBC-subtype association shows up on two
different platforms/cohorts, that's a real signal, not a TCGA/RSEM artifact.

OUTPUT:
    ../brca_metabric/clinical_patient.csv   (patientId, ER_STATUS, PR_STATUS, HER2_STATUS, OS_STATUS, OS_MONTHS)
    ../brca_metabric/expression.csv         (sampleId, patientId, <gene columns>)

HOW TO RUN:
    cd ~/Desktop/tnbc-project/code
    python3 fetch_metabric.py
"""
import csv
import json
import os
import time
import urllib.request

API = "https://www.cbioportal.org/api"
STUDY = "brca_metabric"
OUT_DIR = "../brca_metabric"

# Same genes used in biomarker_pipeline.py (MARKER_GENES for subtype
# assignment) plus the 8 genes that were significant after FDR correction
# in results/screen_results.json. Entrez IDs pulled directly from
# brca_tcga/data_mrna_seq_v2_rsem.txt so they match exactly.
GENE_ENTREZ = {
    "AR": 367, "CCNE1": 898, "CDC6": 990, "EGFR": 1956, "FOXA1": 3169,
    "MET": 4233, "VIM": 7431, "ZEB1": 6935,          # marker genes
    "FTO": 79068, "POLD1": 5424, "POLD2": 5425, "POLE": 5426,
    "PRIM2": 5558, "RRM1": 6240, "RRM2": 6241, "TYMS": 7298,  # significant genes
}

PATIENT_ATTRS = {"OS_STATUS", "OS_MONTHS"}
SAMPLE_ATTRS = {"ER_STATUS", "PR_STATUS", "HER2_STATUS"}  # these are sample-level in METABRIC, not patient-level


def api_get(path):
    req = urllib.request.Request(f"{API}{path}", headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read())


def api_post(path, body):
    data = json.dumps(body).encode()
    req = urllib.request.Request(
        f"{API}{path}", data=data, method="POST",
        headers={"Accept": "application/json", "Content-Type": "application/json"},
    )
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read())


def fetch_clinical():
    print("Fetching patient-level clinical data (OS_STATUS, OS_MONTHS)...")
    prows = api_get(f"/studies/{STUDY}/clinical-data?clinicalDataType=PATIENT&pageSize=100000&pageNumber=0")
    print(f"  {len(prows)} raw patient clinical-data rows")

    by_patient = {}
    for r in prows:
        attr = r["clinicalAttributeId"]
        if attr not in PATIENT_ATTRS:
            continue
        pid = r["patientId"]
        by_patient.setdefault(pid, {})[attr] = r["value"]

    print("Fetching sample-level clinical data (ER_STATUS, PR_STATUS, HER2_STATUS)...")
    srows = api_get(f"/studies/{STUDY}/clinical-data?clinicalDataType=SAMPLE&pageSize=100000&pageNumber=0")
    print(f"  {len(srows)} raw sample clinical-data rows")

    for r in srows:
        attr = r["clinicalAttributeId"]
        if attr not in SAMPLE_ATTRS:
            continue
        # merge onto the same patientId record -- METABRIC is ~1 sample/patient
        pid = r["patientId"]
        by_patient.setdefault(pid, {})[attr] = r["value"]

    print(f"  {len(by_patient)} patients with at least one relevant attribute")
    return by_patient


def fetch_expression():
    print(f"Fetching expression for {len(GENE_ENTREZ)} genes (brca_metabric_mrna, microarray)...")
    entrez_ids = list(GENE_ENTREZ.values())
    entrez_to_symbol = {v: k for k, v in GENE_ENTREZ.items()}

    data = api_post(
        "/molecular-profiles/brca_metabric_mrna/molecular-data/fetch",
        {"entrezGeneIds": entrez_ids, "sampleListId": "brca_metabric_mrna"},
    )
    print(f"  {len(data)} raw expression data points")

    by_sample = {}
    for d in data:
        sample_id = d["sampleId"]
        patient_id = d.get("patientId", sample_id.rsplit("-", 1)[0] if "-" in sample_id else sample_id)
        gene = entrez_to_symbol.get(d["entrezGeneId"])
        if gene is None:
            continue
        row = by_sample.setdefault(sample_id, {"sampleId": sample_id, "patientId": patient_id})
        row[gene] = d["value"]

    print(f"  {len(by_sample)} samples with expression data")
    return by_sample


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    clinical = fetch_clinical()
    with open(f"{OUT_DIR}/clinical_patient.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["patientId", "ER_STATUS", "PR_STATUS", "HER2_STATUS", "OS_STATUS", "OS_MONTHS"])
        for pid, attrs in clinical.items():
            writer.writerow([pid, attrs.get("ER_STATUS", ""), attrs.get("PR_STATUS", ""),
                              attrs.get("HER2_STATUS", ""), attrs.get("OS_STATUS", ""),
                              attrs.get("OS_MONTHS", "")])
    print(f"Saved {OUT_DIR}/clinical_patient.csv")

    expression = fetch_expression()
    gene_cols = list(GENE_ENTREZ.keys())
    with open(f"{OUT_DIR}/expression.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["sampleId", "patientId"] + gene_cols)
        for row in expression.values():
            writer.writerow([row["sampleId"], row["patientId"]] + [row.get(g, "") for g in gene_cols])
    print(f"Saved {OUT_DIR}/expression.csv")

    print("\nDONE.")
