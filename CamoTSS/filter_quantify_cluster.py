#!/usr/bin/env python
"""
Yifan Chen
yifan.chen@bcm.edu
9-16-2026

Restrict CamoTSS-Ultima TSS clusters to the multi-TSS set CamoTSS kept
in scTSS_count_two.h5ad, filter them, and re-count their UMIs exactly from the BAM.
Writes a cell x cluster UMIs h5ad file


Script follows:

Step 0 restrict  scTSS_count_two.h5ad          -> clusters CamoTSS kept as multi-TSS
Step 1 filter    fourFeature.csv + ref_TSS.tsv -> filtered cluster set
Step 2 count     filtered clusters + BAM       -> exact cell x cluster UMIs

Step 0: scTSS_count_two.h5ad is a column subset of scTSS_count_all.h5ad, pruned in
        produce_sclevel() so clusters of the same gene sit at least --clusterDistance
        (default 300 bp) apart.

Step 1: Keep clusters whose summit lies within --max-dist bp of an annotated TSS of
        the SAME gene (ref_TSS.tsv). CamoTSS-Ultima bypasses the logistic FP filter,
        so novel TSS calls are not trustworthy and proximity to an annotated TSS is
        the only confidence signal left. Then keep only genes retaining >= --min-tss
        clusters.

Step 2: Reads the BAM directly, so unlike scTSS_count_*.h5ad there is no
        --maxReadCount ceiling; CamoTSS stops fetching at 10000 reads per gene, which
        saturates every highly expressed gene.

Read geometry is the Ultima UG100 / 10x 5' universal convention: single-end records
antisense to the gene, cap site at the read 3' terminus.

fourFeature.csv columns: UMI_count, SD, summit_UMI_count, unencoded_G_percent,
NO.TSS, gene_id, summit_position

Example:

    python filter_quantify_two.py \
      --two-h5ad    camotss_out/count/scTSS_count_two.h5ad \
      --fourfeature camotss_out/count/fourFeature.csv \
      --ref-tss     camotss_out/ref_file/ref_TSS.tsv \
      --bam         merged.tagged.bam \
      --gtf         genes.gtf.gz \
      --meta        barcodes_camotss.tsv \
      --outdir      camotss_out/quantify_two \
      --nproc 12

Omit --bam to stop after filtering.
"""

import argparse
import collections
import gzip
import multiprocessing as mp
import os
import re
import sys

import anndata
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp

DEFAULT_MAX_DIST = 50
DEFAULT_MIN_TSS = 2
DEFAULT_MAPQ = 255
DEFAULT_NPROC = 8


# STEP 0 RESTRICT TO THE CLUSTERS CAMOTSS KEPT AS MULTI-TSS ------------------------------------
def h5ad_cluster_ids(path):
    """clusterIDs held by an scTSS_count_two.h5ad.

    The h5ad is keyed by transcript name,  so clusterID has to be rebuilt from
    the TSS_start / TSS_end columns to match fourFeature.csv.
    """
    ad = anndata.read_h5ad(path, backed="r")
    try:
        var = ad.var
        missing = {"gene_id", "TSS_start", "TSS_end"} - set(var.columns)
        if missing:
            sys.exit("%s var/ is missing required column(s): %s" % (path, sorted(missing)))
        ids = (var.gene_id.astype(str) + "*"
               + var.TSS_start.astype("int64").astype(str) + "_"
               + var.TSS_end.astype("int64").astype(str))
    finally:
        if ad.file is not None:
            ad.file.close()
    return set(ids)


# STEP 1 TSS CLUSTER FILTERING -----------------------------------------------------------------
def nearest_annotated_tss(clusters, ref):
    """Nearest annotated TSS of the SAME gene, per cluster.
    Returns known_TSS and dist_signed = summit_position - known_TSS.
    """
    by_gene = {g: np.sort(v.values) for g, v in ref.groupby("gene_id").TSS} # group TSS by gene_id in ref
    known = np.full(len(clusters), np.nan)

    for i, (gid, pos) in enumerate(zip(clusters.gene_id.values,
                                       clusters.summit_position.values)):
        arr = by_gene.get(gid)
        if arr is None or not len(arr):
            continue
        j = np.searchsorted(arr, pos)
        cand = arr[max(j - 1, 0):min(j + 1, len(arr)) + 1]
        known[i] = cand[np.argmin(np.abs(cand - pos))]
    return known, clusters.summit_position.values - known


def filter_clusters(fourfeature, ref_tss, outdir, two_h5ad=None,
                    max_dist=DEFAULT_MAX_DIST, min_tss=DEFAULT_MIN_TSS):
    """fourFeature.csv (+ optionally an h5ad) -> confident multi-TSS cluster table,
    also written to <outdir>/confident_TSS.csv. This is the output of steps 0 and 1:
    every cluster that is close enough to an annotated TSS to be trusted, in a gene
    that retains at least min_tss of them. """

    ff = pd.read_csv(fourfeature)
    ff = ff.rename(columns={ff.columns[0]: "clusterID"})
    missing = {"gene_id", "summit_position"} - set(ff.columns)
    if missing:
        sys.exit("fourFeature.csv is missing required column(s): %s" % sorted(missing))

    coord = ff.clusterID.str.split("*", n=1).str[1]
    bad = coord.isna()
    if bad.any():
        sys.exit("%d clusterIDs have no '*': %s" % (bad.sum(), ff.clusterID[bad].head().tolist()))
    ff["cluster_start"] = coord.str.split("_").str[0].astype(int)
    ff["cluster_end"] = coord.str.split("_").str[1].astype(int)

    # step 0, before anything else, so the reported filter counts are relative to it
    if two_h5ad is not None:
        keep = h5ad_cluster_ids(two_h5ad)
        n = len(ff)
        ff = ff[ff.clusterID.isin(keep)].copy()
        print("restricted to %s: %d of %d fourFeature clusters kept, %d of %d h5ad "
              "columns matched" % (os.path.basename(two_h5ad), len(ff), n,
                                   ff.clusterID.nunique(), len(keep)), flush=True)
        if ff.clusterID.nunique() < len(keep):
            print("  unmatched h5ad columns are the '<gene>_newTSS' keys, which collapse "
                  "every unannotated cluster of a gene onto one column, so their "
                  "TSS_start/TSS_end need not equal any single cluster", flush=True)
        if not len(ff):
            sys.exit("no fourFeature cluster matched the h5ad")

    ref = pd.read_csv(ref_tss, sep="\t")
    ff["known_TSS"], ff["dist_signed"] = nearest_annotated_tss(ff, ref)
    ff["dist_abs"] = ff.dist_signed.abs()

    # keep clusters within max_dist bp of an annotated TSS
    ff = ff[ff.dist_abs <= max_dist].copy()

    # keep only genes carrying at least min_tss clusters
    ff["n_cluster"] = ff.groupby("gene_id").gene_id.transform("size")
    ff = ff[ff.n_cluster >= min_tss]

    if not len(ff):
        sys.exit("no clusters survived filtering")

    ff = ff.sort_values(["gene_id", "cluster_start"]).reset_index(drop=True)

    os.makedirs(outdir, exist_ok=True)
    ff.to_csv(os.path.join(outdir, "confident_TSS.csv"), index=False)
    return ff


# STEP 2 EXACT RECOUNTING FROM THE BAM ---------------------------------------------------------
_G = {}   # per-worker globals, filled by _init once in each pool process


def _open_text(path):
    """Accept both plain and gzipped GTFs, so sniff the magic bytes."""
    with open(path, "rb") as fh:
        gzipped = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gzipped else open(path, "rt")


def read_meta(meta):
    """Pass-filter cells, in matrix ROW order, plus the barcode column name.

    Accepts the CamoTSS --cellbarcodeFile (tab-separated, cell_id column), a CSV
    with cell_id or cell_barcode, or a headerless one-barcode-per-line list. Any
    further columns are carried through to .obs.
    """
    sep = "\t" if meta.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(meta, sep=sep)
    col = next((c for c in ("cell_id", "cell_barcode") if c in df.columns), None)
    if col is None:
        if df.shape[1] == 1:                    # headerless list: header became a barcode
            df = pd.DataFrame({"cell_id": [df.columns[0]] + list(df.iloc[:, 0].astype(str))})
            col = "cell_id"
        else:
            sys.exit("--meta needs a cell_id or cell_barcode column (found: %s)" % list(df.columns))
    df[col] = df[col].astype(str)
    if df[col].duplicated().any():
        sys.exit("--meta has duplicate barcodes in %s" % col)
    return df, col


def gene_info(gtf):
    """gene_id -> (chrom, start, end, strand, name), from the GTF's gene lines."""
    info = {}
    with _open_text(gtf) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.split("\t", 9)
            if p[2] != "gene":
                continue
            gid = re.search(r'gene_id "([^"]+)', p[8]).group(1)
            nm = re.search(r'gene_name "([^"]+)', p[8])
            info[gid] = (p[0], int(p[3]), int(p[4]), p[6], nm.group(1) if nm else gid)
    return info


def _init(bam, cellidx, clusters, mapq):
    """Runs once per worker. The BAM handle cannot be pickled, so each process
    opens its own here rather than receiving one."""
    _G["bam"] = pysam.AlignmentFile(bam, "rb")
    _G["cellidx"] = cellidx
    _G["clusters"] = clusters
    _G["mapq"] = mapq


def _one_gene(gid):
    """One gene: fetch its span, collapse reads to one position per (CB, UB),
    assign each UMI to the cluster interval containing it.
    Returns (rows, cols, vals) triplets for the sparse matrix."""
    chrom, strand, ivs = _G["clusters"][gid]        # (chrom, strand, [(start, end, col)])
    lo = min(s for s, _, _ in ivs)
    hi = max(e for _, e, _ in ivs)
    cellidx = _G["cellidx"]
    plus = strand == "+"

    best = {}
    for r in _G["bam"].fetch(chrom, max(lo - 200, 0), hi + 200):
        # PCR duplicates are deliberately kept: the (CB, UB) collapse below reduces
        # them to one molecule anyway, and more reads pin down the most-5' position
        # of that molecule better.
        if r.mapping_quality < _G["mapq"]:
            continue
        if r.is_reverse != plus:                    # keep reads ANTISENSE to the gene
            continue
        if not r.has_tag("GX") or r.get_tag("GX") != gid:
            continue
        if not r.has_tag("UB") or not r.has_tag("CB"):
            continue
        ci = cellidx.get(r.get_tag("CB"))
        if ci is None:                              # not a pass-filter cell
            continue
        p = r.reference_start if plus else r.reference_end - 1  # the read's 3' END, which
                                                                # for antisense (R2) reads is
                                                                # the TSS-proximal end
        k = (ci, r.get_tag("UB"))
        q = best.get(k)
        if q is None:
            best[k] = p
        else:
            best[k] = min(q, p) if plus else max(q, p)   # most 5' in transcript orientation

    acc = collections.Counter()
    for (ci, _), p in best.items():
        for s, e, col in ivs:
            if s <= p <= e:
                acc[(ci, col)] += 1
                break                               # intervals are disjoint
    if not acc:
        return np.empty(0, np.int32), np.empty(0, np.int32), np.empty(0, np.int32)
    rows = np.fromiter((k[0] for k in acc), np.int32, len(acc))
    cols = np.fromiter((k[1] for k in acc), np.int32, len(acc))
    vals = np.fromiter(acc.values(), np.int32, len(acc))
    return rows, cols, vals


def count_clusters(ff, bam, gtf, meta, outdir, nproc=DEFAULT_NPROC, mapq=DEFAULT_MAPQ,
                   min_tss=DEFAULT_MIN_TSS, params=None):
    """Exact cell x cluster UMI counts for the table filter_clusters returned."""
    os.makedirs(outdir, exist_ok=True)

    meta_df, col = read_meta(meta)
    cells = list(meta_df[col])                      # matrix ROW order
    cellidx = {c: i for i, c in enumerate(cells)}

    # Strand comes from the GTF, not ref_TSS.tsv. It decides which reads _one_gene
    # keeps and which end of a read is the 5' end, so it has to come from the same
    # annotation the BAM's GX tags were made against.
    ff = ff.copy()
    info = gene_info(gtf)
    ff["chrom"] = ff.gene_id.map(lambda g: info[g][0] if g in info else None)
    ff["strand"] = ff.gene_id.map(lambda g: info[g][3] if g in info else None)
    ff["gene_name"] = ff.gene_id.map(lambda g: info[g][4] if g in info else g)

    n = len(ff)
    ff = ff[ff.strand.notna()]
    if len(ff) < n:
        print("dropped %d clusters whose gene is not in the GTF" % (n - len(ff)))

    # That drop removes SINGLE clusters, so a gene can fall back to one and break
    # the >1 rule step 1 established. Re-apply it, then recompute n_cluster on
    # what actually reaches the matrix.
    ff = ff[ff.groupby("gene_id").gene_id.transform("size") >= min_tss]
    ff = ff.sort_values(["gene_id", "cluster_start"]).reset_index(drop=True)
    ff["n_cluster"] = ff.groupby("gene_id").gene_id.transform("size")
    if not len(ff):
        sys.exit("no clusters left to count")

    # col is assigned LAST, after every filter and the final sort, so ff's row order
    # IS the matrix's column order, i.e. the .var order of the h5ad written below.
    # Always join on clusterID.
    ff["col"] = np.arange(len(ff))
    clusters = {}
    for gid, sub in ff.groupby("gene_id"):
        clusters[gid] = (sub.chrom.iloc[0], sub.strand.iloc[0],
                         list(zip(sub.cluster_start.astype(int),
                                  sub.cluster_end.astype(int),
                                  sub.col.astype(int))))
    gene_ids = sorted(clusters)
    print("Counting %d clusters over %d genes for %d cells"
          % (len(ff), len(gene_ids), len(cells)), flush=True)

    R, C, V = [], [], []
    with mp.Pool(nproc, initializer=_init,
                 initargs=(bam, cellidx, clusters, mapq)) as pool:
        for k, (r, c, v) in enumerate(pool.imap_unordered(_one_gene, gene_ids,
                                                          chunksize=8), 1):
            if len(r):
                R.append(r); C.append(c); V.append(v)
            if k % 2000 == 0:
                print("  %d/%d genes" % (k, len(gene_ids)), flush=True)

    if not R:                                       # np.concatenate raises on []
        M = sp.csr_matrix((len(cells), len(ff)), dtype=np.int32)
    else:
        M = sp.coo_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))),
                          shape=(len(cells), len(ff)), dtype=np.int32).tocsr()

    ff["umi_total"] = np.asarray(M.sum(0)).ravel()
    ff["cells_detected"] = np.asarray((M > 0).sum(0)).ravel()

    # obs = cells (matrix rows), var = clusters (matrix columns). Both orders were
    # fixed above -- meta's row order defined `cells`, ff's row order defined `col`
    # -- so nothing is reindexed here and X stays aligned to both.
    obs = meta_df.set_index(col).rename_axis("cell_id")
    for c in obs.columns:                 # category is what scanpy wants for groupby
        if obs[c].dtype == object and obs[c].nunique() < len(obs) / 2:
            obs[c] = obs[c].astype("category")

    var = ff.set_index("clusterID").rename_axis("cluster_id").drop(columns="col")
    for c in var.columns:                 # object columns cannot be written to h5ad
        if var[c].dtype == object:
            var[c] = var[c].astype(str)

    ad = anndata.AnnData(X=M, obs=obs, var=var)
    ad.uns["params"] = dict(params or {}, mapq=mapq, nproc=nproc,
                            bam=os.path.abspath(bam), gtf=os.path.abspath(gtf),
                            meta=os.path.abspath(meta))
    h5 = os.path.join(outdir, "filtered_clusters_counts.h5ad")
    ad.write_h5ad(h5, compression="gzip")

    print("total UMIs in clusters: %d" % M.sum())
    print("clusters with >0 UMIs: %d / %d" % (int((ff.umi_total > 0).sum()), len(ff)))
    print("wrote %s" % h5)
    return ad


# MAIN -----------------------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    # step 0 and 1
    p.add_argument("--two-h5ad", default=None,
                   help="CamoTSS scTSS_count_two.h5ad. Restricts the analysis to the "
                        "clusters CamoTSS kept as multi-TSS (pruned by --clusterDistance) "
                        "before any other filter. Omit to use every fourFeature cluster")
    p.add_argument("--fourfeature", required=True, help="CamoTSS count/fourFeature.csv")
    p.add_argument("--ref-tss", required=True, help="CamoTSS ref_file/ref_TSS.tsv")
    p.add_argument("--outdir", required=True)
    p.add_argument("--max-dist", type=int, default=DEFAULT_MAX_DIST,
                   help="keep clusters whose summit is within this many bp of an "
                        "annotated TSS of the same gene [default: %(default)s]")
    p.add_argument("--min-tss", type=int, default=DEFAULT_MIN_TSS,
                   help="keep only genes with at least this many surviving clusters "
                        "[default: %(default)s]")
    # step 2, only runs when --bam is given
    p.add_argument("--bam", default=None, help="give this to run the counting step")
    p.add_argument("--gtf", default=None, help="GTF (plain or gzipped) matching the BAM's GX tags")
    p.add_argument("--meta", default=None,
                   help="pass-filter cell metadata, must contain a cell_id column")
    p.add_argument("--nproc", type=int, default=DEFAULT_NPROC)
    p.add_argument("--mapq", type=int, default=DEFAULT_MAPQ,
                   help="STAR/CellRanger set MAPQ 255 for uniquely mapped reads")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    ff = filter_clusters(a.fourfeature, a.ref_tss, a.outdir, two_h5ad=a.two_h5ad,
                         max_dist=a.max_dist, min_tss=a.min_tss)
    print("Filtered to %d clusters over %d genes" % (len(ff), ff.gene_id.nunique()), flush=True)

    if a.bam:
        for req, flag in ((a.gtf, "--gtf"), (a.meta, "--meta")):
            if req is None:
                sys.exit("%s is required when --bam is given" % flag)
        count_clusters(ff, a.bam, a.gtf, a.meta, a.outdir,
                       nproc=a.nproc, mapq=a.mapq, min_tss=a.min_tss,
                       params={"max_dist": a.max_dist, "min_tss": a.min_tss,
                               "two_h5ad": os.path.abspath(a.two_h5ad) if a.two_h5ad else ""})
    else:
        print("no --bam given: stopped after filtering")


if __name__ == "__main__":
    main()
