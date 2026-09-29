=========================================================================
CamoTSS-Ultima: alternative TSS analysis for Ultima UG100 5' scRNA-seq
=========================================================================

A fork of CamoTSS_ that makes the transcript-cluster (TSS) pipeline run on
**single-end, antisense, TSO-trimmed cDNA reads** -- the read geometry produced
by a 10x Genomics 5' universal library sequenced on an Ultima UG 100 and
converted to a synthetic R2 BAM for Cell Ranger.



.. _CamoTSS: https://github.com/StatBiomed/CamoTSS

:Upstream: CamoTSS 0.1.7 (StatBiomed/CamoTSS)
:This fork: 0.1.7+ultima.1
:Maintainer: Yifan Chen <yifan.chen@bcm.edu>
:License: Apache 2.0 (inherited from upstream)

Please report problems **with this fork** on this repository. Issues with
CamoTSS itself belong upstream_.

.. _upstream: https://github.com/StatBiomed/CamoTSS/issues


Scope
=====

- ``--mode TC`` is the supported and tested path.
- ``--mode CTSS`` and ``--mode TC+CTSS`` are **not usable on this data** and are
  left exactly as upstream. The CTSS window scan requires the 14-16 bp TSO soft
  clip to identify unencoded G; those reads carry no TSO, so every cluster
  returns an empty window list.
- Tested on human GRCh38 (GENCODE v44), Cell Ranger ``multi`` output.


============================================================================== ========================================================================================================
Upstream assumption                                                            Ultima UG100 / 10x 5' universal R2 BAM
============================================================================== ========================================================================================================
Reads arrive mated; R1 is the fragment 5' end                                  single-end records, no mate
cDNA read is sense to the gene                                                 antisense to the gene (Converted by Ultima virtual pair-end reads)
13 bp TSO is inside the read, as a 14-16 bp soft clip encoding the unencoded G TSO, barcode reads, and 5' Gs are trimmed
Cap site is the read 5' terminus                                               cap site is the read **3'** terminus (median +3 bp from an annotated TSS)
============================================================================== ========================================================================================================

All changes to the upstream pipeline are confined to ``CamoTSS/utils/get_counts.py``.
This fork also adds one downstream tool, ``CamoTSS/filter_quantify_cluster.py``
(see `Filtering and exact re-counting`_).



Limitations
=============

**Clusters are no longer false-positive filtered by a trained model.** The
logistic regression in upstream CamoTSS scores four features: ``UMI_count``,
``SD``, ``summit_UMI_count`` and ``unencoded_G_percent``. The last one is
derived from the TSO soft clip, which does not exist in Ultima virtual split reads.

Consequences for  analysis:

- What this fork produces is unsupervised clustering of UMI-deduplicated read
  3' termini, not model-scored cap sites. There is no learned component left in
  the ``TC`` pipeline.
- For UMIs have their 3' terminus more than 100 bp from any annotated
  TSS, those are either genuinely unannotated starts or non-cap 5' ends
  (internal priming, degradation), and this build cannot tell them apart.


Installation
============

.. code-block:: bash

   pip install -U git+https://github.com/yifan-chen-bcm/CamoTSS-Ultima

.. warning::

   Do not install this alongside upstream ``CamoTSS`` in the same environment.
   The distribution names differ (``CamoTSS-Ultima`` vs ``CamoTSS``) but both
   provide the same importable ``CamoTSS`` package and the same ``CamoTSS``
   console script, so they will silently overwrite each other. Uninstall one
   before installing the other.

Verify which build you have:

.. code-block:: bash

   CamoTSS --version
   # CamoTSS (CamoTSS-Ultima) 0.1.7+ultima.1 -- fork of CamoTSS 0.1.7




Usage
=====

.. code-block:: bash

   CamoTSS --mode TC \
     --gtf   gencode.v44.gtf.gz \
     --bam   merged.tagged.bam \
     -c      barcodes_camotss.tsv \
     -r      GRCh38.primary_assembly.genome.fa \
     -o      camotss_out \
     --nproc 12 --minCount 50 --maxReadCount 10000


Output
======

============================== =========================================================
File                           Contents
============================== =========================================================
``count/fourFeature.csv``      every candidate cluster with its four features. **This is
                               where your filtering now happens**
``count/afterfiltered.csv``    clusters surviving the (bypassed) filter -- currently
                               identical to the above
``count/scTSS_count_all.h5ad`` cell x TSS-cluster counts, all clusters
``count/scTSS_count_two.h5ad`` restricted to genes with >=2 clusters separated by
                               ``--clusterDistance``; this is the input for differential
                               TSS usage
``count/fetch_reads.pkl``      per-gene ``(position, CB, cigar)`` tuples after UMI collapse
============================== =========================================================

Clusters named ``<gene_id>_newTSS`` did not match any annotated transcript TSS.
See *Limitations* before counting them. ``CamoTSS-filter`` drops them unless their
summit is within ``--max-dist`` of an annotated TSS.

Filtering and exact re-counting
===============================

``CamoTSS-filter`` (``CamoTSS/filter_quantify_cluster.py``) is run **after**
``CamoTSS --mode TC``. It replaces the bypassed logistic filter with an
annotation-based one, then re-counts UMIs for the surviving clusters directly
from the BAM.

Why it exists
-------------

- **Filter only highly confident TSS** The script filters and only retain highly confident TSS withtin ``--max-dist`` (default 50bp) of an annotated TSS site.
- **Saturated counts.** CamoTSS stops fetching at ``--maxReadCount`` reads per
  gene (default 10000), so ``scTSS_count_*.h5ad`` under-counts highly expressed
  genes. This tool reads the BAM itself and has no read ceiling.

What it does
------------

==========  ============================================  ==========================================
Step        Input                                         Result
==========  ============================================  ==========================================
0 restrict  ``count/scTSS_count_two.h5ad`` (optional)     only the multi-TSS clusters CamoTSS kept
1 filter    ``count/fourFeature.csv`` +                   clusters whose summit is within
            ``ref_file/ref_TSS.tsv``                      ``--max-dist`` bp of an annotated TSS of
                                                          the same gene, in genes that keep at
                                                          least ``--min-tss`` clusters
2 count     filtered clusters + BAM + GTF + cell list     exact cell x cluster UMI matrix
==========  ============================================  ==========================================

Step 2 only runs when ``--bam`` is given. Omit ``--bam`` to stop after filtering.

For counting, a read is kept if it has MAPQ >= ``--mapq``, is antisense to the
gene, has a ``GX`` tag equal to the gene, and has ``CB`` and ``UB`` tags with a
``CB`` in ``--meta``. Reads are collapsed to one molecule per ``(CB, UB)`` at its
most 5' position (the read's 3' terminus). Each molecule is assigned to the
cluster interval that contains that position.
Usage
-----

Filter only:

.. code-block:: bash

   CamoTSS-filter \
     --two-h5ad    camotss_out/count/scTSS_count_two.h5ad \
     --fourfeature camotss_out/count/fourFeature.csv \
     --ref-tss     camotss_out/ref_file/ref_TSS.tsv \
     --outdir      camotss_out/quantify_two

Filter and re-count:

.. code-block:: bash

   CamoTSS-filter \
     --two-h5ad    camotss_out/count/scTSS_count_two.h5ad \
     --fourfeature camotss_out/count/fourFeature.csv \
     --ref-tss     camotss_out/ref_file/ref_TSS.tsv \
     --bam         merged.tagged.bam \
     --gtf         gencode.v44.gtf.gz \
     --meta        barcodes_camotss.tsv \
     --outdir      camotss_out/quantify_two \
     --nproc 12


Arguments
---------

================= ========== ================================================================
Argument          Default    Description
================= ========== ================================================================
``--fourfeature`` required   CamoTSS ``count/fourFeature.csv``
``--ref-tss``     required   CamoTSS ``ref_file/ref_TSS.tsv``
``--outdir``      required   output directory
``--two-h5ad``    none       CamoTSS ``count/scTSS_count_two.h5ad``. Restricts the analysis to
                             the clusters CamoTSS kept as multi-TSS before any other filter.
                             Omit it to start from every cluster in ``fourFeature.csv``
``--max-dist``    50         maximum distance (bp) from a cluster summit to an annotated TSS
                             of the same gene
``--min-tss``     2          minimum number of surviving clusters a gene must keep
``--bam``         none       BAM to re-count from. Passing it turns on step 2
``--gtf``         none       GTF (plain or gzipped) matching the BAM's ``GX`` tags. Required
                             with ``--bam``; the gene strand is taken from here
``--meta``        none       pass-filter cells. Required with ``--bam``. Accepts the CamoTSS
                             ``-c`` barcode file (``cell_id`` column), a CSV with ``cell_id``
                             or ``cell_barcode``, or a headerless one-barcode-per-line list.
                             Extra columns are copied to ``.obs``
``--nproc``       8          worker processes for counting
``--mapq``        255        minimum MAPQ (STAR / Cell Ranger give 255 to unique alignments)
================= ========== ================================================================

Output
------

========================================== =======================================================
File                                       Contents
========================================== =======================================================
``<outdir>/confident_TSS.csv``             filtered clusters: the ``fourFeature.csv`` columns plus
                                           ``cluster_start``, ``cluster_end``, ``known_TSS``,
                                           ``dist_signed`` (summit - known TSS), ``dist_abs`` and
                                           ``n_cluster``
``<outdir>/filtered_clusters_counts.h5ad`` written with ``--bam`` only. ``X``: cells x clusters
                                           UMI counts (sparse int32). ``obs``: the ``--meta``
                                           table, in its row order. ``var``: the cluster table,
                                           indexed by ``cluster_id``, plus ``chrom``, ``strand``,
                                           ``gene_name``, ``umi_total`` and ``cells_detected``.
                                           ``uns['params']``: the run settings
========================================== =======================================================

Clusters whose gene is missing from the GTF are dropped before counting, and the
``--min-tss`` rule is applied again afterwards. As a result, the h5ad can have
fewer clusters than ``confident_TSS.csv``. Always join the two on ``cluster_id``
/ ``clusterID``, never on row position.


Citation
========

Please cite the original CamoTSS paper.

  Hou, R., Hon, C.C. & Huang, Y. CamoTSS: analysis of alternative transcription
  start sites for cellular phenotypes and regulatory patterns from 5' scRNA-seq
  data. *Nat Commun* **14**, 7240 (2023).
  https://doi.org/10.1038/s41467-023-42636-1
