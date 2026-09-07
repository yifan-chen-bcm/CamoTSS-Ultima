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

All changes are confined to ``CamoTSS/utils/get_counts.py``



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
See *Known limitations* before counting them.



Citation
========

Please cite the original CamoTSS paper. This fork adds no new method.

  Hou, R., Hon, C.C. & Huang, Y. CamoTSS: analysis of alternative transcription
  start sites for cellular phenotypes and regulatory patterns from 5' scRNA-seq
  data. *Nat Commun* **14**, 7240 (2023).
  https://doi.org/10.1038/s41467-023-42636-1
