# Metadata standards

## Sample metadata

## Phenotype metadata
### Required columns

 - **phenotype_id** -  id of the molecular trait that has been quantified. This can be either the gene id (RNA-eq eQTLs), probe id (microarray eQTLs), transcript id (full-length transcript usage QTLs), splice junction id (Leafcutter), exon id (exon-level QTLs) or any other molecular trait that has been quantified.
 - **quant_id** - Currently only used for txrevise to normalise promoter/splicing/3'end usage.
-  **group_id** - Used for transcript usage, exon expression and splicing phenotypes. Overlapping phenotypes whose relative expression is quantified belong to the same group (e.g. alternative spliced exons form clusters in Leafcutter).
- **gene_id** - Ensembl gene id. Should be NA if the this information is not available (e.g. novel junction clusters in LeafCutter).
- **gene_start** - Start coordinate of the gene
- **gene_end** - End coordinate of the gene
-   **strand** - Strand of the gene
-   **gene_name** - Gene name extracted from Ensembl biomart.
-   **gene_type** - Gene type (protein coding, lincRNA, etc) extracted from Ensembl biomart.
-  **gene_version** - Ensembl gene version
-  **phenotype_pos** - Genomic position used to determine the centre point of the *cis* window for QTL mapping. 

### Optional columns
*   **phenotype_gc_content** - Percentage GC content of the quantified phenotype. Currently used as a covariate by cqn when normalising gene-level and exon-level counts. For genes, this is exracted directly from Ensembl biomart. For exons it is calculated using bedtools nuc command. 
**phenotype_length** - Length of the phenotype in base pairs. Currently used for gene-level and exon-level counts. Used by cqn and TPM normalisation techniques.

### Definition of the the *cis* window
The center point of the *cis* window is defined by the **phenotype_pos** column in the phenotype metadata file:

 - *gene-level counts* - start of the gene as defined in Ensembl biomart.
 - *exon-level counts* - center point of the exon.
 - *microarray probes* - start of the gene.
 -  *transcript usage* - start of the gene.
 - *txrevise* - start of the gene.
 - leafcutter - center point of the intron cluster.


<!--stackedit_data:
eyJoaXN0b3J5IjpbMTMyNzU1NTk3M119
-->