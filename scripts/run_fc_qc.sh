Rscript feature_counts_qc.R\
 -c data/counts/Alasoo_test_data/Alasoo_merged_gene_counts.txt\
 -s data/sample_metadata/Alasoo_2018.tsv\
 -p data/annotations/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -o ./results/RNA_QC_script_results_Alasoo\
 -m ./mbv\
 --build_html TRUE