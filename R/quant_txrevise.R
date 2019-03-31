constructTxreviseRowData <- function(phenotype_ids, transcript_meta){
  
  #Split phenotype ids into components
  event_metadata = dplyr::data_frame(phenotype_id = phenotype_ids) %>%
    tidyr::separate(phenotype_id, c("gene_id", "txrevise_grp", "txrevise_pos", "transcript_id"), sep = "\\.", remove = FALSE) %>%
    dplyr::mutate(group_id = paste(gene_id, txrevise_pos, sep = "."), quant_id = paste(gene_id, txrevise_grp, txrevise_pos, sep = ".")) %>%
    dplyr::select(phenotype_id, quant_id, group_id, gene_id)
  
  #Extract gene metadata
  gene_meta = extractGeneMetadataFromBiomartFile(transcript_meta)
  
  row_data = dplyr::left_join(event_metadata, gene_data, by = "gene_id")
  return(row_data)
}
