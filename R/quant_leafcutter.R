#Construct metadata df from leafcutter event ids
leafcutterConstructMeta <- function(leafcutter_ids){
  intron_metadata = dplyr::data_frame(phenotype_id = leafcutter_ids) %>%
    tidyr::separate(phenotype_id, c("chromosome","intron_start","intron_end","group_id"), sep = ":", remove = FALSE) %>%
    dplyr::mutate(intron_start = as.integer(intron_start), intron_end = as.integer(intron_end)) %>%
    dplyr::group_by(group_id) %>%
    dplyr::mutate(group_start = min(intron_start), group_end = max(intron_end)) %>%
    dplyr::ungroup()
  return(intron_metadata)
}

#' Annotate leafcutter introns with gene metadata
#'
#' @param leafcutter_ids vector of leafcutter intron ids
#' @param intron_annotation_path Path to leafcutter intron annotations file
#' @param transcript_meta transcript metadata data frame from importBiomartMetadata() function
#'
#' @return data frame of intron metadata
#' @export
leafcutterAnnotateIntrons <- function(leafcutter_ids, intron_annotation_path, transcript_meta){
  
  #Extract gene metadata from transcript metadata
  gene_meta = extractGeneMetadataFromBiomartFile(transcript_meta)
  
  #Import coordinates of gencode introns
  gencode_introns = readr::read_tsv(intron_annotation_path, 
                                    col_names = c("chromosome", "intron_start", "intron_end", "gene_name", "gene_id", 
                                                  "strand", "transcript_id", "intron_number","biotype", "tags"), 
                                    col_types = "ciiccccicc") %>%
    dplyr::select(chromosome, intron_start, intron_end, gene_name, gene_id, strand) %>%
    dplyr::filter(!(gene_id %like% "PAR_Y")) %>% #Remove PAR_Y genes
    dplyr::arrange(chromosome, intron_start, intron_end, gene_name) %>%
    dplyr::distinct() %>%
    tidyr::separate(gene_id, c("gene_id", "gene_version"), sep = "\\.")
  
  #Assign Leafcutter clusters to genes
  leafcutter_clusters = leafcutterConstructMeta(leafcutter_ids)
  match_introns = dplyr::left_join(leafcutter_clusters, gencode_introns, 
                                   by = c("chromosome", "intron_start", "intron_end")) %>%
    dplyr::select(group_id, gene_id, gene_name) %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::distinct()
  
  #Count genes
  gene_count = dplyr::group_by(match_introns, group_id) %>%
    dplyr::mutate(gene_count = length(group_id)) %>%
    dplyr::summarise(gene_id = paste(gene_id, collapse = ";"), 
                     gene_name = paste(gene_name, collapse = ";"), 
                     gene_count = gene_count[1]) %>%
    dplyr::ungroup()
  
  #Construct final metadata table
  #Gene metadata
  required_gene_meta_columns = c("phenotype_id","quant_id","group_id","gene_id",
                                 "chromosome","gene_start","gene_end","strand",
                                 "gene_name","gene_type","gene_gc_content",
                                 "gene_version","phenotype_pos")
  
  #Merge clusters with genes
  gene_features = dplyr::select(gene_meta, gene_id, gene_start, gene_end, gene_type, strand, 
                                gene_gc_content, gene_version, phenotype_pos)
  leafcutter_meta = dplyr::left_join(leafcutter_clusters, gene_count, by = "group_id") %>%
    dplyr::left_join(gene_features, by = "gene_id") %>%
    dplyr::mutate(gene_count = ifelse(is.na(gene_count), 0, gene_count)) %>%
    dplyr::rename(gene_pos = phenotype_pos) %>%
    dplyr::mutate(phenotype_pos = as.integer(ceiling((group_start + group_end)/2))) %>% #Use the center of the cluster as phenotype_pos
    dplyr::mutate(gene_type = ifelse(is.na(gene_type), "leafcutter", gene_type)) %>%
    dplyr::mutate(quant_id = group_id) %>%
    dplyr::select(required_gene_meta_columns, everything())
  
  return(leafcutter_meta)
}