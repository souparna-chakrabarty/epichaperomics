# Credits: Tai Wang

do_GSEA <- function (top_input, pvcutoff, gmt, group_name, Protein.ID.Column=NULL, minimumGSSize=5, maximumGSSize=10) {
  library (clusterProfiler)
  minimumGSSize = as.integer(minimumGSSize)
  maximumGSSize = as.integer(maximumGSSize)
  
  prot.id.col <- eval(parse(text=paste0("top_input$",Protein.ID.Column)))
  print("prot.id.col:")
  print(head(prot.id.col))
  genelist <- top_input [!prot.id.col=="" | is.na(prot.id.col),]
  gn.total <- eval(parse(text=paste0("genelist$",Protein.ID.Column)))
  genelist <- genelist[genelist$P.Value <= pvcutoff,]
  
  #DE.output.test1$top2 [DE.output.test1$top2$P.Value<0.1,] -> genelist ## for debug
  #genelist <- genelist [!(genelist$entrezgene_id=="" | is.na(genelist$entrezgene_id)),]
  #entrezgene.col <- "entrezgene_id" ## for debug
  #gn.total <- eval(parse(text=paste0("DE.output.test1$top2$",entrezgene.col)))
  
  genelist <- genelist[order (genelist$logFC, decreasing = T),]
  gl <- genelist$logFC
  
  #if (is.null(entrezgene.col)) {
  #  gn <- rowid2name (row.names(genelist), targetformat = "entrezgene",useonly1st = T)
  #  genelist <- data.frame (genelist, gn, stringsAsFactors = F)
  #  names (gl) <- gn
  #}
  
  #if (!is.null(entrezgene.col)) {
  gn <- eval(parse(text=paste0("genelist$",Protein.ID.Column)))
  #gn.total <- eval(parse(text=paste0("top_input$",entrezgene.col)))
  names (gl) <- gn
  #}
  
  ## Treat proteinGroups that has more than 1 Protein ID:
  ## This will add entry with multiple Protein IDs into GMT databases, the ont terms is based on pooled ont terms of the multiple EntrezGenes
  
  gene_split <- strsplit (names (gl), split = ";|\\|")
  pos_split <- sapply (gene_split, length) > 1
  p <- names (gl) [pos_split] # Identify proteingroups with multiple Protein IDs
  isinGMT <- any (gmt$gene %in% unlist (strsplit (p, split = ";|\\|")))
  
  if (any(pos_split) & isinGMT) {
    ont.sel <- sapply (p, function (x) { #For each PG, find the ont terms and pool terms, adding the multiple EntrezGenes IDs (e.g 654483;107282092;552900) as Gene name in GMT
      y <- unlist (strsplit(x, split=";|\\|"))
      print("about to set ont.sel")
      ont.sel <- unique (as.character (gmt [gmt$gene %in% y, "ont"]))
      return (ont.sel)
    }) 
    
    library (qdapTools)
    df <- list2df(ont.sel, col1="ont",col2="gene")
    #gmt <- gmt[!gmt$gene %in% unique (unlist (strsplit(df$gene, split=";"))),]
    gmt.2 <- rbind (gmt,df)
    gmt.2 <- gmt.2[order(gmt.2$ont),]
    #  gmt.3 <- gmt.2[!gmt.2$gene %in% unlist (strsplit(p, split=";")), ]
  } else {
    gmt.2 <- gmt
  }
  
  ## If there are duplicated gene names in the geneset, keep the one with the highest abs(t)
  order.abs.gl <- order (abs (gl), decreasing = T)
  pos_dup <-  gl[order.abs.gl] %>% names %>% duplicated
  val.gl.dup <- gl[order.abs.gl] [pos_dup]
  gl <- gl [!gl %in% val.gl.dup] ## This assumes the values are unique to each gene
  name.gl.dup <- names (gl[order.abs.gl] [pos_dup])
  
  #gmt.3 <- gmt.2 [!grepl (";",gmt.2$gene),]
  # droplevels(gmt.3) -> gmt.3
  
  gmt.2 <- gmt.2 [gmt.2$gene %in% gn.total,] ###!!!!!!!!!!!

  #print("************************ minGSSize")
  #print(minimumGSSize)
  #print("************************ maxGSSize")
  #print(maximumGSSize)
  #print("************************ typeof(minGSSize)")
  #print(typeof(minimumGSSize))
  #print("************************ typeof(maxGSSize)")
  #print(typeof(maximumGSSize))

  em2 <- clusterProfiler::GSEA(gl, TERM2GENE = gmt.2, pAdjustMethod = "BH", minGSSize = minimumGSSize, maxGSSize = maximumGSSize, seed = TRUE, pvalueCutoff = 1)
  
  #print("************************ min(em2$setSize)")
  #print(min(em2$setSize))
  #print("************************ max(em2$setSize)")
  #print(max(em2$setSize))
  #print("************************ em2 %>% group_by(setSize)%>% summarise(total_count=n(), .groups = 'drop')")
  #print(em2 %>% group_by(setSize)%>% summarise(total_count=n(), .groups = 'drop'))

  #em2 <- enricher(names (gl), TERM2GENE = gmt.3, universe = gn.total, pAdjustMethod = "BH",minGSSize = 1, maxGSSize = 9999,  pvalueCutoff = 1)
  #em2.readable <- setReadable(em2, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
 
  #y <- as.data.frame (em2.readable)
  y <- as.data.frame (em2)

  y$Group <- NA
  y$Group <- group_name
  
  ce <- strsplit(y$core_enrichment, split ="\\/")
  
  # Convert UniProt ID from core_enrichment into Shown.ID and pass into 'core_enrichment_shown_id' column
  prot.id.col <- eval(parse(text=paste0("top_input$",Protein.ID.Column)))
  
  core_enrichment_shown_id <- sapply (ce, function (x){
    pos <- match (x, prot.id.col)
    gs <- top_input$Shown.ID[pos]
    if (length (gs)==0) {gs==""}
    return (paste0(gs, collapse = "/"))
  })
  
  y$core_enrichment_shown_id <- core_enrichment_shown_id
  print("----------------------------> y$setSize")
  print(y$setSize)
  #save(y, file = "/Users/wilson/Work/Chiosis - Chaperomics/y.RData")
  
  return (list(em=em2, genelist=genelist, gl=gl, df=y))
}
