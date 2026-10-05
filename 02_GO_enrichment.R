library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
# GO enrichment
clustergo <- function(gene, universe, ont = "BP",
                      pvalue = 0.05, qvalue = 0.05, minGSSize = 3) {
  
  gene_id <- suppressMessages(
    bitr(gene, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  )
  
  universe_id <- suppressMessages(
    bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)
  )
  
  enrichGO(
    gene = unique(gene_id$ENTREZID),
    universe = unique(universe_id$ENTREZID),
    OrgDb = org.Hs.eg.db,
    ont = ont,
    readable = TRUE,
    pvalueCutoff = pvalue,
    qvalueCutoff = qvalue,
    minGSSize = minGSSize
  ) |>
    as.data.frame()
}



background <- fread(
  file.path(data_dir, "Report.tsv")
)

background_genes <- unique(background$PG)



# GO enrichment
go_bp <- clustergo(
  gene = protein, universe = background_genes,
  ont = "BP", pvalue = 1, qvalue = 1
)

go_cc <- clustergo(
  gene = protein, universe = background_genes,
  ont = "CC", pvalue = 1, qvalue = 1
)

go_mf <- clustergo(
  gene = protein, universe = background_genes,
  ont = "MF", pvalue = 1, qvalue = 1
)