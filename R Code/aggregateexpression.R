Idents(object.sct.filtered) <- 'predicted.celltype.l2'
Mono_cells <- subset(object.sct.filtered, idents = c("CD14 Mono", "CD16 Mono"))
table(Mono_cells$newer.ident)
#HC N10fs/N10fs 
#7344         297

pseudo_Mono <- AggregateExpression(Mono_cells, assays = "RNA", return.seurat = T, group.by = c("new.ident", "predicted.celltype.l2"))
Cells(pseudo_Mono)

pseudo_Mono$celltype.stim <- paste(pseudo_Mono$new.ident, pseudo_Mono$predicted.celltype.l2, sep = "_")

Idents(pseudo_Mono) <- "predicted.celltype.l2"

bulk.Mono.de <- FindAllMarkers(object = pseudo_Mono, 
                               test.use = "DESeq2")