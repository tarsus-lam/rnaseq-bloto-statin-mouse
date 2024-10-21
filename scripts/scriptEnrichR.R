# Load necessary library
library(biomaRt)

# Connect to the Ensembl database using the biomaRt package
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Convert Ensembl gene IDs to Entrez Gene IDs
ensembl_to_entrez <- getLDS(attributes = c("ensembl_gene_id_version", "entrezgene_id"), mart = ensembl)

# Load your gene list (replace "gene_list.txt" with your file path)
gene_list <- readLines("DifferentialExpression/gene_list.txt")

# Query the conversion table for your gene list
conversion_result <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), 
                           filters = "ensembl_gene_id", 
                           values = gene_list, 
                           mart = ensembl)

# Print the conversion result
print(conversion_result)

write.table(conversion_result['entrezgene_id'], file = "entrez_gene_list.txt", quote = FALSE, col.names = FALSE, row.names = FALSE)
