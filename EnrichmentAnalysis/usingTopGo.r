library(topGO)

geneID2GO <- readMappings(file = system.file("examples/RiceMapFunc", package = "topGO")) # loading background GO information from special formated file
geneNames <- names(geneID2GO) # give gene names of your interest
myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList)= geneNames
GOdata <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO)
GOdata

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
 resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
 resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
 classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher",topNodes=40)
write.table(allRes,"RiceFunc.txt")
