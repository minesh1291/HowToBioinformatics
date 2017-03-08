library("topGO")

makeGOResult= function(GOdata,WriteFile,myInterestedGenes){
  
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  result.E.Fish <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS tests")
  result.E.KS <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  result.C.Fish <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  result.C.KS <- getSigGroups(GOdata, test.stat)
  # whichAlgorithms()
  # whichTests()
  
  allGO = usedGO(object = GOdata) 
  # use it in GenTable as follows:
  
  allRes <- GenTable(GOdata, 
                     classicFisher = result.C.Fish, 
                     elimFisher = result.E.Fish,
                     classicKS = result.C.KS, 
                     elimKS = result.E.KS,
                     numChar=100,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(allGO)
  )
  
  class(allRes)#data.frame
  allRes$FoldEnrichment = allRes$Significant/allRes$Expected
  # sum(is.na(allRes$classicFisher))
  allRes$FDR.classicFisher = p.adjust(allRes$classicFisher,method = "BH")
  allRes$FDR.classicKS = p.adjust(allRes$classicKS,method = "BH")
  allRes$FDR.elimFisher = p.adjust(allRes$elimFisher,method = "BH")
  allRes$FDR.elimKS = p.adjust(allRes$elimKS,method = "BH")
  
  bg.geneList=genesInTerm(GOdata, whichGO=allRes$GO.ID)
  bg.geneList.v2 = lapply(bg.geneList,function(x) paste(x,collapse=", ") )
  
  int.geneList.v1 = lapply(bg.geneList,function(x) x[x %in% myInterestedGenes] )
  int.geneList.v2 = lapply(int.geneList.v1,FUN = function(x) paste(x,collapse=", ") )
  int.geneList.v2 =do.call(rbind.data.frame, int.geneList.v2)
  
  allRes$GeneSymbols=int.geneList.v2[allRes$GO.ID,1]
  
  write.csv(allRes,WriteFile)
}

####################


makeHs.udefined.BG.GOResult= function(allGenes,geneID2GO,myInterestedGenes,OntoType,WriteFile){
  #myInterestedGenes = HUGO_Symbols
  #OntoType = 'BP', CC, MF 
  #WriteFile = csv filename with path
  
  ## GO to Symbol mappings (only the BP ontology is used)
  xx <- annFUN.gene2GO(OntoType, gene2GO =  geneID2GO)
  # head(xx)
  allGenes_Hs <- unique(unlist(xx))
  # allGenes <- IAG_names[IAG_names %in% allGenes_Hs]
  allGenes <- allGenes[allGenes %in% allGenes_Hs]
  # myInterestedGenes <- sample(allGenes, 500)
  geneList <- factor(as.integer(allGenes %in% myInterestedGenes ))
  names(geneList) <- allGenes
  
  GOdata <- new("topGOdata",
                ontology = OntoType,
                allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO) 
  
  makeGOResult(GOdata,WriteFile,myInterestedGenes)
  
}

#######################
#Usage:

map.file="GOmapFile.txt"
geneID2GO <- readMappings(file = map.file)

geneNames <- names(geneID2GO)
allGenes=geneNames
top.Genes <- as.character(read.csv("topGenes.txt")[,1])

makeHs.udefined.BG.GOResult(allGenes = geneNames
                            ,myInterestedGenes = top.Genes
                            ,OntoType = "BP"
                            ,WriteFile = "top.Genes.BP.csv"
                            ,geneID2GO = geneID2GO)

