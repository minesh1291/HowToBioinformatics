library("topGO")
mappingDB="org.Hs.eg.db"
library("org.Hs.eg.db")


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
  allRes$FDR.classicFisher = p.adjust(allRes$classicFisher,method = "fdr")
  allRes$FDR.classicKS = p.adjust(allRes$classicKS,method = "fdr")
  allRes$FDR.elimFisher = p.adjust(allRes$elimFisher,method = "fdr")
  allRes$FDR.elimKS = p.adjust(allRes$elimKS,method = "fdr")
  
  bg.geneList=genesInTerm(GOdata, whichGO=allRes$GO.ID)
  bg.geneList.v2 = lapply(bg.geneList,function(x) paste(x,collapse=", ") )
  
  int.geneList.v1 = lapply(bg.geneList,function(x) x[x %in% myInterestedGenes] )
  int.geneList.v2 = lapply(int.geneList.v1,FUN = function(x) paste(x,collapse=", ") )
  int.geneList.v2 =do.call(rbind.data.frame, int.geneList.v2)
  
  print(class(int.geneList.v2))
  print(head(int.geneList.v2))
  
  allRes$GeneSymbols=int.geneList.v2[,1] #[allRes$GO.ID,1]
  
  print(class(allRes$GeneSymbols))
  print(head(allRes$GeneSymbols))
  
  
  write.csv(allRes,WriteFile)
}




makeHs.GOResult= function(myInterestedGenes,OntoType,WriteFile){
  #myInterestedGenes = HUGO_Symbols
  #OntoType = BP, CC, MF 
  #WriteFile = csv filename with path
  
  ## GO to Symbol mappings (only the BP ontology is used)
  xx <- annFUN.org(OntoType, mapping = mappingDB, ID = "symbol")
  # head(xx)
  allGenes <- unique(unlist(xx))
  # myInterestedGenes <- sample(allGenes, 500)
  geneList <- factor(as.integer(allGenes %in% myInterestedGenes ))
  names(geneList) <- allGenes
  
  GOdata <- new("topGOdata",
                ontology = OntoType,
                allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.org, 
                mapping = mappingDB,
                ID = "symbol") 
  
  makeGOResult(GOdata,WriteFile,myInterestedGenes)
  
}

####################


makeHs.udefined.BG.GOResult= function(allGenes,myInterestedGenes,OntoType,WriteFile){
  #myInterestedGenes = HUGO_Symbols
  #OntoType = BP, CC, MF 
  #WriteFile = csv filename with path
  
  ## GO to Symbol mappings (only the BP ontology is used)
  xx <- annFUN.org(OntoType, mapping = mappingDB, ID = "symbol")
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
                annot = annFUN.org, 
                mapping = mappingDB,
                ID = "symbol") 
  
  makeGOResult(GOdata,WriteFile,myInterestedGenes)
  
}
