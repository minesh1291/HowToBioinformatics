makeGOResult= function(GOdata,WriteFile){
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  pvalFis <- score(resultFis)
  
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  resultKS <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "Fisher test", cutOff = 0.01)
  resultElim <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  
  pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
  
  allGO = usedGO(object = GOdata) 
  # use it in GenTable as follows:
  
  allRes <- GenTable(GOdata, classicFisher = resultFis, elimFisher = resultElim, KS = resultKS, weight = resultWeight,
                     orderBy = "weight", ranksOf = "classicFisher", topNodes = length(allGO)
  )
  
  write.csv(allRes,WriteFile)
  
}

##################################################
## Using the org.XX.xx.db annotations
##################################################

makeHs.GOResult= function(myInterestedGenes,OntoType,WriteFile){
  #myInterestedGenes = HUGO_Symbols
  #OntoType = BP, CC, MF 
  #WriteFile = csv filename with path
  
  ## GO to Symbol mappings (only the BP ontology is used)
  xx <- annFUN.org(OntoType, mapping = "org.Hs.eg.db", ID = "symbol")
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
                mapping = "org.Hs.eg.db",
                ID = "symbol") 
  
  makeGOResult(GOdata,WriteFile)
  
}

####################


makeHs.udefined.BG.GOResult= function(allGenes,myInterestedGenes,OntoType,WriteFile){
  #myInterestedGenes = HUGO_Symbols
  #OntoType = BP, CC, MF 
  #WriteFile = csv filename with path
  
  ## GO to Symbol mappings (only the BP ontology is used)
  xx <- annFUN.org(OntoType, mapping = "org.Hs.eg.db", ID = "symbol")
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
                mapping = "org.Hs.eg.db",
                ID = "symbol") 
  
  makeGOResult(GOdata,WriteFile)
  
}

