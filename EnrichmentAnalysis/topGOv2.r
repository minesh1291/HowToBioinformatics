# source("http://bioconductor.org/biocLite.R")
# biocLite()
# biocLite(c("topGO", "org.Hs.eg.db"))

library("topGO")
library("org.Hs.eg.db")
# save(makeGOResult,file = "makeGOResult.v1.RData")
makeGOResult= function(GOdata,WriteFile){
  
  test.stat <- new("elimCount", testStatistic = GOFisherTest, name = "Fisher test")
  result.E.Fish <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS tests")
  result.E.KS <- getSigGroups(GOdata, test.stat)
  
  test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  result.C.Fish <- getSigGroups(GOdata, test.stat)
   
  # pvalFis <- score(resultFis)
  
  test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")
  result.C.KS <- getSigGroups(GOdata, test.stat)
  
  
  # whichAlgorithms()
  # whichTests()
  
  # pvalWeight <- score(resultWeight, whichGO = names(pvalFis))
  
  allGO = usedGO(object = GOdata) 
  # use it in GenTable as follows:
  
  allRes <- GenTable(GOdata, 
                     classicFisher = result.C.Fish, 
                     elimFisher = result.E.Fish, 
                     classicKS = result.C.KS, 
                     elimKS = result.E.KS, 
                     orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(allGO)
  )
  
  write.csv(allRes,WriteFile)
  # biocLite("hgu95av2.db")
  # library("hgu95av2.db")
  # printGenes(object = GOdata,whichTerms = c("GO:0002320", "GO:0006820", "GO:0006821")
  #            ,chip="hgu95av2.db", file = "tmp.genes.csv")
  
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

