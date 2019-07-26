library(PADOG)
library(limma)
library(Biobase)
getdataaslist = function(x) {
  x = get(x, envir=parent.frame())
  exp = experimentData(x)
  dataset = exp@name
  disease = notes(exp)$disease
  dat.m = exprs(x)
  ano = pData(x)
  design = notes(exp)$design
  annotation = paste(x@annotation, ".db", sep = "")
  targetGeneSets = notes(exp)$targetGeneSets
  list = list(dataset, disease, dat.m, ano, design, annotation, targetGeneSets)
  names(list) = c("dataset", "disease", "dat.m", "ano", "design", "annotation", 
                  "targetGeneSets")
  return(list)
}
files = data(package = "KEGGdzPathwaysGEO")$results[, "Item"]
data(list = files, package = "KEGGdzPathwaysGEO", envir=environment())
nG = array(dim=length(files))
nDe = array(dim=length(files))
k=1
for(i in 1:length(files)){
  file = files[i]
  list = getdataaslist(file)
  esetm = list$dat.m
  group = list$ano$Group
  paired = list$design
  if(paired == "Not Paired"){
    paired = FALSE
  }else{
    paired = TRUE
  }
  block = list$ano$Block
  annotation = list$annotation
  
  if (!is.null(annotation)) {
    aT1 = filteranot(esetm, group, paired, block, annotation)
    esetm = esetm[rownames(esetm) %in% aT1$ID, ]
    rownames(esetm) <- aT1$ENTREZID[match(rownames(esetm), 
                                          aT1$ID)]
  }
  geneNames = rownames(esetm)
  G = factor(group)
  if (paired) {
    design <- model.matrix(~0 + G + block)
    colnames(design) <- substr(colnames(design), 2, 100)
  }else{
    design <- model.matrix(~0 + G)
    colnames(design) <- levels(G)
  }
  fit <- lmFit(esetm, design)
  cont.matrix <- makeContrasts(contrasts = "d-c", levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  pValue=fit2[["p.value"]]
  pValue = array(pValue)
  pValue = p.adjust(pValue,"fdr")
  difE = array()
  nG = array()
  nDe = array()
  j=1
  for(i in 1:length(pValue))
    if(pValue[i]<0.01)
    {
      difE[j]=geneNames[i]
      j=j+1
    }
  nG[k]=length(geneNames)
  nDe[k]=j
  k = k+1
  cat(length(geneNames),file = paste('D:/NIR/bernmix_method/datasets/all.txt'),sep = "\n", append = TRUE)
  cat(j,file = paste('D:/NIR/bernmix_method/datasets/de.txt'),sep = "\n", append = TRUE)
  cat(geneNames, file = paste('D:/NIR/bernmix_method/datasets/',file,'all.txt'))
  cat(difE, file = paste('D:/NIR/bernmix_method/datasets/',file,'_de.txt'))
}