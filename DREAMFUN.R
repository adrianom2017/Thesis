
#With the information in KOexp and Cont, the functions returns a list for all experiments and corresponding 
#control in which KOgene was knocked out. (i.e. if gene10 has been KO in experiment 1 and 2, list[[1]] contains
#a data.frame with control data in the first x row(s) and experiment data in the y following row(s))
#Furthermore column 1 specifies the KO gene, columne 2 the experiment number, columne 3 shows the 
#control / experiment label
ContExpMerge = function(KOgene, KOexp, Cont, expression_data, rne, rnc){
  
  #find index of control experiments, find index of experiments apply indices to expression matrix
  idx = which(KOexp$DeletedGenes == KOgene) #find indices for all experiments in which gene x was perturbed
  exp = unique(KOexp$Experiment[idx]) #find all the experiments in which this gene was subject
  
  M = lapply(exp, function(x){
    idxc = rnc[which(Cont$Experiment == x)]
    idxe = rne[which(KOexp$Experiment == x & KOexp$DeletedGenes == KOgene)]
    tmp = as.data.frame(rbind(expression_data[idxc,], expression_data[idxe,]))
    
    experiment = x
    label = c(rep("Control", length(idxc)), rep("Knock_out", length(idxe)))
    
    tmp = cbind(KOgene, experiment, label, tmp)
  }
    )
}

#Returns a matrix with true positive, true negative, false positive, false negative, precision, recall (TPR) and FNR values at different threshold
#The threshold steps can be specified
#with the variable step which is set to 0.01 by default.
#mes interfered network
#ref gold standard network
PrecisionRecall = function(mes,ref, step = 0.01){
  thresholds = seq(0,1, step)
  c = sapply(thresholds, function(x){
    tmpmes = mes
    tmpmes[which(mes>=x)] = 1
    tmpmes[which(mes<x)] = 0
    c(length(which(ref & tmpmes)), length(which(!ref & !tmpmes)), length(which(!ref & tmpmes)), length(which(ref & !tmpmes))) #tp, tn, fp, fn
  })
  
  prec = apply(c, 2, function(x) x[1]/(x[1]+x[3]))
  recall = apply(c, 2, function (x) x[1]/(x[1]+x[4])) # == tpr
  fpr = apply(c, 2, function (x) x[3]/(x[3]+x[2]))
  
  tmp = as.data.frame(cbind(thresholds, t(c), prec, recall, fpr))
  colnames(tmp) = c("Threshold", "TP", "TN", "FP", "FN", "Precision", 'TPR', 'FPR')
  return(tmp)
}

#Returns an approximation of the area under the curve
AUC = function(x,y){
  deltax = sapply(2:(length(x)), function(u) x[u-1]-x[u])
  deltay = sapply(2:(length(y)), function(u) (y[u]+y[u-1])/2) #for area use mean height between points
  sum(abs(deltax)*abs(deltay))
}

Score = function(gold, pred, threshold, repl){
  print(paste0("Threshold: ", threshold))
  
  P = sum(gold, na.rm = T) #total number of positives in gold standard
  N = sum(gold == 0, na.rm = T) #total number of negatives in gold standard
  
  #Limit predictions to 100'000
  p = unlist(pred)
  g = unlist(gold)
  
  if (length(p) > 99999) {
    #DREAM challenge only scored the best 100'000 edges
    idx = order(p, decreasing = T)
    p = p[idx[1:99999]]
    g = g[idx[1:99999]]
  }
  
  #Find NA values in gold and remove them from both lists
  idx = which(is.na(gold))
  g = gold[-idx]
  p = pred[-idx]
  
  #Set all remaining NA values in the prediction matrix to repl
  idx = which(is.na(p))
  p[idx] = repl
  
  #tp, tn, fp, fn
  c = c(length(which(g & p)), length(which(!g & !p)), length(which(!g & p)), length(which(g & !p)))
  
  prec = c[1]/(c[1]+c[3])
  recall = c[1]/P # == tpr
  fpr = c[3]/N
  
  tmp = as.data.frame(cbind(threshold, t(c), prec, recall, fpr))
  colnames(tmp) = c("Threshold", "TP", "TN", "FP", "FN", "Precision", 'TPR', 'FPR')
  return(tmp)
}

bootstrap = function(files, path){
  #extract dimension
  print(paste0("Try to load file: ", path, files[1]))
  load(paste0(path,files[1]))
  d = dim(as(nem.res$graph, "matrix"))
  d = d[1]
  
  A = array(dim = c(d,d,1))
  for(ff in files){
    load(paste0(path, ff))
    A = abind(A, as(nem.res$graph, "matrix"), along = 3)
  }
  apply(A, 1:2, function(x) mean(x, na.rm = TRUE)) #remove NA values introduced by initialization of A
}

bootstrap_eg = function(files, path, dens = pdensity, hyp = hyper){
  
  #extract dimension
  print(paste0("Try to load file: ", path, files[1]))
  load(paste0(path,files[1]))
  d = dim(nem.res$pos)
  
  A = array(dim = c(d[1],d[2],1))
  for(ff in files){
    load(paste0(path, ff))
    
    g = as(nem.res$graph, "matrix")
    res = nem(dens, models = g, inference = "search", control = hyp)
    
    A = abind(A, res$pos[[1]], along = 3)
  }
  apply(A, 1:2, function(x) mean(x, na.rm = TRUE)) #remove NA values introduced by initialization of A
}


labelToName = function(gene, id){
  idx = sapply(gene, function(x) which(x == id$V1))
  return(id$V2[idx])
}

labelToEntrez = function(gene, id){
  idx = sapply(gene, function(x) which(x == id$V1))
  return(id$V3[idx])
}

nameToEntrez = function(gname, map){
  label = names(map)
  idx = sapply(gname, function(x) which(x == label))
  return(unlist(map[idx]))
}