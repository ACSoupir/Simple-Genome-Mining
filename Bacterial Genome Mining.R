rm(list=ls())
getwd()

system.time({
  
  
  #load packages for importing fasta
  #library(seqinr)
  library(Biostrings)
  
  #find fasta files in the 
  anno1 = list.files(path="Annotations/", pattern="*.fasta")
  
  anno2 = as.character()
  for(i in 1:length(anno1)){
    anno2[i] = paste("Annotations/", anno1[i], sep="")
  }
  
  seqCount = data.frame(matrix(ncol = 3, nrow = length(anno1)))
  colnames(seqCount) = c("TotalNumber", "HypoRemoved", "NoDuplicate")
  row.names(seqCount) = anno1
  
  anno2.list = list()
  for(i in 1:length(anno2)){
    s = readAAStringSet(anno2[i])
    
    seq_names = names(s)
    rexp <- "^(\\w+)\\s?(.*)$"
    t = data.frame(SeqID=sub(rexp,"\\1",seq_names), ProtName=sub(rexp,"\\2",seq_names))
    
    sequence = paste(s)
    df = data.frame(t[2], sequence)
    #print(nrow(df))
    
    anno2.list[[i]] = df
    names(anno2.list)[i] <- paste(anno1[i])
    
    seqCount[i,1] = nrow(df)
  }
  
  anno3.list = list()
  for(i in 1:length(anno2.list)){
    #remove hypothetical proteins so only identifiable ones are compared
    realProt = subset(anno2.list[[i]], ProtName != "hypothetical protein")
    
    realProt = realProt[order(realProt$ProtName),]
    
    #print(i)
    #print(nrow(realProt))
    seqCount[i,2] = nrow(realProt)
    
    #removes duplicated sequences
    realProt = realProt[!duplicated(realProt[,1]),]
    
    #print(nrow(realProt))
    seqCount[i,3] = nrow(realProt)
    
    anno3.list[[i]] = realProt
    names(anno3.list)[i] <- paste(anno1[i])
  }
  
  allSamples = anno3.list[[i]]
  colnames(allSamples)[2] = anno1[1]
  for(i in 2:length(anno3.list)){
    allSamples = merge(allSamples, anno3.list[[i]], by = "ProtName", all = TRUE)
    colnames(allSamples)[i+1] = anno1[i]
  }
  
  #sort gene matrix by gene name
  row.names(allSamples) = allSamples$ProtName
  allSamples = allSamples[,-1]
  allSamples = allSamples[order(row.names(allSamples)),]
  
  #write gene matrix to csv
  #write.csv(allSamples, file = "Clean and Merged Proteins.csv")
})

charKey = read.csv("Characteristic Key.csv")
charKey = charKey[order(charKey$Genome),]

#chars = as.vector(charKey$IAA)
#chars = as.vector(charKey$ACCd)
#chars = as.vector(charKey$F.graminearum)
#chars = as.vector(charKey$F.accuminatum)
#chars = as.vector(charKey$F.oxysporum)
chars = as.vector(charKey$F.proliferatum)

system.time({
  rm(y)
  rm(n)
  for(i in 1:length(allSamples)){
    if(is.na(chars[i])){
      next
    }
    if(chars[i] == 1){
      if(exists("y") == 0){
        y = as.data.frame(allSamples[,i])
        row.names(y) = row.names(allSamples)
        colnames(y)[ncol(y)] = colnames(allSamples)[i]
      } else {
        y[ncol(y)+1] = allSamples[,i]
        colnames(y)[ncol(y)] = colnames(allSamples)[i]
      }
    }
    
    if(chars[i] == 0){
      if(exists("n") == 0){
        n = as.data.frame(allSamples[,i])
        row.names(n) = row.names(allSamples)
        colnames(n)[ncol(n)] = colnames(allSamples)[i]
      } else {
        n[ncol(n)+1] = allSamples[,i]
        colnames(n)[ncol(n)] = colnames(allSamples)[i]
      }
    }
  }
  
  y = y[complete.cases(y),]
  n = n[complete.cases(n),]
  
  dTable = merge(y, n, by = 0, all = TRUE)
  
  dTable = dTable[!complete.cases(dTable),]
  row.names(dTable) = dTable$Row.names
  dTable = dTable[order(row.names(dTable)),]
  dTable = dTable[,-1]
  
})

#write.csv(dTable, file = "Different Genes in IAA Isolates.csv")
#write.csv(dTable, file = "Different Genes in ACCd Isolates.csv")
#write.csv(dTable, file = "Different Genes in F graminearum suppressing Isolates.csv")
#write.csv(dTable, file = "Different Genes in F accuminatum suppressing Isolates.csv")
#write.csv(dTable, file = "Different Genes in F oxysporum suppressing Isolates.csv")
write.csv(dTable, file = "Different Genes in F proliferatum suppressing Isolates.csv")








