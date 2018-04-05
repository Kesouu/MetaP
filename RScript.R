## >>>
### Part 2 : Data manipulation
### Remove columns or patterns in Prophane summary
### a. Function - Remove columns with "unique_spectra_of_" and remove "tr|" in the "hit.list" column
## >>>

supp_colonnes_df <- function(my_data, ...){
  part1 <- select(my_data, -contains("unique.spectra.of."),-one_of(c("considered.hits", "excluded.hits","maximal.sequence.length", "lowest.percent.identity"))) 
  tmp.hitlist = gsub("tr|", '', part1[,"hit.list"], fixed=TRUE)
  head(tmp.hitlist)
  part1[,"hit.list"] = tmp.hitlist
  return(part1)
}

## >>>
### b.  Merge Prophane summary and CAZyme file by hit.list column
## >>>

merge_func <- function(my_data1, my_data2, ...){
  part1 <- merge(my_data1, my_data2, by="hit.list",all.y=TRUE)
  return(part1)
}

## >>>
### Part 3 : Replicate
### b. Function in case of 3 replicates - if 2 values are available else keep values, or paste 0 if <2 values are available. Renormalize the NSAF value between [0-1]
### Then, mean 3 replicates (or 2 if 3rd value is 0)
## >>>

Norm3Replicate <- function(my_data,  pattern, colname,...){
  part1 <- select(my_data,matches(pattern))
  part2 <- ifelse((apply(part1==0,1,sum)>=2),0,part1[,1])
  part3 <- ifelse((apply(part1==0,1,sum)>=2),0,part1[,2])
  part4 <- ifelse((apply(part1==0,1,sum)>=2),0,part1[,3])
  part5 <- cbind(A=part2,B=part3,C=part4)
  part5 <- data.frame(part5)
  part6 <- apply(part5,2,sum)
  part6 <- data.frame(part6)
  part7 <- part5$A/part6[1,]
  part8 <- part5$B/part6[2,]
  part9 <- part5$C/part6[3,]
  part10 <- as.data.frame(cbind(part7,part8,part9))
  part10 <- replacer(part10, 0,NA)
  part11 <- apply(part10,1,mean,na.rm=T)
  part11 <- data.frame(part11)
  colnames(part11) <- c(colname)
  Merge_summary_final_part1 <- cbind.data.frame(Merge_summary_final_part1,part11)
  return(Merge_summary_final_part1)
}

### c. Function in case of 2 replicates - if 2 values are available else keep values, or paste 0 if <2 values are available. Renormalize the NSAF value between [0-1]
### Then, mean 3 replicates (or 2 if 3rd value is 0)
## >>>

Norm2Replicate <- function(my_data,  pattern, colname,...){
  part1 <- select(my_data,matches(pattern))
  part1 <- data.frame(part1)
  #part1 <- replacer(part1, 0,NA)
  part2 <- apply(part1,1,mean) #na.rm=T
  part2 <- data.frame(part2)
  colnames(part2) <- c(colname)
  Merge_summary_final_part1 <- cbind.data.frame(Merge_summary_final1,part2)
  return(Merge_summary_final_part1)
}
