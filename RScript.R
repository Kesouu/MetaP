## >>>
### 2. Remove columns or patterns in Prophane summary
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
### d.  Merge Prophane summary and CAZyme file by hit.list column
## >>>

merge_func <- function(my_data1, my_data2, ...){
  part1 <- merge(my_data1, my_data2, by="hit.list",all.y=TRUE)
  return(part1)
}
