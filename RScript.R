## >>>
### 2. Remove columns with fixed patterns in Prophane summary
### a. Function - Remove columns with "unique_spectra_of_"
## >>>

supp_colonnes_df <- function(my_data, ...){
  part1 <- select(my_data, -contains("unique.spectra.of."))
  part2 <- part1[, !names(part1) %in% c("considered.hits", "excluded.hits", 
                                        "maximal.sequence.length", "lowest.percent.identity")]
##  >>>
### b. Function - Remove "tr|" in the "hit.list" column
##  >>>
  tmp.hitlist = gsub("tr|", '', part2[,"hit.list"], fixed=TRUE)
  head(tmp.hitlist)
  part2[,"hit.list"] = tmp.hitlist
  return(part2)
}

## >>>
### d.  Merge Prophane summary and CAZyme file by hit.list column
## >>>

merge_func <- function(my_data1, my_data2, ...){
  part1 <- merge(my_data1, my_data2, by="hit.list",all.y=TRUE)
  return(part1)
}
