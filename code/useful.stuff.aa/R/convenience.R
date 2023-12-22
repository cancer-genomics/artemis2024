#https://stackoverflow.com/questions/23995384/read-and-rbind-multiple-csv-files
#' @export
squish <- function(file_list){
  dat <- 
  do.call(rbind,
          lapply(file_list, data.table::fread))

}


#' @export
get_split <- function(split_result,index){
  sapply( split_result, "[", index )


}

