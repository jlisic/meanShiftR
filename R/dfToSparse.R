#' Convert data frame to sparse matrix
#' 
#' \code{dfToSparse} performs one-hot encoding to characters and factor types
#' within a data frame.  The encoding is returned as a sparse matrix.
#'
#' @param x A data frame to convert into a sparse matrix.  Only character and
#'   factor types are converted through one-hot encoding.
#' @param transpose A boolean determining if you should query by columns instead of
#'   the default of rows (only for sparse matricies).
#' @return A sparse matrix with one-hot encoded columns, numeric and integer
#'   columsn are preserved in this sparse matrix.  Column names are returned.
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' neighbors <- dfToSparse(x_df) 
#' @useDynLib meanShiftR, .registration = TRUE
#' @export
dfToSparse <-
function(
  x,
  transpose = FALSE
) {

  max_cat <- apply(a,2,max)
  max_cat_cumsum <- c(0, cumsum(max_cat))
  n <- NROW(a)
  p <- NCOL(a)
  
  # transpose
  new_i <- rep(0,n*p)
  new_p <- cumsum(rep(p,n)) 
  
#  for( i in 1:n ) {
#    a_row <- a[i,]
#  
#    for( j in 1:p ) {
#      if( a_row[j] != 0 ) {
#        new_i = c(new_i, max_cat_cumsum[j] + (a_row[j] -1) )
#      }
#    }
#  }

  type <- c()


  r.result <- .C("R_dfToSparse",
    as.integer(x ),            # 1 Not used yet  
    as.integer( n ),          # 2 number of rows 
    as.integer( p ),          # 3 number of columns 
    as.integer( new_i ),      # 4 column indexes
    as.integer( new_p ),      # 5 pointers for sparse matrix 
    as.integer( max_cat ),    # 6 number of categories per column
    as.integer( max_cat_cumsum ), # cumulative sume of the number of categories per column
    as.integer( type ), # 7 Not used yet
    NAOK=FALSE
  )
  
  # create a sparse matrix
  out_sparse <- Matrix(1.1,nrow=2,ncol=1,sparse=TRUE)
  out_sparse@i <- r.result[[4]]
  out_sparse@p <- as.integer(c(0,new_p))
  out_sparse@x <- rep(1,n*p)
  out_sparse@Dim <- as.integer(c(sum(max_cat),n))

  return(out_sparse)
}
