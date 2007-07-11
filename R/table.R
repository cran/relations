relation_table <-
function(x)
{
    if(!is.relation(x))
        stop("Argument 'x' must be a relation.")
    structure(as.data.frame(x),
              class = c("relation_table", "data.frame"))
}

print.relation_table <-
function(x, ...)
{
    y <- as.data.frame(x)
    if (length(row.names(y)) == 0L)
      print.data.frame(y)
    else {
      y <- as.matrix(format(y, justify = "left"))
      rownames(y) <- rep.int("", nrow(y))
      print(y, quote = FALSE)
    }
    invisible(x)
   
}
