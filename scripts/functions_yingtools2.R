
#Functions from yingtools2
#------------------------
#I was not able to install package 'yingtools2' due to outdated environment,
#but these function seems to work. 


#' Extract Phyloseq tax_table
#'
#' Creates data.frame from tax_table, storing the rownames as variable "otu". The opposite of set.tax function.
#'
#' @param phy phyloseq object containing tax_data
#' @return Dataframe containing tax data
#' @export
get.tax <- function(phy) {
  requireNamespace(c("phyloseq"),quietly=TRUE)
  phyloseq::tax_table(phy) %>% data.frame(stringsAsFactors=FALSE) %>% rownames_to_column("otu") %>% as_tibble()
}



#' Convert data frame to phyloseq tax_table
#'
#' Use this on data.frames with tax data. The opposite of get.tax function. Make sure it contains the variable "otu".
#' @param tdata dataframe to be converted back to tax_table.
#' @return formatted tax_table.
#' @export
set.tax <- function(tdata) {
  requireNamespace(c("phyloseq"),quietly=TRUE)
  tt <- tdata %>% column_to_rownames("otu") %>%
    as.matrix() %>% phyloseq::tax_table()
  return(tt)
}



#' Convert Phyloseq to Melted OTU x Sample Data
#'
#' Creates OTU+Sample-level data, using phyloseq object (ID=otu+sample)
#'
#' Essentially gives back the OTU table, in melted form, such that each row represents a certain OTU for a certain sample.
#' Adds sample and taxonomy table data as columns. Uses the following reserved varnames: otu, sample, numseqs, pctseqs.
#' Note that phyloseq has a similar function, \code{psmelt}, but that takes longer.
#' The \code{get.otu.melt} now works by performing operations via data table, making it about 30x faster than before.
#'
#' @param phy phyloseq object containing sample data
#' @param filter.zero Logical, whether or not to remove zero abundances. Default \code{TRUE}.
#' @param tax_data Logical, whether or not to join with \code{tax_data}. Default \code{TRUE}.
#' @param sample_data Logical, whether or not to join with \code{sample_data}. Default \code{TRUE}.
#'
#' @return Data frame melted OTU data
#' @export
#' @examples
#' library(phyloseq)
#' get.otu.melt(cid.phy)
get.otu.melt <- function(phy,filter.zero=TRUE,sample_data=TRUE,tax_data=TRUE) {
  requireNamespace(c("phyloseq","data.table"),quietly=TRUE)
  # supports "naked" otu_table as `phy` input.
  otutab = as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) {
    otutab <- t(otutab)
  }
  otudt = data.table::data.table(otutab, keep.rownames = TRUE)
  data.table::setnames(otudt, "rn", "otu")
  # Enforce character otu key
  # note that .datatable.aware = TRUE needs to be set for this to work well.
  otudt[, otuchar:=as.character(otu)]
  otudt[, otu := NULL]
  data.table::setnames(otudt, "otuchar", "otu")
  # Melt count table
  mdt = data.table::melt.data.table(otudt, id.vars = "otu", variable.name = "sample", variable.factor=FALSE, value.name = "numseqs")
  if (filter.zero) {
    # Remove zeroes, NAs
    mdt <- mdt[numseqs > 0][!is.na(numseqs)]
  } else {
    mdt <- mdt[!is.na(numseqs)]
  }
  # Calculate relative abundance
  mdt[, pctseqs := numseqs / sum(numseqs), by = sample]
  if(tax_data & !is.null(phyloseq::tax_table(phy, errorIfNULL=FALSE))) {
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table::data.table(as(phyloseq::tax_table(phy, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    data.table::setnames(taxdt, "rn", "otu")
    # Enforce character otu key
    taxdt[, otuchar := as.character(otu)]
    taxdt[, otu := NULL]
    data.table::setnames(taxdt, "otuchar", "otu")
    # Join with tax table
    data.table::setkey(taxdt, "otu")
    data.table::setkey(mdt, "otu")
    mdt <- taxdt[mdt]
  }
  if (sample_data & !is.null(phyloseq::sample_data(phy, errorIfNULL = FALSE))) {
    # If there is a sample_data, join with it.
    sampledt = data.table::data.table(as(phyloseq::sample_data(phy, errorIfNULL = TRUE), "data.frame"),keep.rownames=TRUE)
    data.table::setnames(sampledt, "rn", "sample")
    # Enforce character sample key
    sampledt[, samplechar := as.character(sample)]
    sampledt[, sample := NULL]
    data.table::setnames(sampledt, "samplechar", "sample")
    # Join with tax table
    data.table::setkey(sampledt, "sample")
    data.table::setkey(mdt, "sample")
    mdt <- sampledt[mdt]
  }
  mdt <- mdt %>% as_tibble() %>% select(sample,otu,everything())
  return(mdt)
}



#' Log Epsilon Tranformation
#'
#' Use this transformation for plotting log data including 0. You can't use regular log transformation because it can't take zero.
#'
#' The transformation used is \eqn{\log{(|x|+\frac{epsilon}{8})} - \log(\frac{epsilon}{8})}, where epsilon is the parameter controlling the scale. The 1/8 portion is to make distances between ticks equal, so it's visually pleasing.
#' @param epsilon This parameter controls scaling. Think of this as the value of the first axis tick after zero. Default is 0.001.
#' @return Tranformation function to be plugged into ggplot.
#' @examples
#' values <- c(0,10^(-10:0))
#' d <- data.frame(x=1:length(values),y=values)
#' g <- ggplot(d,aes(x=x,y=y,label=y)) + geom_point() + geom_line() + geom_text()
#' g1 <- g + scale_y_continuous(breaks=values) + ggtitle("untransformed")
#' g2 <- g + scale_y_continuous(trans=log_epsilon_trans(0.0001)) + ggtitle("scale_trans, epsilon=0.0001")
#' g3 <- g + scale_y_continuous(trans=log_epsilon_trans(10^-6.)) + ggtitle("scale_trans, epsilon=0.0000001")
#' g4 <- g + scale_y_continuous(trans=log_epsilon_trans(10^-10)) + ggtitle("scale_trans, epsilon=0.0000000001")
#' gridExtra::grid.arrange(g1,g2,g3,g4,nrow=2)
#' @author Ying Taur
#' @export
log_epsilon_trans <- function(epsilon=0.001) {
  requireNamespace("scales",quietly=TRUE)
  trans <- function(x) {
    # if (is.null(epsilon)) {return(x)}
    sign(x)*(log(abs(x)+epsilon/8)-log(epsilon/8))
  }
  inv <- function(y) {
    # if (is.null(epsilon)) {return(y)}
    sign(y)*epsilon/8*(exp(abs(y))-1)
  }
  scales::trans_new(paste0("log_epsilon-",format(epsilon)),
                    transform = trans,
                    inverse = inv,
                    breaks=log_epsilon_trans_breaks(epsilon),
                    format=pretty_number,
                    domain=c(-Inf,Inf))
}

#' Breaks for Log Epsilon Tranformation
#'
#' This is used by scant_trans as default method for breaks. Will fill in logs of 10.
#' @param epsilon scaling parameter used in [log_epsilon_trans()]
#' @return break function returning break values.
#' @export
log_epsilon_trans_breaks <- function(epsilon) {
  # if (is.null(epsilon)) {return(scales::extended_breaks())}
  function(x) {
    firsttick <- round(log(epsilon,10))
    lasttick <- floor(log(x[2],10))
    x <- c(0,10^(firsttick:lasttick))
    by <- ceiling(length(x) / 5)
    x[seq(1,length(x),by=by)]
  }
}





#' Pretty Numeric Format (Non-scientific)
#'
#' Use to format axes in ggplot with non-scientific notation. Good for abundances!
#'
#' Note,
#' @param x number vector to be formatted.
#' @return Expression for x, in non-scientific notation.
#' @examples
#' x <- c(12,23.456789,1111e-7,230000022.11111,0.001234567)
#' pretty_number(x)
#'
#' dtime <- as.difftime(x,units="secs")
#' pretty_number(dtime)
#' @rdname pretty_number
#' @export
pretty_number <- function(x,...) UseMethod("pretty_number")
#' @rdname pretty_number
#' @export
pretty_number.default <- function(x,digits=2) {
  sapply(x,function(y) format(y,scientific=FALSE,trim=TRUE,big.mark=",",digits=digits))
}
#' @rdname pretty_number
#' @export
pretty_number.difftime <- function(x,...) {
  zero <- as.POSIXct(0,origin="1970-01-01")
  sapply(x,function(d) {
    diff <- d+zero-zero
    num <- diff %>% as.numeric() %>% pretty_number.default(...)
    units <- units(diff)
    paste(num,units)
  })
}




