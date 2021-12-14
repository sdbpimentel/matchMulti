#' Create Balance Table
#'
#' Given an unmatched sample of treated and control units and (optionally) a
#' matched sample from the same data, produces a table with pre- and post-match
#' measures of covariate balance.
#'
#'
#' @param df.orig a data frame containing the data before matching
#' @param df.match an optional data frame containing the matched sample. Must
#'   have all variable names to be balanced.
#' @param treatment name of the binary indicator for treatment status
#' @param var.names List of variable names to calculate balance for.  If NULL,
#'   use all variables found in the df.orig data.frame.
#' @param treat.wts optional weights for treated units in the original sample
#' @param ctrl.wts optional weights for control units in the original sample
#' @param mt.wts optional weights for treated units in the matched sample
#' @param mc.wts optional weights for treated units in the matched sample
#' @param verbose a logical value indicating whether detailed output should be
#'   printed.
#' @param include.tests Include tests of imbalance on covariates (TRUE/FALSE).
#'
#' @return A data.frame of balance measures, with one row for each covariate in
#'   \code{df.orig} except \code{treatment}, and columns for treated and control
#'   means, standardized differences in means, p-values from a 2-sample t-test,
#'   and p-values from either Fisher's exact test (if the covariate is binary)
#'   or a Wilcoxon signed rank test otherwise.  If \code{df.match} is specified
#'   there are twice as many columns, one set for the pre-match samples and one
#'   set for the post-match samples.
#'
#' @importFrom plyr aaply
#'
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#'   Springer-Verlag.
#'
#'   Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#'   Springer-Verlag.
#' @keywords internal
#' @export balanceTable
balanceTable <- function(df.orig, df.match = NULL,
                         treatment, 
                         var.names = NULL,
                         #cat.vars = NULL,
                         treat.wts = NULL, ctrl.wts = NULL, mt.wts = NULL, mc.wts = NULL, 
                         include.tests = FALSE,
                         verbose = FALSE){
  #if(is.null(cat.vars)) cat.vars <- rep(FALSE, ncol(df.orig))
  
  #Validate input
  stopifnot( is.data.frame(df.orig) )
  stopifnot( is.null(df.match) || is.data.frame(df.match) )
  
  if (!(treatment %in% colnames(df.orig))) {
    stop(paste0('Treatment variable "', treatment,'" not found'))
  }	
  
  if ( !is.null( var.names ) ) {
    stopifnot( is.character(var.names) )
    if ( length( var.names ) == 0 ) {
      stop( "Cannot have no variables to balance on. Set var.names to NULL to calc balance on everything." )
    }
    stopifnot( all( var.names %in% names(df.orig ) ) )
    df.orig = df.orig[ c(treatment, var.names) ]
  }
  
  
  
  if(!is.null(df.match) ) {
    if ( ! all( colnames(df.orig) %in% colnames(df.match) ) ) {
      stop( 'df.match must have tx column and all columns selected for balance in df.orig')
      df.match = df.match[ colnames(df.orig) ]
    }
  }
  
  
  #if(!is.null(cat.vars) && length(cat.vars) != ncol(df.orig)){
  #	stop('cat.vars must have exactly one entry for each column in df.orig')
  #}
  if(any(is.na(df.orig[[treatment]])) || (!is.null(df.match) && any(is.na(df.match[[treatment]])))){
    stop('NAs are present in the treatment variable')
  }
  
  #non.numeric <- colnames(df.orig)[laply(df.orig, inherits, what = c('character','factor'))]
  
  cov.orig <- handleNA(df.orig, verbose = verbose)
  if(!is.null(df.match)){		
    cov.match <- handleNA(df.match, verbose = verbose)
    cov.match <- resolve.cols(cov.match, cov.orig)
    cov.orig <- resolve.cols(cov.orig, cov.match)
  } else {
    cov.match <- NULL
  }
  
  #search through and find all binaries.
  binary.ind <- laply(cov.orig, is.binary)
  treat.ind <- colnames(cov.orig) == treatment
  
  sdiff.out <- plyr::aaply(colnames(cov.orig)[which(!treat.ind)], 1, sdiff, 
                           treatment = treatment, orig.data = cov.orig, match.data = cov.match, 
                           treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                           mt.wts = mt.wts, mc.wts = mc.wts, .drop = FALSE)
  rownames(sdiff.out) <- colnames(cov.orig)[which(!treat.ind)]
  
  #TODO: figure out whether to keep weight arguments and incorporate them or drop them
  
  if ( include.tests ) {
    
    t.test.out <- plyr::aaply(colnames(cov.orig)[which(!treat.ind)], 1, 
                              ttest.balance, 
                              treatment = treatment, orig.data = cov.orig, match.data = cov.match, 
                              treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                              mt.wts = mt.wts, mc.wts = mc.wts, .drop = FALSE)
    rownames(t.test.out) <-  colnames(cov.orig)[which(!treat.ind)]
    
    
    if( any(!binary.ind)){
      wilc.test.out <- plyr::aaply(colnames(cov.orig)[which(!binary.ind & !treat.ind)], 1, 
                                   wilc.balance,
                                   treatment = treatment, orig.data = cov.orig, match.data = cov.match, 
                                   treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                                   mt.wts = mt.wts, mc.wts = mc.wts, .drop = FALSE)
      rownames(wilc.test.out) <- colnames(cov.orig)[which(!binary.ind & !treat.ind)]
    }
    
    if( any(binary.ind[-which(treat.ind)])){
      fisher.test.out <- plyr::aaply(colnames(cov.orig)[which(binary.ind & !treat.ind)], 1, 
                                     fisher.balance, 
                                     treatment = treatment, 
                                     orig.data = cov.orig, match.data = cov.match, 
                                     treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                                     mt.wts = mt.wts, mc.wts = mc.wts, .drop = FALSE) 
      rownames(fisher.test.out) <- colnames(cov.orig)[which(binary.ind & !treat.ind)]
    }
    
    if (all(binary.ind)){
      test.out <- fisher.test.out
    } else if (all(!binary.ind[-which(treat.ind)])) {
      test.out <- wilc.test.out	
    } else {
      test.out <- as.data.frame(matrix(nrow = nrow(sdiff.out), ncol = ncol(fisher.test.out)))
      rownames(test.out) <- rownames(sdiff.out)
      if (ncol(test.out) == 2){
        colnames(test.out) <- c('Fisher/Wilcox Pvalue Before', 'Fisher/Wilcox Pvalue After')
      } else {
        colnames(test.out) <- 'Fisher/Wilcox Pvalue'
      }
      test.out[which(rownames(sdiff.out) %in% rownames(fisher.test.out)),] <- fisher.test.out
      test.out[which(rownames(sdiff.out) %in% rownames(wilc.test.out)),] <- wilc.test.out
    }
    
    out.tab <- cbind(sdiff.out,t.test.out, test.out)
    as.data.frame( out.tab )
  } else {
    as.data.frame( sdiff.out )
  }
}
