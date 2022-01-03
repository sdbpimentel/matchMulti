

balanceTable.internal <- function(data,
                                  treatment, 
                                  school.id = NULL,
                                  var.names = NULL,
                                  treat.wts = NULL, ctrl.wts = NULL, 
                                  include.tests = FALSE,
                                  verbose = FALSE) {
  
  vars = setdiff(colnames(data), treatment )
  if ( !is.null( school.id ) ) {
    vars = setdiff( vars, school.id )
  }
  
  sdiff.out <- plyr::laply(vars, sdiff, 
                           treatment = treatment, 
                           orig.data = data, 
                           treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                           .drop = FALSE)
  rownames(sdiff.out) <- vars
  
  
  #TODO: figure out whether to keep weight arguments and incorporate them or
  #drop them
  
  if ( include.tests ) {
    
    t.test.out <- plyr::laply(vars, agg.balance, 
                              treatment = treatment, 
                              school.id = school.id,
                              data = data, 
                              treat.wts = treat.wts, ctrl.wts = ctrl.wts,
                              .drop = FALSE)
    rownames(t.test.out) <-  vars
    colnames(t.test.out) <- "Agg PValue"
    
    CRVE.out <- plyr::laply(vars, CRVE.balance, 
                            treatment = treatment, 
                            school.id = school.id,
                            data = data, 
                            treat.wts = treat.wts, ctrl.wts = ctrl.wts,
                            .drop = FALSE)
    rownames(CRVE.out) <-  vars
    colnames(CRVE.out) <- "CRVE PValue"
    
    out.tab <- cbind(sdiff.out, t.test.out, CRVE.out)
    as.data.frame( out.tab )
  } else {
    as.data.frame( sdiff.out )
  }
}





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
#' @param school.id Identifier for groups (for example schools); need to pass if
#'   p-values for balance statistics are desired.
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
                         school.id = NULL,
                         var.names = NULL,
                         #cat.vars = NULL,
                         treat.wts = NULL, ctrl.wts = NULL, 
                         mt.wts = NULL, mc.wts = NULL, 
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
    if ( !is.null( school.id ) ) {
      df.orig = df.orig[ c(treatment, school.id, var.names) ]
    } else {
      df.orig = df.orig[ c(treatment, var.names) ]
    }
  }
  
  
  
  if(!is.null(df.match) ) {
    if ( ! all( colnames(df.orig) %in% colnames(df.match) ) ) {
      stop( 'df.match must have tx column and all columns selected for balance in df.orig')
    }
    df.match = df.match[ colnames(df.orig) ]
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
  
  main_table = balanceTable.internal(data = cov.orig,
                                     treatment = treatment, school.id = school.id, var.names = var.names, 
                                     treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                                     include.tests = include.tests, verbose = verbose )
  
  # Add in match balance if present.
  if ( !is.null( df.match ) ) {
    match_table = balanceTable.internal(data = cov.match,
                                     treatment = treatment, school.id = school.id, var.names = var.names, 
                                     treat.wts = treat.wts, ctrl.wts = ctrl.wts, 
                                     include.tests = include.tests, verbose = verbose )
    
    colnames(main_table) = paste0( "Before ", colnames(main_table) )
    colnames(match_table) = paste0( "After ", colnames(match_table) )
    main_table = cbind( main_table, match_table )
  }

  return( main_table )
}


