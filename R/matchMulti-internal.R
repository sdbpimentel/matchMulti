#' Extract School-Level Covariates
#' 
#' Given a vector of variables of interest for students in a single school,
#' extracts a single value for the school
#' 
#' If the input is numeric, \code{agg} returns the mean; if the input is not
#' numeric, an error will be thrown unless all values are the same, in which
#' case the single unique value will be returned.
#' 
#' @param x a vector containing student-level observations for a school. If it
#' is a factor it must contain only a single level.
#' @return A single value of the same type as the input vector.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export agg
agg <-
function(x){
	if(length(unique(x)) == 1) return(x[1])
	return(mean(x, na.rm = TRUE))	
}



#' Collect Matched Samples
#' 
#' After students and schools have both been matched separately, assembles the
#' matched student samples corresponding to the school match into a single
#' dataframe of student-level data.
#' 
#' 
#' @param student.matches a list of lists object produced by
#' \code{matchStudents}, with each element of the second list containing a
#' dataframe composed of a matched sample for a different treated-control
#' school pairing.
#' @param school.match a dataframe, produced by \code{matchSchools}, with two
#' columns, one containing treated school IDs and the other containing matched
#' control school IDs.
#' @param school.id the name of the column storing the unique school identifier
#' (in the dataframes stored in \code{student.matches})
#' @param treatment the name of the column storing the binary treatment status
#' indicator (in the dataframes stored in \code{student.matches})
#' @return a dataframe containing the full set of matched samples for the
#' multilevel match.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export assembleMatch
assembleMatch <-
function(student.matches, school.match, school.id, treatment){
	#school.match <- school.match.list$matches[,1,drop = FALSE]
	#final.match <- student.matches[[1]][[1]]
	#final.match$pair.id <- rep(NA, nrow(final.match))
	#final.match <- final.match[-c(1:nrow(final.match)),]
	final.match <- NULL
	for(i in 1:nrow(school.match)){
		bind.obj <- student.matches[[school.match[i,1]]][[school.match[i,2]]]
		bind.obj$pair.id <- i
		if (is.null(final.match)) {
			final.match <- bind.obj
		} else {
			final.match <- rbind(final.match, bind.obj)			
		}
	}	
	return(final.match)
}

#' @export
#'
#' @rdname pairmatchelastic
elastic <- function (mdist, n = 0, val = 0) {
    st <- max(as.numeric(c(rownames(mdist), colnames(mdist))))
    h <- matrix(val, dim(mdist)[1], n)
    if(n > 0){
	    colnames(h) <- (st + 1):(st + n)
	}
    cbind(mdist, h)
}




fisher.balance <- function(varname, treatment, orig.data, match.data = NULL){
		
	orig.v <- orig.data[,which(colnames(orig.data) == varname)]
	
	tab.orig <- table(orig.v, orig.data[[treatment]])
	if(nrow(tab.orig) == 1) {
		t.orig <- list('p.value' = 1)
	} else {
		t.orig <- fisher.test(tab.orig)	
	}
	
	if(is.null(match.data)) return(c('Fisher Test Pvalue' = t.orig$p.value))

	match.v  <- match.data[,which(colnames(match.data) == varname)]	
	
	tab.match <- table(match.v, match.data[[treatment]])
	if(nrow(tab.match) == 1) {
		t.match <- list('p.value' = 1)
	} else {
		t.match <- fisher.test(tab.match)	
	}
	return(c('Fisher Test Pvalue Before' = t.orig$p.value, 'Fisher Test Pvalue After' = t.match$p.value))
}




#' Handle Missing Values
#' 
#' Preprocesses a dataframe of matching covariates so the Mahalanobis distance
#' can be calculated.
#' 
#' Preprocessing involves three main steps: (1) converting factors to matrices
#' of dummy variables (2) for any variable with NAs, adding an additional
#' binary variable indicating whether it is missing (3) imputing all NAs with
#' the column mean.  This follows the recommendations of Rosenbaum in section
#' 9.4 of the referenced text.
#' 
#' @param X a matrix or dataframe of covariates to be used for matching
#' @param verbose logical value indicating whether detailed output should be
#' provided.
#' @return a matrix containing the preprocessed data.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2010). \emph{Design of Observational
#' Studies}.  Springer-Verlag.
#' @keywords internal
#' @importFrom plyr laply
#' @export handleNA
handleNA <-
function(X, verbose = FALSE){
	if (is.data.frame(X)) {
        X.chars <- which(plyr::laply(X, class) == "character")
        if (verbose && length(X.chars) > 0) {
            print("character variables found in X, converting to factors")
            for (i in X.chars) {
                X[, i] <- factor(X[, i])
            }
        }
        X.factors <- which(plyr::laply(X, class) == "factor")
        for (i in which(plyr::laply(X, function(x) any(is.na(x))))) {
			if (verbose) print(paste('Missing values found in variable ', colnames(X)[i] ,'; imputing and adding missingness indicator'))        	
            if (i %in% X.factors) {
                X[, i] <- addNA(X[, i])
            }
            else {
                X[[paste(colnames(X)[i], "NA", sep = "")]] <- is.na(X[, 
                  i])
                X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
            }
        }
        for (i in rev(X.factors)) {
            X <- cbind(X[, -i], model.matrix(as.formula(paste("~", 
                colnames(X)[i], "-1")), data = X))
        }
    } else {
		if(is.null(colnames(X))) colnames(X) <- 1:ncol(X)
        for (i in c(1:ncol(X))) {
            if (any(is.na(X[, i]))) {
                X <- cbind(X, is.na(X[, i]))
                colnames(X)[ncol(X)] <- paste(colnames(X)[i], 
                  "NA", sep = "")
                X[which(is.na(X[, i])), i] <- mean(X[, i], na.rm = TRUE)
            }
        }
    }
	X
}




#' Check if a variable is binary
#' 
#' Examines a vector that is not coded as a logical to see if it contains only
#' 0s and 1s.
#' 
#' 
#' @param x A vector.
#' @return a logical value, \code{TRUE} if the vector contains only 0s and 1s
#' and \code{FALSE} otherwise.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export is.binary
is.binary <-
function(x){
	if (max(x) != 1 || min(x) != 0) return(FALSE)
	if (length(unique(x)) == 2) return(TRUE)
	#for now leave things with imputed means continuous
	#if (length(unique(x)) == 3){
	#	other.val <- setdiff(unique(x),c(0,1))
		#check if this other value is an imputed mean (due to NAs) - if so round and treat as binary	
	#}
	return(FALSE)
}




#' Compute School Distance from a Student Match
#' 
#' Defines a distance between two schools whose students have been matched
#' based on the size of the resulting matched sample and on the student-level
#' covariate balance.
#' 
#' The distance is computed by (1) subtracting the harmonic mean of the treated
#' and control counts in the matched sample from \code{largeval} (2) adding
#' \code{largeval} for each covariate among \code{studentvars} that has an
#' absolute standardized difference exceeding 0.2.  This encourages the school
#' match to choose larger schools with better balance.
#' 
#' @param matchFrame dataframe containing all matched students.
#' @param treatFrame dataframe containing all students from the treated school.
#' @param ctrlFrame dataframe containing all students from the control school.
#' @param student.vars names of variables on which to evaluate balance in the
#' matched sample.  Must be present in the column names of each of
#' \code{matchFrame}, \code{treatFrame} and \code{ctrlFrame}.
#' @param treatment name of the treatment variable. Must be present in the
#' column names of each of \code{matchFrame}, \code{treatFrame} and
#' \code{ctrlFrame}.
#' @param largeval a large penalty value to be added to the distance for each
#' student-level imbalance.
#' @return a numeric distance.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export match2distance
match2distance <-
function(matchFrame, treatFrame, ctrlFrame, student.vars, treatment, largeval){
	treat.v <- matchFrame[[treatment]]
	#compute harmonic mean of sample sizes and invert it
	harm.mean <- 2/(1/sum(treat.v) + 1/sum(1-treat.v) )
	#check balance and count imbalances
	if (!is.null(student.vars)) {
		sdiff.student <- sapply(student.vars, function(x) sdiff(x, treatment = treatment, orig.data= rbind(treatFrame, ctrlFrame), match.data = matchFrame)[6])
		return(largeval - harm.mean + largeval*sum(abs(sdiff.student) > 0.2))
	} else {
		return(largeval - harm.mean)		
	}
}



#' Match Schools on Student-based Distance
#' 
#' Takes in a school distance matrix created using information from the
#' first-stage student match and matches schools optimally, potentially
#' 
#' The \code{school.fb} argument encodes a refined covariate balance
#' constraint: the matching algorithm optimally balances the interaction of the
#' variables in the first list element, then attempts to further balance the
#' interaction in the second element, and so on.  As such variables should be
#' added in order of priority for balance.
#' 
#' @param dmat a distance matrix for schools, with a row for each treated
#' school and a column for each control school.
#' @param students a dataframe containing student and school covariates, with a
#' different row for each student.
#' @param treatment the column name of the binary treatment status indicator in
#' the \code{students} dataframe.
#' @param school.id the column name of the unique school ID in the
#' \code{students} dataframe.
#' @param school.fb an optional list of character vectors, each containing a
#' subset of the column names of \code{students}.  Each element of the list
#' should contain all the names in previous elements (producing a nested
#' structure).
#' @param penalty a numeric value, treated as the cost to the objective
#' function of excluding a treated school.  If it is set lower, more schools
#' will be excluded.
#' @param verbose a logical value indicating whether detailed output should be
#' printed.
#' @param tol a numeric tolerance value for comparing distances.  It may need
#' to be raised above the default when matching with many levels of refined
#' balance.
#' @return a dataframe with two columns, one containing treated school IDs and
#' the other containing matched control school IDs.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @importFrom plyr ddply
#' @importFrom rcbsubset rcbsubset
#' @keywords internal
#' @export matchSchools
matchSchools <-
function(dmat, students, treatment, school.id,
 school.fb, penalty, verbose, tol){
 	#DIST_TOLERANCE <- 1e-3
	
	balance.covs <- unique(c(unlist(school.fb)))
	school.df <- students2schools(students, c(balance.covs, treatment), school.id)
	
	#reorder school.df to ensure schools are in same order as implied by dimnames on dmat
	reord <- match(unlist(dimnames(dmat)), as.character(school.df[[school.id]]))
	if (any(is.na(reord))) stop('matchSchools found schools not handled by matchStudents')
	school.df <- school.df[reord,,drop = FALSE]
	if(is.null(school.fb)){
		match.out <- rcbsubset::rcbsubset(dmat,
		                       exclude.penalty = penalty, tol = tol) #	DIST_TOLERANCE)
	} else {
		match.out <- rcbsubset::rcbsubset(dmat, fb.list = school.fb, 
		                       treated.info = school.df[school.df[[treatment]] == 1,
		                                                ,drop = FALSE], 
		                       control.info = school.df[school.df[[treatment]] == 0,
		                                                ,drop = FALSE], 
		                       exclude.penalty = penalty, tol = tol)# tol =  	DIST_TOLERANCE)
	}
	out.frame <- cbind('TreatID' = school.df[[school.id]][which(school.df[[treatment]] == 1)][as.numeric(rownames(match.out$matches))], 'CtrlID' = school.df[[school.id]][which(school.df[[treatment]] == 0)][match.out$matches])
	
	return(out.frame)
}

#NB: when empty matches are formed between students, an NA goes into the distance
#matrix.  This works fine downstream because rcbsubset trims out NAs from
#distance matrix as forbidden matches because they are not finite.


#' Compute Student Matches for all Pairs of Schools
#' 
#' Iterates over all possible treated-control school pairs, optionally computes
#' and stores an optimal student match for each one, and generates a distance
#' matrix for schools based on the quality of each student match.
#' 
#' The \code{penalty.qtile} and \code{min.keep.pctg} control the rate at which
#' students are trimmed from the match.  If the quantile is high enough no
#' students should be excluded in any match; if the quantile is very low the
#' \code{min.keep.pctg} can still ensure a minimal sample size in each match.
#' 
#' @param students a dataframe containing student covariates, with a different
#' row for each student.
#' @param treatment the column name of the binary treatment status indicator in
#' the \code{students} dataframe.
#' @param school.id the column name of the unique school ID in the
#' \code{students} dataframe.
#' @param match.students logical value.  If \code{TRUE}, students are matched
#' within school pairs and some students will be excluded.  If \code{FALSE},
#' all students will be retained in the matched sample for each school pair.
#' @param student.vars column names of variables in \code{students} on which to
#' match students and assess balance of student matches in evaluating match
#' quality.
#' @param school.caliper matrix with one row for each treated school and one
#' column for each control school, containing zeroes for pairings allowed by
#' the caliper and \code{Inf} values for forbidden pairings.  When \code{NULL}
#' no caliper is imposed.
#' @param verbose a logical value indicating whether detailed output should be
#' printed.
#' @param penalty.qtile a numeric value between 0 and 1 specifying a quantile
#' of the distribution of all student-student matching distances.  The
#' algorithm will prefer to exclude treated students rather than form pairs
#' with distances exceeding this quantile.
#' @param min.keep.pctg a minimum percentage of students in the smaller school
#' in a pair which must be retained, even when treated students are excluded.
#' @return A list with two elements: \item{student.matches}{ a list with one
#' element for each treated school.  Each element is a list with one element
#' for each control school, and each element of these secondary lists is a
#' dataframe containing the matched sample for the corresponding
#' treated-control pairing. } \item{schools.matrix}{ a matrix with one row for
#' each treated school and one column for each control school, giving matching
#' distances based on the student match. }
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export matchStudents
matchStudents <-
function(students, treatment, school.id, match.students, student.vars, school.caliper = NULL, verbose, penalty.qtile, min.keep.pctg){	
	if(match.students){
		X <- students[student.vars]
		X <- handleNA(X, verbose = verbose)
		varying <- apply(X,2, function(x) length(unique(x)) > 1)
		if(!all(varying) && verbose){
			print('Constant-value variables found, they will not be used to calculate Mahalanobis distance.')
		} 
		X <- X[,which(varying),drop = FALSE]
		#TODO: examine performance of this step, consider optimization
		maha.mat <- smahal(students[[treatment]],X)
		#TODO: decide how to handle calipers
		student.penalty <- quantile(maha.mat, penalty.qtile)
	}
	
	#TODO: add documentation to show treatment should be 0/1
	t.schools <- unique(students[[school.id]][students[[treatment]]==1])
	c.schools <- unique(students[[school.id]][students[[treatment]]==0]) 
	nt.schools <- length(t.schools)
	nc.schools <- length(c.schools)
	schoolmatch.mat <- matrix(nrow = nt.schools, ncol = nc.schools)
	matches.list <- list()
	for (i in 1:nt.schools) {
		matches.list[[t.schools[i]]] <- list()
		schli <- students[students[[school.id]] == t.schools[i],,drop = FALSE]
 		for (j in 1:nc.schools) {
			schlj <- students[students[[school.id]] == c.schools[j],,drop = FALSE] 
			if (match.students) {
				#grab correct slice of maha distance
				treat.idx <- c(students[[school.id]] == t.schools[i])[students[[treatment]] == 1]
				ctrl.idx <- c(students[[school.id]] == c.schools[j])[students[[treatment]] == 0]
				maha.slice <- maha.mat[which(treat.idx), which(ctrl.idx),drop = FALSE]
				min.keep <- floor((1-min.keep.pctg)*min(dim(maha.slice)))
				match.out <- pairmatchelastic(maha.slice, n = min.keep, val = student.penalty)
				if(is.null(match.out)) return(NULL)
				matches.list[[t.schools[i]]][[c.schools[j]]] <- rbind(schli[as.numeric(rownames(match.out)),,drop = FALSE], schlj[match.out,,drop = FALSE])
			} else {
				matches.list[[t.schools[i]]][[c.schools[j]]] <- rbind(schli,schlj)
			}
			schoolmatch.mat[i,j] <- match2distance(matches.list[[t.schools[i]]][[c.schools[j]]], schli, schlj, student.vars, treatment, largeval = nrow(students)/2)
		}
	}
	if(!is.null(school.caliper)) schoolmatch.mat <- schoolmatch.mat + school.caliper

	rownames(schoolmatch.mat) <- t.schools
	colnames(schoolmatch.mat) <- c.schools

	return(list('student.matches' = matches.list, 'schools.matrix' = schoolmatch.mat))
}





#' Optimal Subset Matching without Balance Constraints
#'
#' Conducts optimal subset matching as described in the reference.
#'
#' \code{pairmatchelastic} is the main function, which conducts an entire match.
#' \code{elastic} is a helper function which augments the original distance
#' matrix as described in the reference. 
#' 
#' The original versions of these functions were written by Paul Rosenbaum and
#' distributed in the supplemental material to the paper: "Optimal Matching of
#' an Optimally Chosen Subset in Observational Studies," Paul R. Rosenbaum,
#' Journal of Computational and Graphical Statistics, Vol. 21, Iss. 1, 2012.
#'
#' @param mdist distance matrix with rows corresponding to treated units and
#'   columns corresponding to controls.
#' @param n maximum number of treated units that can be excluded.
#' @param val cost of excluding a treated unit (i.e. we prefer to exclude a
#'   treated unit if it increases the total matched distance by more than
#'   \code{val}).
#' @return \code{elastic} returns an augmented version of the input matrix
#'   \code{mdist}.  \code{pairmatchelastic} returns a matrix of 1 column whose
#'   values are the column numbers of matched controls and whose rownames are
#'   the row numbers of matched treated units.
#' @author Paul R. Rosenbaum (original forms), modifications by Luke Keele and
#'   Sam Pimentel
#'   
#' @references Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally
#'   Chosen Subset in Observational Studies."  Journal of Computational and
#'   Graphical Statistics, 21.1, 57-71.
#'   
#' @export pairmatchelastic
pairmatchelastic <- function (mdist, n = 0, val = 0) {
    ro <- dim(mdist)[1]
    co <- dim(mdist)[2]
    k <- ro + co
    mdist <- elastic(mdist, n = n, val = val)
    m <- rcbsubset::rcbsubset(mdist)
	if(is.null(m)) return(NULL)
	mt <- m
	mt$matches <- m$matches[which(m$matches <= co),,drop=FALSE]
    mt$matches
}




#' Ensure Dataframes Share Same Set Columns
#' 
#' Takes in two dataframes.  For each column name that is in the second frame
#' but not in the first frame, a new column of zeroes is added to the first
#' frame.
#' 
#' 
#' @param df1 a dataframe.
#' @param df2 a dataframe.
#' @return a dataframe
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @export resolve.cols
resolve.cols <-
function(df1, df2){
	new.cols <- setdiff(colnames(df2), colnames(df1))
	for(newc in new.cols){
		df1[[newc]] <- rep(0, nrow(df1))
	}
	df1
}



#' Balance Measures
#'
#' Balance assessment for individual variables, before and after matching
#'
#' The \code{sdiff} function computes the standardized difference in means. The
#' other functions perform different kinds of balance tests: \code{t.balance}
#' does the 2-sample t-test, \code{fisher.balance} does Fisher's exact test for
#' binary variable, and \code{wilc.balance} does Wilcoxon's signed rank test.
#'
#' @aliases sdiff ttest.balance fisher.balance wilc.balance
#' @param varname name of the variable on which to test balance
#' @param treatment name of the binary indicator for treatment status
#' @param orig.data a data frame containing the data before matching
#' @param match.data an optional data frame containing the matched sample
#'
#' @return a labeled vector.  For \code{sdiff}, the vector has six elements if
#'   \code{match.data} is provided: treated and control means and standardized
#'   differences before and after matching.  If \code{match.data} is not
#'   provided, the vector has only the three elements corresponding to the
#'   pre-match case.
#'
#'   For the other functions, if \code{match.data} is provided, the vector
#'   contains p-values for the test before and after matching. Otherwise a
#'   single p-value is given for the pre-match data.
#' 
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#'
#'   Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#'   Springer-Verlag.
#'
#'   Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#'   Springer-Verlag.
#' @keywords internal
#' @export sdiff
sdiff <-
function(varname, treatment, orig.data, match.data = NULL) {
		
	orig.v <- orig.data[,which(colnames(orig.data) == varname)]
	m1 <- mean(orig.v[orig.data[[treatment]] == 1], na.rm = TRUE)
	m0 <- mean(orig.v[orig.data[[treatment]] == 0], na.rm = TRUE)
	
	s1 <- sqrt(var(orig.v[orig.data[[treatment]] == 1], na.rm= TRUE))
	if(length(orig.v[orig.data[[treatment]] == 1]) == 1) s1 <- 0
	s0 <- sqrt(var(orig.v[orig.data[[treatment]] == 0], na.rm = TRUE))
	if(length(orig.v[orig.data[[treatment]] == 0]) == 1) s0 <- 0
	sdiff.before <- (m1 -m0)/sqrt(0.5*(s1^2 + s0^2))
	if(s1 == 0 && s0 == 0) sdiff.before <- ifelse(m1 == m0, 0, Inf)
	if(is.null(match.data)) return(c('Treated Mean' = m1, 'Control Mean' = m0,'SDiff' = sdiff.before))
		
	#if match data is also provided
	match.v  <- match.data[,which(colnames(match.data) == varname)]
	m1.m <- mean(match.v[match.data[[treatment]] == 1], na.rm= TRUE)
	m0.m <- mean(match.v[match.data[[treatment]] == 0], na.rm = TRUE)	

	sdiff.after <- (m1.m - m0.m)/sqrt(0.5*(s1^2 + s0^2))
	if(s1 == 0 && s0 == 0) sdiff.after <- ifelse(m1.m == m0.m, 0, Inf)

	return(c('Treated Mean Before' = m1, 'Control Mean Before' = m0,'SDiff Before' = sdiff.before, 'Treated Mean After' = m1.m, 'Control Mean After' = m0.m, 'SDiff After' = sdiff.after))
}

## The original version of the following function was written by Paul Rosenbaum and distributed in the supplemental material to the paper: "Optimal Matching of an Optimally Chosen Subset in Observational Studies," Paul R. Rosenbaum, Journal of Computational and Graphical Statistics, Vol. 21, Iss. 1, 2012.


#' Robust Mahalanobis Distance
#' 
#' Computes robust Mahalanobis distance between treated and control units.
#' 
#' For an explanation of the robust Mahalanobis distance, see section 8.3 of
#' the first reference.  This function was written by Paul Rosenbaum and
#' distributed in the supplemental material to the second reference.
#' 
#' @param z vector of treatment indicators (1 for treated, 0 for controls).
#' @param X matrix of numeric variables to be used for computing the
#' Mahalanobis distance.  Row count must match length of \code{z}.
#' @return a matrix of robust Mahalanobis distances, with a row for each
#' treated unit and a column for each control.
#' @author Paul R. Rosenbaum.
#' @references Rosenbaum, Paul R. (2010). \emph{Design of Observational
#' Studies}.  Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2012) "Optimal Matching of an Optimally Chosen Subset in
#' Observational Studies."  Journal of Computational and Graphical Statistics,
#' 21.1, 57-71.
#' @keywords internal
#' @export smahal
smahal <-
function(z,X){
    X<-as.matrix(X)
    n<-dim(X)[1]
    rownames(X)<-1:n
    k<-dim(X)[2]
    m<-sum(z)
    for (j in 1:k) X[,j]<-rank(X[,j])
    cv<-cov(X)
    vuntied<-var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
    out<-matrix(NA,m,n-m)
    Xc<-X[z==0,]
    Xt<-X[z==1,]
    rownames(out)<-rownames(X)[z==1]
    colnames(out)<-rownames(X)[z==0]
    icov<-ginv(cv)
    for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
    out
}




#' Aggregate Student Data into School Data
#' 
#' Takes a dataframe of student-level covariates and aggregates selected
#' columns into a dataframe of school covariates.
#' 
#' Aggregation is either done by taking averages or by selecting the unique
#' factor value when a school has only one value for a factor.  As a result,
#' \code{school.covs} should only include variables that are numeric or do not
#' vary within schools.
#' 
#' @param students a dataframe of students.
#' @param school.cov a character vector of column names in \code{students} that
#' should be aggregated by school.
#' @param school.id the name of the column in \code{students} containing the
#' unique school identifier.
#' @return a dataframe of aggregated data, with one row for each school and
#' columns in \code{school.covs} and \code{school.id}.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @keywords internal
#' @importFrom plyr ddply colwise
#' @export students2schools
students2schools <-
function(students, school.cov, school.id){
	school.df <- students[c(school.id, school.cov)]
	school.df <- plyr::ddply(school.df, school.id, plyr::colwise(agg))
	school.df	
}





#' Outcome analysis.
#' 
#' Calculates confidence interval via grid search.
#' 
#' 
#' @param beta Confidence interval value
#' @param obj a multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @param alpha Level of test for confidence interval, default is .05 for 95\%
#' CI.
#' @param alternative Direction of test.
#' @return The endpoint of an estimated confidence interval.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export ci_func
ci_func <- function(beta, obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, alpha, alternative="less"){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     match.data$adj.y <- match.data[[out.name]] - z*beta
     strata.ranks <- tapply(match.data$adj.y, match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     #One-sided p-value for Test of the Sharp Null
     if(alternative=="less"){
     	 pval <- pnorm(T.c/sqrt(var.T.c))
     } else {
     	pval <- 1 - pnorm(T.c/sqrt(var.T.c))
     }
     pval - alpha
}



#' Outcome analysis.
#' 
#' Calculates Hodges-Lehmann point estimate for treatment effect.
#' 
#' 
#' @param beta Point estimate value
#' @param obj A multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @return A point estimate for constant-additive treatment effect.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export pe_func
pe_func <- function(beta, obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     z <- match.data[[treat.name]]
     match.data$adj.y <- match.data[[out.name]] - z*beta
     strata.ranks <- tapply(match.data$adj.y, match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     E.T <- T.c/var.T.c
     E.T - 0
     }



#' Outcome analysis.
#' 
#' Calcualtes p-values for test of sharp null for treatment effect.
#' 
#' 
#' @param obj A multiMatch object
#' @param out.name Name of outcome covariate
#' @param schl_id_name Name of school (group) identifier
#' @param treat.name Name of treatment indicator
#' @param wt Logical flag for whether p-value should weight strata by size.
#' @return A p-value for constant-additive treatment effect.
#' @author Luke Keele, Penn State University, \email{ljk20@@psu.edu}
#' 
#' Sam Pimentel, University of Pennsylvania, \email{spi@@wharton.upenn.edu}
#' @references Rosenbaum, Paul R. (2002). \emph{Observational Studies}.
#' Springer-Verlag.
#' 
#' Rosenbaum, Paul R. (2010). \emph{Design of Observational Studies}.
#' Springer-Verlag.
#' @keywords internal
#' @export pval_func
pval_func <- function(obj, out.name = NULL, schl_id_name = NULL, treat.name = NULL, wt=TRUE){
	 match.data <- obj$matched
     #No of Strata
     n.s <- dim(unique(match.data[paste(schl_id_name)]))[1]
     Q.s <- matrix(NA, n.s, 1)

    #Rank Within Paired Strata
     strata.ranks <- tapply(match.data[[out.name]], match.data$pair.id, rank, simplify=TRUE)
     match.data$ranks <- matrix(unlist(strata.ranks), nrow(match.data), 1)
     z <- match.data[[treat.name]]
     sub.trt <- match.data[z==1,]
     sub.ctrl <- match.data[z==0,]
     n.s1 <- tapply(sub.trt$ranks, sub.trt$pair.id, length)
     n.s2 <- tapply(sub.ctrl$ranks, sub.ctrl$pair.id, length)
     d <- (sum(n.s1 + n.s2))
     q1 <- as.vector(tapply(sub.trt$ranks, sub.trt$pair.id, mean))
     q2 <- as.vector(tapply(sub.ctrl$ranks, sub.ctrl$pair.id, mean))

     #Constant Weights
     if(wt==FALSE){
     	     ws.1 <- 1
     Q.s.c <- (q1- q2)*ws.1
     T.c <- sum(Q.s.c)
     var.T.c <- sum(Q.s.c^2)
     #One-sided p-value for Test of the Sharp Null
     pval.c1 <- pnorm(T.c/sqrt(var.T.c))
     pval.c2 <- 1 - pnorm(T.c/sqrt(var.T.c))
     pval.c <- min(pval.c1, pval.c2)
     return(pval.c)	
     } else {
     d <- (sum(n.s1 + n.s2))
     ws.2 <- (n.s1 + n.s2) /d
     Q.s.p <- ((q1- q2))*ws.2
     T.p <- sum(Q.s.p)
     var.T.p <- sum(Q.s.p^2)
     #One-sided p-value for Test of the Sharp Null
     pval.p1 <- pnorm(T.p/sqrt(var.T.p))
     pval.p2 <- 	1 -  pnorm(T.p/sqrt(var.T.p))
     pval.p <- min(pval.p1, pval.p2)
     return(pval.p)
     }
 	
}


	.Random.seed <-
c(403L, 2L, -1639498799L, -754364541L, -1683490286L, 379350436L, 
1078499687L, 1589511869L, -1219400140L, 1029600666L, -1054625307L, 
-1532756609L, -429804282L, 1011830416L, 407488307L, -758999855L, 
-565162368L, 951015262L, 562803881L, 1605322587L, 421379434L, 
1469141004L, -1433597873L, 1321623813L, -68085316L, 1084279954L, 
2060973965L, 2062577607L, -2012048210L, 478958856L, -2033024245L, 
2000236009L, 1706228408L, -1293755034L, 1336923649L, 582598803L, 
2125527874L, 927839476L, -1437490953L, 1989155053L, -2031543420L, 
-188140118L, -1244812523L, -2004854353L, -866887242L, 2103017440L, 
-1064677597L, 938057793L, 1061976880L, 1083126094L, 1971594169L, 
1201165771L, -2110282118L, -573764740L, 1136788415L, -2062193515L, 
-1788546388L, -1049513854L, -1909653411L, 128169815L, -1554011586L, 
249765496L, 1682373979L, 1157557433L, 765232392L, -2046851146L, 
-1282986767L, 798518051L, -1036391054L, -2046359612L, -215374969L, 
1790472541L, 1039150932L, 1492841850L, 1313624261L, 1529311775L, 
-1725696602L, 293566256L, 630947283L, 2092141553L, 2004020512L, 
-1676088386L, 446817225L, -366267525L, 194187594L, -2095546964L, 
696780975L, 604803685L, -2058728932L, 1020224626L, -1908319507L, 
-13383897L, -307207346L, -1165857880L, 468209579L, -1123722039L, 
965029080L, 1030975878L, -1737340831L, 1539811827L, -838018974L, 
-1478766892L, 1623234519L, 1206699469L, 1920000356L, -1483823798L, 
-1293736779L, -245554161L, 496155926L, -980187200L, 236729219L, 
-717192799L, -181523056L, -661201298L, 128616473L, -225127253L, 
-286824294L, 2095544412L, 1554969951L, 419233589L, 509694668L, 
1867124770L, 2089564797L, -1457205385L, -1493172066L, 887242584L, 
-160929093L, -932383655L, -200970264L, -1157940714L, 1680555409L, 
-1330469309L, -424659630L, -1087366684L, -443459545L, 1879090301L, 
139790708L, -397851302L, 395809701L, -915651649L, 808284230L, 
1316650704L, -41056013L, 1106894353L, 977610816L, -218833890L, 
287095529L, 474882075L, -547624918L, -1539015860L, 106535183L, 
-2063250363L, -880734724L, 781958226L, -661283379L, -1072283769L, 
-719979922L, -2019001784L, -1068528565L, 102952361L, 828637432L, 
-1211997146L, 1746527553L, -489008429L, -801875326L, -184319820L, 
-1357549513L, -84384979L, -925018044L, -1435441814L, 194859477L, 
1395001199L, 1554169462L, -2046660192L, -678234525L, -633532671L, 
-1864827280L, 428546574L, -1019833735L, -1407533557L, -1033105606L, 
1078490684L, -1952621185L, -662901035L, 1473937132L, 1339576386L, 
2145105821L, -1993193961L, -1167115394L, 778000696L, -1222969829L, 
530798457L, 1527763400L, 73960950L, -716148687L, 1295694435L, 
-2088475086L, -415272060L, -1477946937L, -647530339L, 270872596L, 
1485522618L, 1252579077L, -1018039073L, 1721955686L, -706707216L, 
-1034845165L, 774998449L, -1676786592L, 1761710590L, 1991259529L, 
10089147L, -1307418998L, 1757452652L, -1144463633L, -424810203L, 
1826686940L, -613685838L, -1760136019L, -1295592345L, -598054002L, 
67349608L, -208369557L, 254906889L, 1082383512L, -1424896058L, 
397517601L, -1193152589L, -845330654L, 122612540L, -534358158L, 
2068021040L, 8320548L, -375641032L, 256272602L, 28016016L, 1716247292L, 
1635820436L, -322871486L, 606012928L, 63172156L, -382725552L, 
-1737023838L, 310528840L, 683926028L, -652781348L, 1451366994L, 
-1006569344L, 263177044L, 426429928L, 1324797690L, 1829915472L, 
473205804L, -1990182540L, 135665298L, 555179616L, 2021881484L, 
1470830432L, 1360339682L, -508768088L, -115024692L, -1410756100L, 
861071282L, -690730512L, 1699132868L, 1884713816L, 677367674L, 
1310254512L, -1470587460L, 42445908L, -287387710L, 143791904L, 
581804988L, 2065948240L, 2088999394L, -454777496L, 1846793548L, 
-1011598084L, 1735500530L, 1864289024L, -2000130444L, -625869816L, 
-1338373702L, 1747594192L, -938289684L, 554270612L, 1893715666L, 
592480L, 440948044L, -159477472L, 1225221762L, -746093080L, -1089905684L, 
2008455740L, 2106795634L, 1215978352L, -817088220L, 1752410680L, 
1672522138L, -484037168L, 373558716L, -317610860L, 343330946L, 
202732992L, 931170108L, -1766057840L, 2003831586L, -141458936L, 
-1071546100L, -323081828L, -876897774L, -264151168L, -2059242732L, 
-2028627672L, 309927994L, 351228880L, 1366722156L, 658141748L, 
-1132702062L, 1986661728L, -308141748L, -1318226976L, -132857438L, 
-251010648L, 1113161356L, 1191278012L, 1608262258L, 1822501872L, 
-1932223932L, -1360388008L, 242746L, -569992080L, -213449668L, 
1878318356L, -1675959614L, -356149536L, -1660265348L, -938671408L, 
-506601054L, 1329566760L, -243755252L, -397507140L, -868636878L, 
1152552640L, 156231732L, -1843482936L, 1454598330L, -1869331824L, 
1485416172L, -347385260L, -363917678L, 101668128L, 1589124556L, 
266596064L, -1356022590L, 1472781352L, -286219348L, -2131687748L, 
-837466126L, 724447792L, -453813340L, 1533060536L, 147785562L, 
-1108862832L, -812324996L, -2126256876L, -1189002942L, -253035264L, 
550822332L, 2078219600L, 359270690L, -1328647608L, -1936499700L, 
-228560292L, 1837266386L, 278380032L, -971468460L, -114064152L, 
-401843974L, 1234957008L, 87514924L, 635485428L, -1429523054L, 
444885344L, -991123444L, 910878944L, 812869218L, -2023886680L, 
-1132128180L, -722774148L, -52028110L, 1594330736L, -1060075708L, 
793332696L, -1542647942L, 1015727152L, -44687428L, -1386456492L, 
491475010L, -1216145760L, -1624749380L, 142051280L, -1163925150L, 
123259624L, 1339692236L, 290753020L, -698093070L, 390664320L, 
1101753844L, 1973127816L, -1332558022L, -971671216L, 1819883628L, 
-2073311852L, -73585582L, 172298720L, 1653801164L, 339421856L, 
1169305986L, -1423813400L, 1222720236L, 871706940L, 140405874L, 
472096368L, -1365912028L, 1934565432L, -1161054950L, -327487280L, 
526052668L, -947761260L, -1058719870L, -693932096L, -1474810308L, 
1955160464L, 1572275234L, 1619271944L, 760702092L, -1238117604L, 
517574802L, 1269383552L, 1466416148L, 2046326568L, -1729675718L, 
-1341103664L, -1801019028L, 379826100L, -349686254L, 1293063136L, 
1675564492L, 810403040L, 796543010L, -914911448L, -1099987444L, 
-418165956L, -1538483214L, 781947632L, 1526263916L, 480774333L, 
755153287L, -1458819728L, -1667917810L, -1116940997L, 327432733L, 
-2125979958L, 411392488L, 70520529L, 706071235L, -855277116L, 
-1743444366L, 871225159L, 2049951505L, -1631995242L, -1543348236L, 
-1915655803L, -1955863345L, -1289495096L, 1865856102L, -465964397L, 
1219836405L, -1615719470L, -1576941536L, -1189273079L, 500939579L, 
-593523060L, 779756314L, -1464425073L, -1290863911L, 844682606L, 
-2054620932L, -720616851L, -86761961L, -454419488L, 1605199806L, 
-1336797813L, 267388429L, 1335614010L, 1589228152L, 346298913L, 
510244051L, -1452428300L, 541332610L, -407421673L, -811045919L, 
-1768145626L, -768803804L, -410224939L, 1266954879L, 965010200L, 
-372573642L, 1445043395L, 1221538053L, 1169303714L, 88178320L, 
-1620215239L, 1231651563L, -1299741732L, -1789860342L, -278701825L, 
-2008449719L, 436186846L, 1504580748L, -1828543395L, 139385447L, 
66323344L, 2007132078L, -1968167333L, -1010058691L, -1985344342L, 
-1029915320L, -715685519L, 389048035L, -549179932L, -1049190126L, 
-63136537L, -902059343L, 1998860534L, -434078956L, 1328796709L, 
-1670958801L, 1633568488L, -29264954L, -1465014669L, 140282069L, 
-246115662L, 1916836480L, 2089130601L, 1941491419L, 1440262700L, 
-1607507398L, 532032751L, -2085208199L, -1684618226L, 1111947740L, 
-1374102963L, 2022402487L, 1487374144L, -279708130L, -1452439317L, 
-581961619L, 1438126298L, 1596600472L, 1731331841L, -409563981L, 
2050467796L, -654900510L, 849601399L, 1919815361L, 1919927878L, 
-2082522748L, 1397458741L, 2129392927L, -795419528L, -2113369770L, 
-1494915101L, -527699547L, -1724898238L, -592578000L, 300044889L, 
1373684683L, -1413454788L, 375122538L, 50636959L, 264765737L, 
1460868542L, -1834261332L, 1014407677L, 1677961543L, 773481008L, 
868213582L, -1986271365L, -18853667L, 923179274L, -715359192L, 
1903702161L, -155499389L, -692296316L, 2088666034L, -1136849657L, 
-50243759L, 759174870L, -541159500L, -120266811L, 420978447L, 
558308744L, -1395050330L, 97658963L, 833573301L, -1914169582L, 
1603796192L, 1813041481L, 1359440763L, 1023088204L, 511092698L, 
-1091521073L, -1790966887L, 1594594222L, -44672836L, 807203117L, 
562889431L, 862333728L, 1361529982L, 555758155L, 1018986573L, 
1429380730L, -1368963400L, -882868875L)
