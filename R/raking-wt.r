


#' Survey Raking in R
#'
#' Performs raking of survey weights, using a technique very similar to the SAS IHB Macro
#'
#' @param formula a formula object, with the input/start weights on the left of the ~ operator and the raking variables on the right, with each raking variable separated by the + operator. Two-way variable raking can be done by using the * operator.
#' @param data data frame of respondents in which to interpret the variables named in the formula.
#' @param targets a single data frame or a list of data frames containing the target/calibration totals/figures. The raking variables listed in the formula must match the column names of those data frames. The target data frames must be listed in the same order as their corresponding raking variables in the formula.
#' @param namesTotals vector, where each element is the name of the column that contains the target/calibration totals in the corresponding data frame listed in targets.
#' @param maxiter (default 50) maximum number of iterations to attempt for convergence.
#' @param epsilon (default 1) convergence is declared if the maximum of the differences between marginal total of the raked weights and the target/calibration totals is less than or equal to epsilon. If epsilon<1 it is taken to be a fraction of the total of the sinput/start weights in the formula.
#' @param freq Logical. If true (default), frequency tables are displayed (this is useful to see if cells should be collapsed).
#' @param freqWarn a warning message will be displayed if one or more cells contained fewer than this number (default 15) of observations.
#' @param verbose controls the amount of output displayed on the screen.
#' @param trim logical, if true large weights will be trimmed. This is done using a technique very similar to the IGCV method of Izrael et al. (2009)
#' @param trim.method "median" (default) or "mean"; see trim.max
#' @param trim.max (default 6) each weight w is trimmed if w > trim.max * trim.method(weights)
#' @param trim.group if NULL (default), the overall mean/median of the weights is used when trimming; if trim.group is the name of a column of data, then the weights will be trimmed based on which group of trim.group they belong to.
#' @examples
#' # Packages
#' library(PracTools)
#' library(sampling)
#' library(dplyr)
#'
#' # Population data
#' data(MDarea.pop)
#'
#' # population counts for two strata variables
#' table(MDarea.pop$BLKGROUP, MDarea.pop$Gender)
#'
#' # take a randomly stratified sample based on BLKGROUP and Gender
#' set.seed(20230128)
#' sample = sampling::strata(
#'   data=MDarea.pop,
#'   stratanames=c("BLKGROUP","Gender"),
#'   size=sample(30:100, size=12, replace=TRUE),
#'   method="srswor"
#' )
#'
#' #### example 1: rake only on Gender
#' # population totals
#' targets_Gender = data.frame(
#'   Gender=1:2,
#'   popcount=as.numeric(table(MDarea.pop$Gender)))
#'
#' # initial weights
#' sample_ex1 = sample %>%
#'   mutate(wgt.init = 1)
#'
#' # raking
#' rake_ex1 = raking(
#'   wgt.init ~ Gender,
#'   data=sample_ex1,
#'   targets=list(targets_Gender),
#'   namesTotals="popcount"
#' )
#' sample_ex1$wgt.rake = rake_ex1$weights

#' # check that sum of weights match population targets
#' sample_ex1 %>%
#'   group_by(Gender) %>%
#'   summarise(wgttotals = sum(wgt.rake))
#'
#'
#'
#' #### example 2: rake on Gender and Geographic area using marginal totals
#' # population totals
#' targets_Gender = data.frame(
#'   Gender=1:2,
#'   popcount=as.numeric(table(MDarea.pop$Gender)))
#' targets_BLKGRP = data.frame(
#'   BLKGROUP=1:6,
#'   popcount=as.numeric(table(MDarea.pop$BLKGROUP)))
#'
#' # initial weights
#' sample_ex2 = sample %>%
#'   mutate(wgt.init = 1)
#'
#' # raking
#' rake_ex2 = raking(
#'   wgt.init ~ Gender + BLKGROUP,
#'   data=sample_ex2,
#'   targets=list(targets_Gender, targets_BLKGRP),
#'   namesTotals="popcount"
#' )
#' sample_ex2$wgt.rake = rake_ex2$weights
#'
#' # check that sum of weights match population targets
#' sample_ex2 %>%
#'   group_by(Gender) %>%
#'   summarise(wgttotals = sum(wgt.rake))
#' sample_ex2 %>%
#'   group_by(BLKGROUP) %>%
#'   summarise(wgttotals = sum(wgt.rake))
#'
#'
#'
#' #### example 3: rake on all Gender-GeographicArea bins
#' # (equivalent to post-stratification)
#' # population totals
#' targets_GenderBLK = as.data.frame(table(MDarea.pop$BLKGROUP, MDarea.pop$Gender)) %>%
#'   rename(BLKGROUP = Var1,
#'          Gender   = Var2,
#'          popcount = Freq)
#'
#' # initial weights
#' sample_ex3 = sample %>%
#'   mutate(wgt.init = 1)
#'
#' # raking
#' rake_ex3 = raking(
#'   wgt.init ~ Gender*BLKGROUP,
#'   data=sample_ex3,
#'   targets=list(targets_GenderBLK),
#'   namesTotals="popcount"
#' )
#' sample_ex3$wgt.rake = rake_ex3$weights
#'
#' # check that sum of weights match population targets
#' sample_ex3 %>%
#'   group_by(Gender,BLKGROUP) %>%
#'   summarise(wgttotals = sum(wgt.rake))
#' @references
#' * Izrael, Hoaglin & Battaglia (2000), ["A SAS Macro for Balancing a Weighted Sample"](http://www2.sas.com/proceedings/sugi25/25/st/25p258.pdf), Proceedings of the 25th Annual SAS Users Group International Conference, Paper 258
#' * Izrael, Battaglia & Frankel (2009), ["Extreme Survey Weight Adjustment as a Component of Sample Balancing (a.k.a. Raking)"](https://support.sas.com/resources/papers/proceedings09/247-2009.pdf), SAS Global Forum, Paper 247-2009
#' @export
raking = function(
    formula, data, targets, namesTotals, maxiter=50, epsilon=1, freq=TRUE, freqWarn=15,
		verbose=c("some", "none", "full"), trim=FALSE, trim.method=c("median", "mean"), trim.max=6, trim.group=NULL)
{

  # call matching & input verification
  call <- match.call()
  verbose <- match.arg(verbose)
  trim.method <- match.arg(trim.method)
  if (!plyr::is.formula(formula)) stop("'formula' must be a formula")
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  if (!is.list(targets)) stop("'targets' must be a single data frame or a list of data frames")
  if (is.data.frame(targets)) targets <- list(targets)
  for (i in 1:length(targets))
      {if (!is.data.frame(targets[[i]])) stop("each element of 'targets' must be a data.frame")}

  if (!is.numeric(maxiter) || length(maxiter)!=1 || maxiter<=0 || maxiter%%1!=0)
    stop("'maxiter' must be a positive integer")
  if (!is.numeric(epsilon) || length(epsilon)!=1 || epsilon<=0 || (epsilon>=1 && epsilon %%1!=0))
    stop("'epsilon' must be a proportion (i.e., between 0 and 1) or a positive integer")
  if (!is.logical(freq)) stop("'freq' must be TRUE or FALSE")
  if (!is.numeric(freqWarn) || length(freqWarn)!=1 || freqWarn<=0 || freqWarn%%1!=0)
    stop("'freqWarn' must be a positive integer")

  if (!is.logical(trim)) stop("'trim' must be TRUE or FALSE")
  if (!is.numeric(trim.max) || length(trim.max)!=1 || trim.max<=0)
    stop("'trim.max' must be a positive number")
  if (trim && !is.null(trim.group) && is.na(match(trim.group, names(data)))) stop("'", trim.group, "' is not valid column of 'data'")

  # extracting the names of the target datasets
  namesTargets <- paste(deparse(call[[4]]), collapse="")
  #namesTargets <- "list(targetsFRA_w1_sexAge.dat, targetsFRA_w1_region.dat, targetsFRA_w1_educ.dat)"
  namesTargets <- gsub("[[:space:]]+", "", namesTargets)
  namesTargets <- unlist(strsplit(namesTargets, "\\(|,|\\)"))
  namesTargets <- namesTargets[-1]

  # extracting the variable names from the formula
  formula <- gsub("[[:space:]]+", "", formula)
  tmp <- strsplit(as.character(formula), "~|\\+")
  nameWts <- tmp[[2]]
  if (is.na(match(nameWts, names(data)))) stop("'", nameWts, "' is not a valid column of 'data'")
  levels <- strsplit(tmp[[3]], "\\*")
  if(length(levels) != length(targets)) stop("The number of raking variables listed in the formula doesn't match the number of listed data frames in 'targets'")

  # initializations
  nbTargets <- length(targets)
  if (length(namesTotals)==1) namesTotals <- rep(namesTotals, nbTargets)
  if (epsilon < 1) epsilon <- epsilon * sum(data[ , nameWts])
  wtsRake <- NULL

  # combining the values of levels into a single column
  for (i in 1:nbTargets)
      {
      # checking that levels[[i]] and namesTotals[i] are valid columns of targets[[i]]
      tmp <- match(c(levels[[i]], namesTotals[i]), names(targets[[i]]))
      if (any(is.na(tmp))) stop("One or more of the following are not valid column(s) of '",
      			namesTargets[i], "': ", paste(c(levels[[i]], namesTotals[i]), collapse=", "))

      # combining levels[[i]] into a single column of targets[[i]]
  	tmp <- as.character(interaction(targets[[i]][ , levels[[i]]], drop=T))
      targets[[i]] <- cbind(" "=tmp, targets[[i]])

      # checking that levels[[i]] are valid columns of data
      tmp <- match(levels[[i]], names(data))
      if (any(is.na(tmp))) stop("One or more of the following are not valid column(s) of 'data': ", paste(levels[[i]], collapse=", "))

      # combining levels[[i]] into a single column, doing so creating the wtsRake dataset (part 1)
      tmp <- as.character(interaction(data[ , levels[[i]]], drop=T))
      wtsRake <- cbind(wtsRake, tmp)
      }

  # putting the wtsRake data frame together (part 2)
  if (is.null(trim.group))
     {
     wtsRake <- data.frame(wtsRake, data[ , nameWts])
     colnames(wtsRake) <- c(paste("cell", 1:nbTargets, sep=""), "weights")
     } else {
   	      wtsRake <- data.frame(wtsRake, data[ , nameWts], data[ , trim.group])
     		  colnames(wtsRake) <- c(paste("cell", 1:nbTargets, sep=""), "weights", trim.group)
     	       }

  # printing frequency tables
  if (freq)
     {
     for (i in 1:nbTargets)
          {
          tmp <- as.data.frame(table(wtsRake[ , i]))
          names(tmp)[1] <- ""
          cat("Frequency table for ", paste(levels[[i]], collapse=" x "), ": \n", sep="")
  	    print(tmp)
  	    if (any(tmp$Freq < freqWarn)) message("One or more cells with frequency < ", freqWarn)
  	    cat("\n\n")
  	    }
     }

  # raking the weights
  iter <- 1; converged <- F
  if (verbose!="none") cat("cv of weights before raking =", sqrt(var(wtsRake$weights))/mean(wtsRake$weights), "\n")

  while (iter <= maxiter)
        {
        if (verbose=="full") cat("Iteration #" , iter, ":\n")
        for (i in 1:nbTargets)
            {
             # summing weights in each cell
             tmp <- split(wtsRake$weights, wtsRake[ , i])
             tmp <- sapply(tmp, sum)
             tmp <- data.frame(names(tmp), tmp, row.names=NULL)

             # computing the inflating/deflating factor for each cell
             factors <- merge(tmp, targets[[i]], by=1)
             factors <- data.frame(factors[ , 1], factors[ , namesTotals[[i]]] / factors[ , 2])
             colnames(factors) <- c(paste("cell", i, sep=""), "factors")

             # multiplying the weights by the above inflating/deflating factors
             wtsRake <- plyr::join(wtsRake, factors, by=paste("cell", i, sep=""))
             wtsRake$weights <-  wtsRake$weights * wtsRake[ , "factors"]
             wtsRake <- wtsRake[ , -ncol(wtsRake)]

             if (verbose=="full") cat("\t cv after calibrating on", paste(levels[[i]], collapse=" x "),
              		"=", sqrt(var(wtsRake$weights))/mean(wtsRake$weights), "\n")
            }
         #if (verbose=="some") cat("cv after iteration #", iter, "=", sqrt(var(wtsRake$weights))/mean(wtsRake$weights), "\n")

         # weight trimming
         if (trim)
         {
     	    # comptuting the upper bound(s)
       	if (is.null(trim.group)) {
       		if (trim.method=="median") bound <- median(wtsRake$weights) * trim.max else
       								   bound <- mean(wtsRake$weights) * trim.max
       		wtsRake <- cbind(wtsRake, bound)
       		} else {
      		     	tmp <- split(wtsRake$weights, wtsRake[ , trim.group]	)
        		    if (trim.method=="median") bounds <- sapply(tmp, median) * trim.max else
        		    							   bounds <- sapply(tmp, mean) * trim.max
      		     	bounds <- data.frame(bounds, names(bounds), row.names=NULL)
      		     	colnames(bounds) <- c("bound", trim.group)
      		     	wtsRake <- plyr::join(wtsRake, bounds, by=trim.group)
      		     	}

     		# saving the values of the weights before trimming
       	wtsRake$weights.ori <- wtsRake$weights

     		needed <- T; iter.trim <- 0
          while(iter.trim <= 10)
          {
             # checking if trimming is needed
  		   needed <- any(wtsRake$weights > wtsRake$bound)
              if (!needed)
                 {
                 if (verbose=="some") cat("Trimming for iteration	#" , iter, ":\n")
                 if (iter.trim==0) { if (verbose!="none") cat("\t No weights needed to be trimmed \n")
                 	  } else {
                 	  	     if (verbose!="none")
                 	  	     	{
                           	cat("\t After", iter.trim, "trimming cycle(s),", sum(wtsRake$trimmed), "observations had their weights trimmed \n")
              	   	     	if (verbose=="full") cat("\t")
              	   	     	cat("\t Sum of weights for those", sum(wtsRake$trimmed), "observations before trimming =",
              	   	     			 sum(wtsRake$weights.ori[wtsRake$trimmed==1]), "\n")
              	   	     	if (verbose=="full") cat("\t")
              	   	     	cat("\t Sum of weights for those", sum(wtsRake$trimmed), "observations after trimming =",
              	   	     			sum(wtsRake$weights[wtsRake$trimmed==1]), "\n")
              	   	     	}
              	   	     wtsRake$trimmed <- NULL; wtsRake$weights.ori <- NULL
              	   	     }
                  break
              	}

              # actual trimming
   			wtsRake$weights <- pmin(wtsRake$weights, wtsRake$bound)

     		    # identifying the weights that were trimmed
     		    wtsRake$trimmed <- ifelse(wtsRake$weights==wtsRake$bound, 1, 0)

     	        for (i in 1:nbTargets)
     	            {
   	             sumOriWts <- split(wtsRake$weights.ori, wtsRake[ , i])
                   sumOriWts <- sapply(sumOriWts, sum)
                   sumOriWts <- data.frame(names(sumOriWts), sumOriWts, row.names=NULL)

      	             sumNotTrim <- wtsRake[wtsRake$trimmed==0, ]
     	             sumNotTrim <- split(sumNotTrim$weights, sumNotTrim[ , i])
                   sumNotTrim <- sapply(sumNotTrim, sum)
                   sumNotTrim <- data.frame(names(sumNotTrim), sumNotTrim, row.names=NULL)

     	             sumTrim <- wtsRake[wtsRake$trimmed==1, ]
     	             sumTrim <- split(sumTrim$weights, as.factor(sumTrim[ , i]), drop=FALSE)
                   sumTrim <- sapply(sumTrim, sum)
                   sumTrim <- data.frame(names(sumTrim), sumTrim, row.names=NULL)
                   # these 3 lines were added on Jun 4, 2020 as the drop=F of slipt no longer works correctly
                   sumTrim = merge(sumTrim,  sumOriWts, by=1, all.y=T)
                   sumTrim <- sumTrim[ , -ncol(sumTrim)]
                   sumTrim[is.na(sumTrim[ , 2]), 2] <- 0

                   factorTrim <- merge(sumOriWts, merge(sumTrim, sumNotTrim, by=1), by=1)
                   factorTrim[ , ncol(factorTrim)+1] <- (factorTrim[ , 2] - factorTrim[ , 3]) / factorTrim[ , 4]
                   factorTrim <- data.frame(factorTrim[ , 1], factorTrim[, ncol(factorTrim)])
                   colnames(factorTrim) <- c(paste("cell", i, sep=""), "factors")

                   wtsRake <- plyr::join(wtsRake, factorTrim, by=paste("cell", i, sep=""))
                   wtsRake$weights <-  ifelse(wtsRake$trimmed==1, wtsRake$weights, wtsRake$weights * wtsRake[ , "factors"])
                   wtsRake <- wtsRake[ , -ncol(wtsRake)]
      	         }
      	   iter.trim <- iter.trim + 1
      	} }
      	if (verbose=="some") cat("cv after iteration #", iter, "=", sqrt(var(wtsRake$weights))/mean(wtsRake$weights), "\n")
     	   if (verbose=="some" & trim) cat("\n")

         # checking if converged
         diff <- NULL
         for (i in 1:nbTargets)
            {
             # summing weights in each cell
             tmp <- split(wtsRake$weights, wtsRake[ , i])
             tmp <- sapply(tmp, sum)
             tmp <- data.frame(names(tmp), tmp, row.names=NULL)

             # computing the differences between sum of weights and targets
             tmp <- merge(tmp, targets[[i]], by=1)
             diff[[i]] <- tmp[ , namesTotals[[i]]] - tmp[ , 2]
             }
         if (max(abs(unlist(diff))) <= epsilon) {converged <- TRUE; break}
         iter <- iter + 1
         if (verbose=="full") cat("\n")
        }

  # checks
  if (any(wtsRake$weights <= 0)) stop("One or more weights equal to or less than 0")
  if (any(is.na(wtsRake$weights))) stop("One or more weights with missing value")
  if (!converged) warning("The raking algorithm failed to converged within the maximum number of iterations allowed; see 'maxiter'")
  if (!converged && trim) warning("It's possible that the raking algorithm failed to converged because of weight trimming; try with 'trim=F' or with a larger value for 'trim.max'")

  # output
  list(weights=wtsRake$weights, converged=converged, nbIter=iter)
} # raking = function(...)







