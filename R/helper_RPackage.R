gini <- function (x, weights = rep(1, length = length(x))) 
{
    ox <- order(x)
    x <- x[ox]
    weights <- weights[ox]/sum(weights)
    p <- cumsum(weights)
    nu <- cumsum(weights * x)
    n <- length(nu)
    nu <- nu/nu[n]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

#' adjustCiteMetrics
#' 
#' Use resampling to correct citation inequality measures for marginals bias
#' 
#' @param papers required dataframe containing all citable papers studied (one row per paper). This dataframe must include two columns: the published paper IDs (\code{paperID}) and the publication year of the paper (\code{publishedYear}). 
#' @param cites required dataframe containing all citations made to the papers listed in \code{papers} over a given set of years (one row per citation). This dataframe must include at least two columns: the IDs of the cited papers (\code{citedPaperID}), corresponding to the IDs given in \code{papers}, and the cited paper's publication year (\code{citedPaperYear}). Optionally, this dataframe may include the published year of paper sending the citation (\code{citingYear}) may be included optional.  \code{citingYear} must be included if users wish to filter citations made within \code{citationWindow} years of the cited papers' publication dates.  
#' @param pubID character string, the column of \code{papers} containing publication IDs (default is `paperID').  IDs may be numeric, character, or factor.
#' @param pubYear character string, the column of \code{papers} containing publication years of each paper in numeric format (default is `publishedYear')
#' @param citedID character string, the column of \code{cites} containing the IDs of cited papers (default is `citedPaperID'). IDs may be numeric, character, or factor.
#' @param citedYear character string, the column of \code{cites} containing the years (as numeric data) in which papers were cited (default is `citedPaperYear')
#' @param citingYear character string, the column of \code{cites} containing the year (as numeric data) in which each citing paper was published (default = NULL)
#' @param citationWindow numeric scalar, the number of years covered by the citation window. For example, 2 means the function only analyzes citations made within 2 years after the cited paper is published. (default = NULL)
#' @param quantiles numeric vector in the interval \(0,1\] indicating the citation quantiles to be computed. Default = c(.20, .80), for the percentage of papers accounting for 20 percent and 80 percent of all citations.
#' @param refYear numeric scalar, the reference (or baseline) year among available publication years of papers in \code{papers}. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)
#' @param refPaperCount a numeric scalar identifying the referenced publication count. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)
#' @param refCiteCount a numeric scalar identifying the referenced citation count. Either \code{refYear} or \code{refPaperCount}/\code{refCiteCount} must be provided. (default = NULL)
#' @param sims optional; a numeric scalar identifying the number of times to run the simulation (default = 10)
#' @param lowCitesInclusion logical, TRUE: reports all results; FALSE: excludes any years from the results where the total number of citations is lower than the reference year in at least \code{lowCitesThreshold} percent of simulations (default=FALSE)
#' @param lowCitesThreshold numeric scalar in the interval \[0,1\] identifying threshold to define lowCites (default = 0.1)
#' @param paperDupThreshold numeric scalar in the interval \(0,1\] (default = 0.95)
#' @param periods optional; numeric vectors identifying which years we analyze in the data. If it is NULL (default), the function automatically selects all years in \code{papers}. 
#' @param verbose logical; returns detailed warnings related to lowCites (default=FALSE)
#'
#' @details
#' paragraph ABC
#' 
#' paragraph DEF
#' 
#' @return An object of class STM
#' \describe{
#'  \item{mu}{The corpus mean of topic prevalence and coefficients}
#'  \item{sigma}{Covariance matrix}
#' }
#' 
#' @references
#' ref1 \url{https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' ref2 \url{https://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' 
#' @seealso \code{\link{citeIneq}}
#'
#' @examples 
#' data(papers)
#' data(cites)
#' result <- adjustCiteMetrics(papers = papers, 
#' cites = cites, 
#' pubID = "paperID", 
#' pubYear = "publishedYear", 
#' citedID = "citedID",
#' citedYear = "citedYear",
#' refYear = 1996)
#' @export
adjustCiteMetrics <- function(papers,
                              cites,
                              pubID="paperID",
                              pubYear="publishedYear",
                              citedID="citedPaperID",
                              citedYear="citedPaperYear",
                              citingYear=NULL,
                              citationWindow=NULL,
                              quantiles = c(0.2, 0.8),
                              refYear=NULL, 
                              refPaperCount=NULL, refCiteCount=NULL,
                              sims=10, 
                              lowCitesInclusion=FALSE,
                              lowCitesThreshold=0.1,
                              paperDupThreshold=0.95, 
                              periods=NULL,
                              verbose=FALSE) {
    

    ## This simulator computes the Gini coefficient (Gini) for the
    ## observed citation distribution

    ## It also computes a "subsampled" Gini (GiniRS) by randomly sampling
    ## nR0 papers from the original set of citable papers, then randomly
    ## samples nC0*nS0 citations to those papers from the total set of
    ## citations to the papers in the nR0 subsample.

    ## In this way, GiniRS attempts to correct for any effects of changing
    ## marginal numbers of citations from sending papers to a set of receiving
    ## papers by imposing the same marginals on each of these three quantities 
        
    #######
    ## Error checking for inputs
    
    if (!is.data.frame(papers)) {
        stop("papers must be a data frame with columns containing the publication ID and publication year")
        }
    
    if (((length(pubID)!=1)||!is.character(pubID))) {
        stop("pubID must be a character string identifying the column of papers containing publication IDs")
        }
    
    if ((length(pubYear)!=1)||(!is.character(pubYear))) {
        stop("pubYear must be a character string identifying the column of papers containing publication years")
        }
    
    if (!is.data.frame(cites)) {
        stop("cites must be a data frame with columns containing the IDs of cited papers and years in which they were cited by other papers")
        }
    
    if ((length(citedID)!=1)||(!is.character(citedID))) {
        stop("citedID must be a character string identifying the column of cites containing the IDs of cited papers") 
        }
    
    if ((length(citedYear)!=1)||(!is.character(citedYear))) {
        stop("citedYear must be a character string identifying the column of cites containing the years in which papers were cited") 
        }

    if (!is.null(citingYear)) {
        if ( (length(citingYear)!=1) || (!is.character(citingYear))) {
            stop("citingYear must be NULL or a character string identifying the column of cites containing the year in which each citing paper was published") 
        }
    }
    
    if (!is.null(citationWindow)) {
        if ((length(citationWindow)!=1) || (!is.numeric(citationWindow))) {
            stop("citationWindow must be a numeric scalar or NULL (the default)")
        }
    }

    sims <- as.integer(sims)
    if ((length(sims)!=1)||(sims<=0)) {
        stop("sims must be a positive integer")
    }
    
    if ((!is.numeric(paperDupThreshold)) || (length(paperDupThreshold)!=1) || (paperDupThreshold>1) || (paperDupThreshold<=0)) {
        stop("paperDupThreshold should be a numeric scalar in the interval (0,1]")
    }

    if ((!is.numeric(quantiles)) || any(quantiles>=1) || any(quantiles<=0)) {
        stop("quantiles should be numeric in the interval (0,1)")
    }
    ## Convert quantiles to percent
    quantiles <- 100*quantiles

    ## Create names of quantiles for output
    qInt <- floor(quantiles)
    qDec <- (quantiles)%%1
    qStr <- sprintf('q%02d', qInt)
    for (i in 1:length(qStr)) {
        if (qDec[i]>0) {
            # qStr[i] <- paste0(qStr[i],
            #                   substr(as.character(qDec[i]), 2,
            #                          nchar(as.character(qDec[i]))))
            qStr[i] <- paste0("q", gsub("\\.*0+$", "", as.character(quantiles[i]), perl=TRUE))
        }
    }

    if (!is.null(refYear)) {
        if ((!is.numeric(refYear)) || (length(refYear)!=1)) {
            stop("refYear should be a numeric scalar identifying the referenced year among available pubYear in papers")
        }
    }
    
    if (!is.null(refPaperCount)) {
        if ( (length(refPaperCount)!=1) || (!is.numeric(refPaperCount)) || (refCiteCount<=0)) {
            stop("refPaperCount must be NULL or a numeric scalar identifying the referenced publication count") 
        }
    }
    
    if (!is.null(refCiteCount)) {
        if ( (length(refCiteCount)!=1) || (!is.numeric(refCiteCount)) || (refCiteCount<=0)) {
            stop("refCiteCount must be NULL or a positive numeric scalar identifying the referenced citation count") 
        }
    }
    
    if (!is.logical(lowCitesInclusion)) {
        stop("lowCitesInclusion must be logical")
    }

    if ((!is.numeric(lowCitesThreshold)) || (lowCitesThreshold>1) || (lowCitesThreshold<0)) {
        stop("lowCitesThreshold should be numeric in the interval [0,1]")
    }
    
    if (!is.null(periods) & !is.numeric(periods)) {
        stop("periods should be numeric vector or NULL (the default)")
    }

    if (!is.logical(verbose)) {
        stop("verbose must be logical")
    }

    ## Set up dataframes
    ## papers
    if(pubID %in% colnames(papers)) {
        papers$paperID <- papers[,pubID]
    } else stop(paste("papers does not contain a column named", pubID, ": pubID should refer to a variable in papers"))
    
    if(pubYear %in% colnames(papers)) {
        papers$publishedYear <- papers[,pubYear]
    } else stop(paste("papers does not contain a column named", pubYear, ": pubYear should refer to a variable in papers"))
    
    ## cites
    if(citedID %in% colnames(cites)) {
        cites$citedPaperID <- cites[,citedID]
    } else stop(paste("cites does not contain a column named", citedID, ": citedID should refer to a variable in cites"))
    
    if(citedYear %in% colnames(cites)) {
        cites$citedPaperYear <- cites[,citedYear]
    } else stop(paste("cites does not contain a column named", citedYear, ": citedYear should refer to a variable in cites"))
    
    if (is.null(citingYear)) {
        cites$citingPaperYear <- cites[,citingYear]
    } else if (citingYear %in% colnames(cites)) {
        cites$citingPaperYear <- cites[,citingYear]
    } else stop(paste("cites does not contain a column named", citingYear, ": citingYear should refer to a variable in cites or be left NULL"))
   
     
    ## Filter out "citedPaperID" that are not included in papers$paperID. 
    citesBeforeFilterLength <- length(cites$citedPaperID)
    cites <- cites[cites$citedPaperID %in% papers$paperID, ]
    citesAfterFilterLength <- length(cites$citedPaperID)
    
    if (citesBeforeFilterLength != citesAfterFilterLength) {
        warning("Some citation information in the cites dataframe have been deleted because the cited papers do not exist in papers")
    }
    
    ## Filter cites in a given citation window.
    if (is.null(cites$citingPaperYear) && !is.null(citationWindow)) {
            stop("Citing paper's published year needs to be provided to apply the citation window")
        } else if (length(cites$citingPaperYear)>0) {
            cites <- cites[(cites$citingPaperYear - cites$citedPaperYear) <= citationWindow,]
            } else cites <- cites
                
        
    ## Reference year / Referenced marginals
    ## Set up nR0 and nS0XnC0
    if (!is.null(refYear)) {
        if (!(refYear %in% unique(papers$publishedYear))) {
            stop("refYear must be included in pubYears in papers ")
        }
        
        if (!is.null(refPaperCount)) {
            stop("either refYear or refPaperCount must be provided.")
        }
        
        if (!is.null(refCiteCount)) {
            stop("either refYear or refCiteCount must be provided.")
        }
        
        nR0 <- nrow(papers[papers$publishedYear==refYear,])
        nS0XnC0 <- nrow(cites[cites$citedPaperYear==refYear,])
        } 
    
    if (is.null(refYear)) {
        if (is.null(refPaperCount)) {
            stop("Either refYear or refPaperCount needs to be provided.")
        } 
        
        if (is.null(refCiteCount)) {
            stop("Either refYear or refCiteCount needs to be provided.")
        }
        
        nR0 <- refPaperCount
        nS0XnC0 <- refCiteCount
    }

    ## Set up Years
    if (is.null(periods)) {
        years <- unique(papers$publishedYear)
    } else if (sum(periods %in% unique(papers$publishedYear))==length(periods)) {
        years <- periods
    } else stop("periods should exist in papers$pubYear")
    
    ## Force to integer
    nR0 <- as.integer(nR0)
    nS0XnC0 <- as.integer(nS0XnC0)
    
    resALL <- NULL
    
    for (i in 1:length(years)) {
        
        if (is.null(refYear)) {
            reference <- FALSE 
            } else if (years[i]==refYear) reference <- TRUE else reference <- FALSE
        
        ## Set up papers and cites by year
        papersYear <- papers[papers$publishedYear==years[i],]$paperID
        citesYear <- cites[cites$citedPaperYear==years[i],]$citedPaperID
        
        ## Sum of citable papers
        nR <- length(unique(papersYear))
        nSXnC <- length(citesYear)
        
        ## Set up vectors to store the distribution of citations averaged
        ##  across simulation runs for both the "true" and "subsampled" citations
        cumDistRS <- rep(0, nR0)
        
        ## Set up vectors to hold the simulated Ginis
        giniCrs <- herfRS <- withCitesRS <- rep(NA, sims)
        
        pctListRS <- as.data.frame(matrix(nrow=sims, ncol=length(quantiles)))
        
        for (q in 1:length(quantiles)) {
            colnames(pctListRS)[q] <- paste("pct", quantiles[q], "RS", sep="")
        }
        
        lowCites <- rep(0, sims)
        
        ## Loop over simulation runs to reduce Monte Carlo error
        for (j in 1:sims) {
            
            ## Sample id's of subsample of citable papers
            if (!reference&&(nR0>(paperDupThreshold*nR))) {
                duplicated <- TRUE
                papersDup <- papersYear
                citesDup <- citesYear
                ctr <- 0
                done <- FALSE
                while (!done) {
                    ctr <- ctr +1
                    papers0 <- paste0(papersYear, "dup", ctr)
                    papersDup <- c(papersDup, papers0)
                    cites0 <- paste0(citesYear, "dup", ctr)
                    citesDup <- c(citesDup, cites0)
                    if (length(papersDup)>(nR0*(1+paperDupThreshold))) done <- TRUE
                }
                papersRS <- sample(papersDup, nR0, replace=FALSE)
            } else {
                duplicated <- FALSE
                papersRS <- sample(papersYear, nR0, replace=FALSE)
            }    
            
            ## Set up placeholders for citation distributions
            citeDistRS <- rep(0, nR0)
            
            ## Find eligible cites
            if (duplicated) {
                eligibleCites <- citesDup[citesDup %in% papersRS]
            } else {
                eligibleCites <- citesYear[citesYear %in% papersRS]
            }
            
            ## Subsample eligible cites
            if (nS0XnC0<=length(eligibleCites)) {
                citedPapersRS <- sample(eligibleCites, nS0XnC0, replace=FALSE)
            } else {
                if (verbose) {
                    if (lowCitesInclusion) {
                        warnMsg <- "Using all available cites."
                    } else {
                        warnMsg <- paste("Year", years[i], "reported as NA for corrected inequality metrics.")
                    }
                    warning(paste0("Not enough citations for complete subsample in year ", years[i],".  Needed ",
                                   nS0XnC0,
                                   "; had ",
                                   length(eligibleCites),
                                   ". ", warnMsg))
                }
                citedPapersRS <- eligibleCites
                lowCites[j] <- 1
            }
            
            
            ## Tabulate cites papers from resampling adjustment
            citeCountsRS <- rep(0, nR0)
            citeCountsRS0 <- rev(sort(as.vector(table(citedPapersRS))))
            citeCountsRS[1:length(citeCountsRS0)] <- citeCountsRS0
            cumDistRS <- cumDistRS + citeCountsRS
            
            ## Compute Gini for subsample
            giniCrs[j] <- gini(citeCountsRS)
            
            ## Compute % with citations for subsample
            withCitesRS[j] <- sum(citeCountsRS>0)/nR0
            
            ## Compute smallest % of papers with 20%, 80% of cites for subsamp
            cumpctRS <- 100*cumsum(citeCountsRS)/sum(citeCountsRS)
            
            for (q in 1:length(quantiles)) {
                pctListRS[j,q] <- 100*(sum(cumpctRS<quantiles[q])+1)/nR0    
            }
            
            ## Compute Herfindahl-Hirschmann index for subsample
            sharesRS <- citeCountsRS/sum(citeCountsRS)
            herfRS[j] <- sum(sharesRS^2)
            
        }
        
        
        ## Tabulate cited papers for "true" distibuiton
        citeCounts <- rep(0, nR)
        citeCounts0 <- rev(sort(as.vector(table(citesYear))))
        citeCounts[1:length(citeCounts0)] <- citeCounts0
        
        ## Compute Gini for "true" distribution
        giniC <- gini(citeCounts)
        
        ## Compute % with citations for "true" distribution
        withCites <- sum(citeCounts>0)/nR
        
        ## Compute smallest % of papers with 20%, 80% of cites
        cumpct <- 100*cumsum(citeCounts)/sum(citeCounts)
        
        pctList <- as.data.frame(matrix(nrow=1, ncol=length(quantiles)))
        
        for (q in 1:length(quantiles)) {
            colnames(pctList)[q] <- paste("pct", quantiles[q], sep="")
            pctList[1,q] <- 100*(sum(cumpct<quantiles[q])+1)/nR
        }
        
        ## Compute Herfindahl-Hirschmann index
        shares <- citeCounts/sum(citeCounts)
        herf <- sum(shares^2)
        
        ## Compute low cites percentage
        lowCites <- mean(lowCites)
        
        ## Collect output
        res <- list(giniUC=giniC,
                    everCitedUC=withCites,
                    hhiUC=herf,
                    giniRS=mean(giniCrs),
                    everCitedRS=mean(withCitesRS),
                    hhiRS=mean(herfRS),
                    nR=nR,
                    nR0=nR0,
                    nSXnC=nSXnC,
                    nS0XnC0=nS0XnC0,
                    lowCites=lowCites
                    )

        ## Collect names of all metrics
        availMetrics <- c("everCited", "gini", "hhi", qStr)

        ## Add quantiles to output
        for (q in 1:length(quantiles)) {
            res <- c(res, assign(paste0(qStr[q], "UC"), pctList[1,q]))
            names(res)[length(res)] <- paste0(qStr[q], "UC")
        }
        
        for (q in 1:length(quantiles)) {
            res <- c(res, assign(paste0(qStr[q], "RS"), mean(pctListRS[,q])))
            names(res)[length(res)] <- paste0(qStr[q], "RS")
        }
        
        resALL <- rbind(resALL, c(years[i], unlist(res)))
    }
    
    resALL <- as.data.frame(resALL)
    names(resALL)[1] <- "year"
    
    ## lowCites option: if TRUE (default), deletes lowCites==1; if FALSE, presents all but gives a warning.
    lowCitesYears <- resALL$year[resALL$lowCites>lowCitesThreshold]
    if (isFALSE(lowCitesInclusion)) {
        resALL[resALL$lowCites>lowCitesThreshold, grepl("RS", names(resALL))] <- NA
        warning("Resampled inequality measures for the years ",  paste(lowCitesYears, collapse=", "), " recorded as NA due to smaller number of publication and citation count than the referenced year.  Consider using a different base year or providing a lower refCiteCount. Set verbose=TRUE to see more details.")
    } else {
        warning(paste("Results for the years ", paste(lowCitesYears, collapse=", "), " may be incorrect due to smaller number of publication and citation count than the referenced year. Consider using a different base year or providing a lower refCiteCount. Set verbose=TRUE to see more details."))
    }
    
    ## Assign class to output
    class(resALL)  <- "adjustCiteMetrics"

    ## Add available metrics
    resALL$availMetrics <- availMetrics
    
    return(resALL)
}

#' citeIneq
#' 
#' Extract citation inequality metrics from an \code{adjustCiteMetrics} object after resampling to correct for marginals bias
#' 
#' @param result a \code{adjustCiteMetrics} class object created by \code{adjustCiteMetrics}.
#' @param metric character vector indicating the inequality measures to report; possible choices are `everCited', `gini', `hhi', or user-requested quantiles, such as `q20' or `q80'.  Default is `all', which reports all metrics computed by \code{adjustCiteMetrics}. 
#' @param type character string, the default, `resampled', reports the adjusted inequality measures.  Set to `uncorrected' to see the inequality metrics before sampling adjustment is applied
#' @param showMargins logical, report information of the marginal and reference number of papers and citations (default is FALSE) 
#' 
#' @examples 
#' citeIneq(result) 
citeIneq <- function(result, 
                     metric="all", 
                     type="resampled", 
                     showMargins=FALSE
                     ) {

    ## check that result if an adjustCiteMetrics class
    if (class(result)!="adjustCiteMetrics")
        stop("result must be of class adjustCiteMetrics")

    ## Check for possible type
    if ((type!="resampled")&&(type!="uncorrected"))
        stop("type must be `resampled' or `uncorrected'") 

    ## Check for all metrics
    if (identical(metric, "all")) {
        metric <- result$availMetrics
    }
    
    ## Check called metrics exist in adjustCiteMetrics class
    if (length(setdiff(metric, result$availMetrics))) {
        stop("metric can only request inequality measures reported in the supplied adjustCiteMetrics object")
    }
    
    ## Combination of metric and type
    if ("resampled" %in% type) {
        metricRS <- paste0(metric, "RS")
    } else metricRS <- NULL
    
    if ("uncorrected" %in% type) {
        metricUC <- paste0(metric, "UC")
    } else metricUC <- NULL
    
    ## Apply margins option
    if(isTRUE(showMargins)) {
        marginsList <- c("nR", "nR0", "nSXnC", "nS0XnC0")
    } else marginsList <- NULL
    
    metric <- c("year", metricRS, metricUC, marginsList)
    
    result <- result[metric] 
    result <- as.data.frame(result)
    result <- result[order(result$year),]
    rownames(result) <- 1:nrow(result)
    
    return(result)
}
