

# functions
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}


#' Calculate test statistic for ANCOM
#' 
#' @description
#' Calculates test statistics for differences in OTU abundances between treatment groups.
#' 
#' @param otu_data the OTU dataset.
#' @param n_otu the number of OTUs.
#' @param alpha the significance level at which the tests are to be performed.
#' @param multcorr type of correction for multiple comparisons, see Details.
#' @param Wexact logical, should Wilcoxon tests return exact p-values?
#' @param ncore if ncore>1, then \pkg{doParallel} will be loaded and used.
#' 
#' @details
#' \code{multcorr} can take on values of 1 (no correction), 2 (a less stringent)
#' correction, or 3 (a more stringent correction).
#' 
#' @note
#' This function is intended to be called by \code{\link{ANCOM}}, see the documentation of
#' that function for details on using the method.
#' 
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @importFrom exactRankTests wilcox.exact
#' @importFrom coin kruskal_test
#' @importMethodsFrom coin pvalue
#' @export
#' 
#' 

ancom.detect <- function(otu_data, n_otu, alpha, multcorr, ncore){
  
  ## Detect whether the data are dependent or not
  if( ncol(otu_data) == n_otu+1  ){
    Group     <- otu_data[, ncol(otu_data) ]
    ID        <- rep( 1 , nrow(otu_data) )
    repeated <- FALSE
    fformula  <- formula("lr ~ Group")
  } else if( ncol(otu_data) == n_otu+2  ){
    Group     <- otu_data[, ncol(otu_data)-1 ]
    ID        <- otu_data[, ncol(otu_data)   ]
    repeated <- TRUE
    fformula  <- formula("lr ~ Group | ID")
    
  } else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  ## Detect which test to use, 
  ## Dependent data: Friedman test
  ## Independent data: Wilcoxon Rank-Sum or the Kruskal-Wallis
  ## exactRankTests::wilcox.exact is faster than stats::kruskal.test and stats::wilcox.test
  if( repeated==FALSE ){
    if( length(unique(Group))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  } else{
    tfun <- stats::friedman.test
  }
  
  
  
  ## Parallelized way to get the logratio.mat
  ## Doubles the number of computations to make, so only run the parallel
  ## version if there are multiple cores. Method may also add some computational
  ## overhead, so if only 2 cores, the nested for-loop shoud have advantage
  ## over the parallel loop (though I have not tested that).
  ## For some reason this is taking much longer, do not run the parallel loop as of now.
  if( FALSE ){
    registerDoParallel( cores=ncore )
    
    aa <- bb <- NULL
    logratio.mat <- foreach( bb = 1:n_otu, .combine='rbind', .packages="foreach" ) %:% 
      foreach( aa = 1:n_otu , .combine='c',  .packages="foreach" ) %dopar% {
        if( aa==bb ){
          p_out <- NA
        } else{
          data.pair <- otu_data[,c(aa,bb)]
          lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
          lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
          p_out  <- tfun(formula=fformula, data = lr_dat)$p.value
        }
        p_out
      }
    rownames(logratio.mat) <- colnames(logratio.mat) <- NULL
  } else{
    logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
    for(ii in 1:(n_otu-1)){
      for(jj in (ii+1):n_otu){
        data.pair <- otu_data[,c(ii,jj)]
        lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
        lr_dat <- data.frame( lr=lr, Group=Group, ID=ID )
        
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }
    }  
    ind <- lower.tri(logratio.mat)
    logratio.mat[ind] <- t(logratio.mat)[ind]
  }
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<alpha))
    })
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<alpha))
    })
  }
  return(W)
  }

############################################################
############################################################


#' Run the ANCOM method
#' 
#' @description
#' Runs ANCOM to test for differences in OTU abundances between treatment groups.
#' 
#' @param OTUdat the OTU dataset. See Details for formatting instructions.
#' @param sig the significance level (or FDR) at which the tests are to be performed.
#' @param multcorr type of correction for multiple comparisons, see Details.
#' @param tau a tuning parameter in the Stepwise testing method. See Details.
#' @param theta a tuning parameter in the Stepwise testing method. See Details.
#' @param repeated logical determining whether the data have repeated measures (e.g., longitudinal design).
#' @details
#' The ANCOM method was developed and tested with default values of the two tuning parameters 
#' (\code{tau=0.02} and \code{theta=0.1}). For consistency, users are recommended to leave 
#' these tuning parameters at their default values, unless they wish to explore the performance
#'  of ANCOM for different values of the tuning parameters.
#' 
#' Data should be formatted as follows: each row is a subject, and each column is an OTU.
#' The final column should contain the grouping variable.
#' 
#' To adjust for multiple testing, \code{multcorr} may take on the following values:
#' \itemize{
#' \item{ \code{1}: }{ A stringent correction}
#' \item{ \code{2}: }{ A less stringent correction}
#' \item{ \code{3}: }{ No correction (default)}
#' }
#' The more stringent correction is not available in the shiny application.
#' 
#' @note
#' The function \code{\link{plot_ancom}} will produce plots for objects produced by \code{ANCOM}.
#' 
#' @return
#' The function produces a list with the following elements:
#' \itemize{
#' \item{ \code{W}: }{ values of the test statistics.}
#' \item{ \code{detected}: }{ names of OTUs detected.}
#' \item{ \code{dframe}: }{ the input dataframe.}
#' }
#'  
#' @export
#' 
#' @examples
#' 
#' \dontrun{
#' ## Create and run a small example
#' 
#' nn <- 10
#' pp <- 20
#' sim_otu <- matrix( 0, nrow=nn, ncol=pp+1 )
#' sim_otu <- data.frame(sim_otu)
#' colnames(sim_otu) <- c( paste0("OTU_", letters[1:pp] ), "Group" )
#' sim_otu[,pp+1]    <- c( rep("Control",nn/2), rep("Treatment",nn/2)  )
#' idx_trt <- sim_otu$Group=="Treatment"
#' 
#' for( ii in 1:pp ){
#'   sim_otu[,ii] <- rpois( nn, 1 )
#' }
#' 
#' # Create some significance
#' sim_otu[idx_trt,3] <- rpois( nn/2, 8)
#' sim_otu[idx_trt,7] <- rpois( nn/2, 8)
#' sim_otu[idx_trt,9] <- rpois( nn/2, 8)
#' 
#' ancom.out <- ANCOM( OTUdat = sim_otu, sig = 0.20, multcorr = 2 )
#' ancom.out$W
#' ancom.out$detected
#' }
#' 

ANCOM  <-  function(OTUdat, sig=0.05, multcorr=3, tau=0.02, theta=0.1, repeated=FALSE ){
  
  #OTUdat <-  read.delim(filepath,header=TRUE)
  
  num_col <- ncol( OTUdat )
  
  if( repeated==FALSE ){
    colnames(OTUdat)[ num_col ] <- "Group"    # rename last column as "Group"
    num_OTU      <- ncol(OTUdat) - 1
    
    sub_drop <- data.frame( nm_drop= "N/A" )
    sub_keep <- data.frame( nm_keep= "All subjects" )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "No subjects entirely removed (not a repeated-measures design)" )
    
  } else{
    colnames(OTUdat)[ num_col-1 ] <- "Group"    # rename 2nd last column as "Group"
    colnames(OTUdat)[ num_col   ] <- "ID"    # rename last column as "ID"
    OTUdat$ID    <- factor( OTUdat$ID )
    num_OTU      <- ncol(OTUdat) - 2
    
    ## Drop subjects if missing at a given time point
    crossTab <- table( OTUdat$Group , OTUdat$ID  )==0
    id_drop  <- apply( crossTab, 2, FUN=function(x) any(x)  )
    nm_drop  <- names( which( id_drop ) )
    idx_drop <- OTUdat$ID %in% nm_drop
    OTUdat   <- OTUdat[ idx_drop==FALSE, ]
    
    if( nrow(OTUdat)==0 ){ stop("Too many missing values in data, all subjects dropped") }    
    OTUdat$ID <- droplevels( OTUdat$ID )    
    num_dropped <- sum(id_drop)
    num_retain  <- length(id_drop) - num_dropped
    
    sub_drop <- data.frame( nm_drop=paste(nm_drop, collapse=", " ) )
    sub_keep <- data.frame( nm_keep= paste(levels(OTUdat$ID), collapse=", " ) )
    colnames(sub_drop) <- "Subjects removed"
    colnames(sub_keep) <- "Subjects retained"
    n_summary   <- paste0( "Analysis used ", num_retain, " subjects (", num_dropped, " were removed due to incomplete data)")
  }
  
  OTUdat$Group <- factor( OTUdat$Group )
  OTUdat       <- data.frame( OTUdat[ which(is.na(OTUdat$Group)==FALSE),],row.names=NULL )
  
  W.detected   <- ancom.detect(OTUdat, num_OTU, sig, multcorr, ncore=1 )
  W_stat       <- W.detected
  
  
  ## Per Shyamal (June 4, 2015): 
  ##  If number of OTUs is < 10, then use 'arbitrary' method, reject Ho if 
  ##  W > p - 1, where p = number of OTUs. If number of OTUs > 10, then use
  ##  the stepwise method. Consequently, only one output will be produced,
  ##  instead of detected_arbitrary and detected_stepwise, produce "detected"
  ## Rephrase "arbitrary", since it's not arbitrary, more of an empirical method.
  
  ## Detected using arbitrary cutoff
  # Previous code:
  # detected_arbitrary <- colnames(OTUdat)[ which( W.detected > num_OTU*theta ) ]
  
  if( num_OTU < 10 ){
    detected <- colnames(OTUdat)[which(W.detected > num_OTU-1 )]    
  } else{
    ## Detected using a stepwise mode detection
    if( max(W.detected)/num_OTU >= theta ){
      c.start <- max(W.detected)/num_OTU
      cutoff  <- c.start-c(0.05,0.10,0.15,0.20,0.25)
      
      prop_cut <- rep(0,length(cutoff))
      for(cut in 1:length(cutoff)){
        prop_cut[cut] <- length(which(W.detected>=num_OTU*cutoff[cut]))/length(W.detected)
      } 
      
      del <- rep(0,length(cutoff)-1)
      for( ii in 1:(length(cutoff)-1) ){
        del[ii] <- abs(prop_cut[ii]-prop_cut[ii+1])
      }
      
      if(       del[1]< tau & del[2]<tau & del[3]<tau ){ nu=cutoff[1]
      }else if( del[1]>=tau & del[2]<tau & del[3]<tau ){ nu=cutoff[2]
      }else if( del[2]>=tau & del[3]<tau & del[4]<tau ){ nu=cutoff[3]                                
      }else{ nu=cutoff[4] }
      
      up_point <- min(W.detected[ which( W.detected >= nu*num_OTU ) ])
      
      W.detected[W.detected>=up_point] <- 99999
      W.detected[W.detected<up_point]  <- 0
      W.detected[W.detected==99999]    <- 1
      
      detected <- colnames(OTUdat)[which(W.detected==1)]
      
    } else{
      W.detected <- 0
      detected   <- "No significant OTUs detected"
    }
    
  }
  
  #results_list <- list( W         = W_stat,
  #                      Arbitrary = detected_arbitrary,
  #                      Stepwise  = detected_stepwise )
  #idx0 <- lapply( results_list , FUN=length)
  #results_list[idx0==0] <- "No significant OTUs detected"
  
  results <- list( W=W_stat, detected=detected, dframe=OTUdat, repeated=repeated,
                   n_summary=n_summary, sub_drop=sub_drop, sub_keep=sub_keep)
  class(results) <- "ancom"
  
  return(results)  
  
}
#########################################################################


#' Plot of data from objects of class 'ancom'
#' 
#' @description
#' Produces comparison boxplots of data for objects of class \code{ancom}.
#' 
#' @param object object of class \code{ancom}.
#' @param ncols the number of columns for \code{ggplot} to produce using \code{facet_wrap}.
#'        If \code{ncol=-1}, then the function will attempt to 
#' @param ... space for additional arguments (none corrently)
#' 
#' @details
#' \code{plot_ancom} uses \pkg{ggplot} to produce graphics.
#' 
#' 
#' @import ggplot2
#' 
#' @export
#' 


plot_ancom <- function( object, ncols=-1, ... ){
  
  if( !(class(object)=="ancom") ){
    stop("'object' is not of class ancom")
  }
  
  repeated <- object$repeated
  
  Group  <- OTU <- ID <- NULL
  dframe <- object$dframe
  
  if( repeated==FALSE ){
    colnames(dframe)[ ncol(dframe) ] <- "Group"  
  } else{
    colnames(dframe)[ ncol(dframe)-1 ] <- "Group"
  }
  
  # OTUs that were detected
  Sig_OTU <- object$detected
  
  if( Sig_OTU[1] == "No significant OTUs detected" ){
    ## Plot a message so that SOMETHING is produced
    plot( 1:5 , 1:5 , col="white", xaxt='n', yaxt='n', xlab="", ylab="", frame.plot=FALSE )
    text( 3, 3, labels="No significant OTUs detected" )
    
  } else{
    
    if( repeated==FALSE ){
      W_check <- data.frame( colnames(dframe)[-ncol(dframe)], object$W, row.names=NULL)
      colnames(W_check) <- c("OTU_ID","W")  
    } else{
      W_check <- data.frame( colnames(dframe)[-c(ncol(dframe)-1, ncol(dframe) )], object$W, row.names=NULL)
      colnames(W_check) <- c("OTU_ID","W")  
    }
    
    W_check        <- W_check[which(  W_check$OTU_ID %in% Sig_OTU),]
    W_check        <- W_check[order(-W_check$W),]
    nplot          <- nrow(W_check)
    W_check$OTU_ID <- as.character(W_check$OTU_ID)
    
    # Get the DATA to plot
    for( ii in 1:nplot ){
      dsub <- dframe[ , colnames(dframe) %in% c(W_check$OTU_ID[ii],"Group") ]
      colnames(dsub)=c("OTU","Group")
      OTU_name <- rep( W_check$OTU_ID[ii] , nrow(dsub) )
      pltDat00 <- data.frame( OTU_name , dsub )              
      if( ii==1 ){
        pltDat <- pltDat00
      } else{
        pltDat <- rbind(pltDat, pltDat00 )
      } 
    }
    
    pltDat$OTU      <- log( pltDat$OTU + 1 )
    pltDat$OTU_name <- factor( pltDat$OTU_name , Sig_OTU )
    
    if( ncols<1 ){
      ncols <- min(3, nplot)
    }
    
    gplot <- ggplot( pltDat , aes(x=factor(Group), y=OTU ) ) + 
      facet_wrap( ~ OTU_name , ncol=ncols, scales="free_y") + 
      geom_boxplot()
    
    gplot + labs(x = "Grouping Factor" , y="Log of Abundance"  ) + theme(
      panel.background  = element_rect( fill="white" , colour="black"),
      panel.grid        = element_blank(),
      strip.text        = element_text( size=rel(1.25)),
      strip.background  = element_rect( fill="grey90" , color="black"),
      axis.title        = element_text( size=rel(1.25) , color="black"),
      axis.text         = element_text( size=rel(1.05) , color="black"),
      # legend.position   = c(0.90,0.10),
      legend.position   = "none",
      legend.background = element_rect(colour = "black", fill="white"),
      legend.key        = element_rect(fill = "white"),
      legend.title      = element_text(size = 15 ),
      legend.text       = element_text(size = 12 )
    )
    
  }
  
}


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

psmelt.fast <- function(physeq){
  # Access covariate names from object, if present
  if(!inherits(physeq, "phyloseq")){
    rankNames = NULL
    sampleVars = NULL
  } else {
    # Still might be NULL, but attempt access
    rankNames = rank_names(physeq, FALSE)
    sampleVars = sample_variables(physeq, FALSE) 
  }
  # Define reserved names
  reservedVarnames = c("Sample", "Abundance", "OTU")  
  # type-1a conflict: between sample_data 
  # and reserved psmelt variable names
  type1aconflict = intersect(reservedVarnames, sampleVars)
  if(length(type1aconflict) > 0){
    wh1a = which(sampleVars %in% type1aconflict)
    new1a = paste0("sample_", sampleVars[wh1a])
    # First warn about the change
    warning("The sample variables: \n",
            paste(sampleVars[wh1a], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1a, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the sample variables.
    colnames(sample_data(physeq))[wh1a] <- new1a
  }
  # type-1b conflict: between tax_table
  # and reserved psmelt variable names
  type1bconflict = intersect(reservedVarnames, rankNames)
  if(length(type1bconflict) > 0){
    wh1b = which(rankNames %in% type1bconflict)
    new1b = paste0("taxa_", rankNames[wh1b])
    # First warn about the change
    warning("The rank names: \n",
            paste(rankNames[wh1b], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1b, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the conflicting taxonomic ranks
    colnames(tax_table(physeq))[wh1b] <- new1b
  }
  # type-2 conflict: internal between tax_table and sample_data
  type2conflict = intersect(sampleVars, rankNames)
  if(length(type2conflict) > 0){
    wh2 = which(sampleVars %in% type2conflict)
    new2 = paste0("sample_", sampleVars[wh2])
    # First warn about the change
    warning("The sample variables: \n",
            paste0(sampleVars[wh2], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new2, collapse=", "), "\n",
            "to avoid conflicts with taxonomic rank names.")
    # Rename the sample variables
    colnames(sample_data(physeq))[wh2] <- new2
  }
  # Enforce OTU table orientation. Redundant-looking step
  # supports "naked" otu_table as `physeq` input.
  otutab = otu_table(physeq)
  if(!taxa_are_rows(otutab)){otutab <- t(otutab)}
  ## Speedyseq specific code starts here
  # Convert the otu table to a tibble in tall form (one sample-taxon obsevation
  # per row)
  tb <- otutab %>% 
    as("matrix") %>%
    data.table::as.data.table(keep.rownames = "OTU") %>%
    data.table::melt(id.vars = c("OTU"), variable.name = "Sample", 
                     value.name = "Abundance")
  # Add the sample data if it exists
  if (!is.null(sampleVars)) {
    sam <- sample_data(physeq) %>%
      as("data.frame") %>% 
      data.table::as.data.table(keep.rownames = "Sample")
    tb <- tb[sam, on = .(Sample = Sample)]
  }
  # Add the tax table if it exists
  if (!is.null(rankNames)) {
    tax <- tax_table(physeq) %>%
      as("matrix") %>%
      as.data.frame %>%
      data.table::as.data.table(keep.rownames = "OTU")
    # NOTE: Conversion to data.frame functions to converts taxonomy vars to
    # factors if stringsAsFactors = TRUE, for phyloseq compatibility.
    tb <- tb[tax, on = .(OTU = OTU)]
  }
  # Arrange by Abundance, then OTU names (to approx. phyloseq behavior)
  tb <- tb %>%
    data.table::setorder(-Abundance, OTU)
  # Return as a data.frame for phyloseq compatibility
  tb %>% as.data.frame
}

fungicide.adonis <- function(phyloseq.object){
  distance.matrix <- phyloseq::distance(phyloseq.object, "bray") # create bray-curtis distance matrix
  return(adonis2(distance.matrix~Fungicide, as(sample_data(phyloseq.object), "data.frame"), permutations = 9999)) 
}

fungicide.betadisper <- function(phyloseq.object){
  distance.matrix <- phyloseq::distance(phyloseq.object, "bray") # create bray-curtis distance matrix
  betadisp <- betadisper(distance.matrix, phyloseq.object@sam_data$Fungicide)
  test.beta <- permutest(betadisp)
  return.list <- list(betadisp, test.beta)
  return(return.list)
}

fungicide.anosim <- function(phyloseq.object){
  distance.matrix <- phyloseq::distance(phyloseq.object, "bray")
  return(anosim(distance.matrix, phyloseq.object@sam_data$Fungicide, permutations = 999))
}

fungi.pcoa <- function(phyloseq.object){
  ordination <- ordinate(phyloseq.object, "PCoA", "bray")
  plot1 <- plot_ordination(phyloseq.object, ordination = ordination, type = "samples") 
  plot1.data <- plot1$data
  
  pcoa1.perc.variation.R3 <- as.character(round(ordination$values$Relative_eig[1]*100, 2))
  pcoa2.perc.variation.R3 <- as.character(round(ordination$values$Relative_eig[2]*100, 2))
  
  plot2 <- ggplot() + 
    geom_point(data = plot1.data, aes(x = Axis.1, y = Axis.2, shape = Treatment, fill = Fungicide), alpha = 0.8, size = 2) +
    theme_bw() +
    ylab(paste("PCoA2", "-",pcoa2.perc.variation.R3, "%")) + 
    xlab(paste("PCoA1", "-",pcoa1.perc.variation.R3, "%")) +
    scale_fill_manual(values=cbbPalette) +
    scale_shape_manual(values=c(21, 22, 23)) +
    guides(fill=guide_legend(override.aes=list(shape=21)))
  
  return(plot2)
}

### Functions from https://github.com/hzi-bifo/OligoMM/blob/master/constrained.pcao.functions.R
variability_table <- function(cca){
  
  chi <- c(cca$tot.chi,
           cca$CCA$tot.chi, cca$CA$tot.chi)
  variability_table <- cbind(chi, chi/chi[1])
  colnames(variability_table) <- c("inertia", "proportion")
  rownames(variability_table) <- c("total", "constrained", "unconstrained")
  return(variability_table)
  
}
cca_ci <- function(cca, permutations=5000){
  
  var_tbl <- variability_table(cca)
  p <- permutest(cca, permutations=permutations)
  ci <- quantile(p$F.perm, c(.05,.95))*p$chi[1]/var_tbl["total", "inertia"]
  return(ci)
  
}

library(exactRankTests)
library(nlme)
library(dplyr)

# OTU table should be a matrix/data.frame with each feature in rows and sample in columns. 
# Metadata should be a matrix/data.frame containing the sample identifier. 

# Data Pre-Processing
feature_table_pre_process = function(feature_table, meta_data, sample_var, group_var = NULL, 
                                     out_cut = 0.05, zero_cut = 0.90, lib_cut = 1000, neg_lb){
  feature_table = data.frame(feature_table, check.names = FALSE)
  meta_data = data.frame(meta_data, check.names = FALSE)
  # Drop unused levels
  meta_data[] = lapply(meta_data, function(x) if(is.factor(x)) factor(x) else x)
  # Match sample IDs between metadata and feature table
  sample_ID = intersect(meta_data[, sample_var], colnames(feature_table))
  feature_table = feature_table[, sample_ID]
  meta_data = meta_data[match(sample_ID, meta_data[, sample_var]), ]
  
  # 1. Identify outliers within each taxon
  if (!is.null(group_var)) {
    group = meta_data[, group_var]
    z = feature_table + 1 # Add pseudo-count (1) 
    f = log(z); f[f == 0] = NA; f = colMeans(f, na.rm = T)
    f_fit = lm(f ~ group)
    e = residuals(f_fit)
    y = t(t(z) - e)
    
    outlier_check = function(x){
      # Fitting the mixture model using the algorithm of Peddada, S. Das, and JT Gene Hwang (2002)
      mu1 = quantile(x, 0.25, na.rm = T); mu2 = quantile(x, 0.75, na.rm = T)
      sigma1 = quantile(x, 0.75, na.rm = T) - quantile(x, 0.25, na.rm = T); sigma2 = sigma1
      pi = 0.75
      n = length(x)
      epsilon = 100
      tol = 1e-5
      score = pi*dnorm(x, mean = mu1, sd = sigma1)/((1 - pi)*dnorm(x, mean = mu2, sd = sigma2))
      while (epsilon > tol) {
        grp1_ind = (score >= 1)
        mu1_new = mean(x[grp1_ind]); mu2_new = mean(x[!grp1_ind])
        sigma1_new = sd(x[grp1_ind]); if(is.na(sigma1_new)) sigma1_new = 0
        sigma2_new = sd(x[!grp1_ind]); if(is.na(sigma2_new)) sigma2_new = 0
        pi_new = sum(grp1_ind)/n
        
        para = c(mu1_new, mu2_new, sigma1_new, sigma2_new, pi_new)
        if(any(is.na(para))) break
        
        score = pi_new * dnorm(x, mean = mu1_new, sd = sigma1_new)/
          ((1-pi_new) * dnorm(x, mean = mu2_new, sd = sigma2_new))
        
        epsilon = sqrt((mu1 - mu1_new)^2 + (mu2 - mu2_new)^2 + 
                         (sigma1 - sigma1_new)^2 + (sigma2 - sigma2_new)^2 + (pi - pi_new)^2)
        mu1 = mu1_new; mu2 = mu2_new; sigma1 = sigma1_new; sigma2 = sigma2_new; pi = pi_new
      }
      
      if(mu1 + 1.96 * sigma1 < mu2 - 1.96 * sigma2){
        if(pi < out_cut){
          out_ind = grp1_ind
        }else if(pi > 1 - out_cut){
          out_ind = (!grp1_ind)
        }else{
          out_ind = rep(FALSE, n)
        }
      }else{
        out_ind = rep(FALSE, n)
      }
      return(out_ind)
    }
    out_ind = t(apply(y, 1, function(i) unlist(tapply(i, group, function(j) outlier_check(j)))))
    feature_table[out_ind] = NA
  }
  
  # 2. Discard taxa with zeros  >=  zero_cut
  zero_prop = apply(feature_table, 1, function(x) sum(x == 0, na.rm = T)/length(x[!is.na(x)]))
  taxa_del = which(zero_prop >= zero_cut)
  if(length(taxa_del) > 0){
    feature_table = feature_table[- taxa_del, ]
  }
  
  # 3. Discard samples with library size < lib_cut
  lib_size = colSums(feature_table, na.rm = T)
  if(any(lib_size < lib_cut)){
    subj_del = which(lib_size < lib_cut)
    feature_table = feature_table[, - subj_del]
    meta_data = meta_data[- subj_del, ]
  }
  
  # 4. Identify taxa with structure zeros
  if (!is.null(group_var)) {
    group = factor(meta_data[, group_var])
    present_table = as.matrix(feature_table)
    present_table[is.na(present_table)] = 0
    present_table[present_table != 0] = 1
    
    p_hat = t(apply(present_table, 1, function(x)
      unlist(tapply(x, group, function(y) mean(y, na.rm = T)))))
    samp_size = t(apply(feature_table, 1, function(x)
      unlist(tapply(x, group, function(y) length(y[!is.na(y)])))))
    p_hat_lo = p_hat - 1.96 * sqrt(p_hat * (1 - p_hat)/samp_size)
    
    struc_zero = (p_hat == 0) * 1
    # Whether we need to classify a taxon into structural zero by its negative lower bound?
    if(neg_lb) struc_zero[p_hat_lo <= 0] = 1
    
    # Entries considered to be structural zeros are set to be 0s
    struc_ind = struc_zero[, group]
    feature_table = feature_table * (1 - struc_ind)
    
    colnames(struc_zero) = paste0("structural_zero (", colnames(struc_zero), ")")
  }else{
    struc_zero = NULL
  }
  
  # 5. Return results
  res = list(feature_table = feature_table, meta_data = meta_data, structure_zeros = struc_zero, outliers = out_ind, zero_proportions = taxa_del)
  return(res)
}

# ANCOM main function
ANCOM2 = function(feature_table, meta_data, struc_zero = NULL, main_var, p_adj_method = "BH", 
                 alpha = 0.05, adj_formula = NULL, rand_formula = NULL){
  # OTU table transformation: 
  # (1) Discard taxa with structural zeros (if any); (2) Add pseudocount (1) and take logarithm.
  if (!is.null(struc_zero)) {
    num_struc_zero = apply(struc_zero, 1, sum)
    comp_table = feature_table[num_struc_zero == 0, ]
  }else{
    comp_table = feature_table
  }
  comp_table = log(as.matrix(comp_table) + 1)
  n_taxa = dim(comp_table)[1]
  taxa_id = rownames(comp_table)
  n_samp = dim(comp_table)[2]
  
  # Determine the type of statistical test and its formula.
  if (is.null(rand_formula) & is.null(adj_formula)) {
    # Basic model
    # Whether the main variable of interest has two levels or more?
    if (length(unique(meta_data%>%pull(main_var))) == 2) {
      # Two levels: Wilcoxon rank-sum test
      tfun = exactRankTests::wilcox.exact
    } else{
      # More than two levels: Kruskal-Wallis test
      tfun = stats::kruskal.test
    }
    # Formula
    tformula = formula(paste("x ~", main_var, sep = " "))
  }else if (is.null(rand_formula) & !is.null(adj_formula)) {
    # Model: ANOVA
    tfun = stats::aov
    # Formula
    tformula = formula(paste("x ~", main_var, "+", adj_formula, sep = " "))
  }else if (!is.null(rand_formula)) {
    # Model: Mixed-effects model
    tfun = nlme::lme
    # Formula
    if (is.null(adj_formula)) {
      # Random intercept model
      tformula = formula(paste("x ~", main_var))
    }else {
      # Random coefficients/slope model
      tformula = formula(paste("x ~", main_var, "+", adj_formula))
    }
  }
  
  # Calculate the p-value for each pairwise comparison of taxa.
  p_data = matrix(NA, nrow = n_taxa, ncol = n_taxa)
  colnames(p_data) = taxa_id
  rownames(p_data) = taxa_id
  for (i in 1:(n_taxa - 1)) {
    # Loop through each taxon.
    # For each taxon i, additive log ratio (alr) transform the OTU table using taxon i as the reference.
    # e.g. the first alr matrix will be the log abundance data (comp_table) recursively subtracted 
    # by the log abundance of 1st taxon (1st column) column-wisely, and remove the first i columns since:
    # the first (i - 1) columns were calculated by previous iterations, and
    # the i^th column contains all zeros.
    alr_data = apply(comp_table, 1, function(x) x - comp_table[i, ]) 
    # apply(...) allows crossing the data in a number of ways and avoid explicit use of loop constructs.
    # Here, we basically want to iteratively subtract each column of the comp_table by its i^th column.
    alr_data = alr_data[, - (1:i), drop = FALSE]
    n_lr = dim(alr_data)[2] # number of log-ratios (lr)
    alr_data = cbind(alr_data, meta_data) # merge with the metadata
    
    # P-values
    if (is.null(rand_formula) & is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        tfun(tformula, data = data.frame(x, alr_data, check.names = FALSE))$p.value
      }
      ) 
    }else if (is.null(rand_formula) & !is.null(adj_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE), 
                   na.action = na.omit)
        summary(fit)[[1]][main_var, "Pr(>F)"]
      }
      )
    }else if (!is.null(rand_formula)) {
      p_data[-(1:i), i] = apply(alr_data[, 1:n_lr, drop = FALSE], 2, function(x){
        fit = tfun(fixed = tformula, 
                   data = data.frame(x, alr_data, check.names = FALSE),
                   random = formula(rand_formula),
                   na.action = na.omit)
        anova(fit)[main_var, "p-value"]
      }
      ) 
    }
  }
  # Complete the p-value matrix.
  # What we got from above iterations is a lower triangle matrix of p-values.
  p_data[upper.tri(p_data)] = t(p_data)[upper.tri(p_data)]
  diag(p_data) = 1 # let p-values on diagonal equal to 1
  
  # Multiple comparisons correction.
  q_data = apply(p_data, 2, function(x) p.adjust(x, method = p_adj_method))
  
  # Calculate the W statistic of ANCOM.
  # For each taxon, count the number of q-values < alpha.
  W = apply(q_data, 2, function(x) sum(x < alpha))
  
  # Organize outputs
  out = data.frame(taxa_id, W, row.names = NULL, check.names = FALSE)
  # Declare a taxon to be differentially abundant based on the quantile of W statistic.
  # We perform (n_taxa - 1) hypothesis testings on each taxon, so the maximum number of rejections is (n_taxa - 1).
  out = out%>%mutate(detected_0.9 = ifelse(W > 0.9 * (n_taxa -1), TRUE, FALSE),
                     detected_0.8 = ifelse(W > 0.8 * (n_taxa -1), TRUE, FALSE),
                     detected_0.7 = ifelse(W > 0.7 * (n_taxa -1), TRUE, FALSE),
                     detected_0.6 = ifelse(W > 0.6 * (n_taxa -1), TRUE, FALSE))
  
  # Taxa with structural zeros are automatically declared to be differentially abundant
  if (!is.null(struc_zero)){
    res = data.frame(taxa_id = rownames(struc_zero), W = Inf, detected_0.9 = TRUE, 
                     detected_0.8 = TRUE, detected_0.7 = TRUE, detected_0.6 = TRUE, 
                     row.names = NULL, check.names = FALSE)
    res[match(taxa_id, res$taxa_id), ] = out
  }else{
    res = out
  }
  
  return(res)
}

'%!in%' <- function(x,y)!('%in%'(x,y))

# REFORMAT TAXONOMIES --------------------------------------------------------------------------------------------------
# >>> EXTRACT LAST TAXONOMIC LEVEL ------------------------------------------------------------------
# thanks to https://rdrr.io/github/jerryzhujian9/ezR/src/R/basic.R

blank2na = function(x, na.strings=c('','.','NA','na','N/A','n/a','NaN','nan')) {
  if (is.factor(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    # the levels will be reset here
    x = factor(x)
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else if (is.character(x)) {
    lab = attr(x, 'label', exact = T)
    labs1 <- attr(x, 'labels', exact = T)
    labs2 <- attr(x, 'value.labels', exact = T)
    # trimws will convert factor to character
    x = trimws(x,'both')
    if (! is.null(lab)) lab = trimws(lab,'both')
    if (! is.null(labs1)) labs1 = trimws(labs1,'both')
    if (! is.null(labs2)) labs2 = trimws(labs2,'both')
    if (!is.null(na.strings)) {
      # convert to NA
      x[x %in% na.strings] = NA
      # also remember to remove na.strings from value labels 
      labs1 = labs1[! labs1 %in% na.strings]
      labs2 = labs2[! labs2 %in% na.strings]
    }
    if (! is.null(lab)) attr(x, 'label') <- lab
    if (! is.null(labs1)) attr(x, 'labels') <- labs1
    if (! is.null(labs2)) attr(x, 'value.labels') <- labs2
  } else {
    x = x
  }
  return(x)
}

# In the tax_table add a column naming the highest resolution taxonomy 
# achieved for each OTU, remove _ and add sp. to genera
ReformatTaxonomy <- function(dataframe){
  taxa_table <- as.data.frame(as.matrix(tax_table(dataframe)))
  # remember to do run this function only once on your dataframe
  taxa_table$Species <- as.character(taxa_table$Species)
  taxa_table[taxa_table=="Unclassified"]<- NA
  taxa_table[taxa_table=="Unidentified"]<- NA
  taxa_table[taxa_table==""]<- NA
  #taxa_table[which(is.na(taxa_table$Species) == FALSE), ]$Species <- paste(
   # taxa_table$Species[is.na(taxa_table$Species)==FALSE], "sp.", sep = " ")
  taxa_table <- taxa_table[c(8,1,2,3,4,5,6,7,9,10)]
  taxa_table[] = lapply(taxa_table, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
  lastValue <- function(x) tail(x[!is.na(x)], 1)
  last_taxons<- apply(taxa_table[,1:8], 1, lastValue)
  taxa_table$BestMatch <- last_taxons
  taxa_table[, "BestMatch"] <- gsub("_", " ", taxa_table[, "BestMatch"])
  taxa_table$BestMatch <- ifelse(is.na(taxa_table$Species == TRUE), paste(taxa_table$BestMatch, "sp.", sep = " "), taxa_table$Species)
  taxa_table$Taxonomy <- paste(taxa_table$OTU_ID, taxa_table$BestMatch, sep="-")
  #taxa_table[, "Species"] <- gsub(" sp.", "", taxa_table[, "Species"])
  tax_table(dataframe) <- tax_table(as.matrix(taxa_table))
  return(dataframe)
}

# # In the tax_table add a column naming the highest resolution taxonomy 
# # achieved for each OTU, remove _ and add sp. to genera
# ReformatTaxonomy <- function(dataframe, taxa){
#   taxa_table <- as.data.frame(as.matrix(tax_table(dataframe)))
#   # remember to do run this function only once on your dataframe
#   taxa_table$Genus <- as.character(taxa_table$Genus)
#   taxa_table[taxa_table=="Unclassified"]<- NA
#   taxa_table[taxa_table==""]<- NA
#   taxa_table[which(is.na(taxa_table$Genus) == FALSE), ]$Genus <- paste(
#     taxa_table$Genus[is.na(taxa_table$Genus)==FALSE], "sp.", sep = " ")
#   taxa_table$OTU <- row.names(taxa_table)
#   if (taxa == "ITS"){
#     taxa_table <- taxa_table[c(8,1,2,3,4,5,6,7)]
#   }else{
#     taxa_table <- taxa_table[c(7,1,2,3,4,5,6)]
#   }
#   taxa_table[] = lapply(taxa_table, blank2na, na.strings=c('','NA','na','N/A','n/a','NaN','nan'))
#   lastValue <- function(x) tail(x[!is.na(x)], 1)
#   last_taxons<- apply(taxa_table, 1, lastValue)
#   taxa_table$BestMatch <- last_taxons
#   taxa_table[, "BestMatch"] <- gsub("_", " ", taxa_table[, "BestMatch"])
#   taxa_table$Taxonomy <- paste(taxa_table$OTU, taxa_table$BestMatch, sep="-")
#   taxa_table[, "Genus"] <- gsub(" sp.", "", taxa_table[, "Genus"])
#   tax_table(dataframe) <- tax_table(as.matrix(taxa_table))
#   return(dataframe)
# }



