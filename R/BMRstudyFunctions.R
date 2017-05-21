#BMRstudyFunctions.R

## Function for data cleaning that uses OpenTree of Life taxonomic name resolution. 
ttol_queryOttIds <- function(ttol, nrequest=500, tree=TRUE){
  if(tree){
    ttol_names <- ttol$tip.label
  } else {
    ttol_names <- ttol
  }
  ttol_names <- gsub("_", " ", ttol_names) %>% gsub("'", "", .)
  ttol_names <- lapply(as.character(ttol_names), function(x) strsplit(x, " ", fixed=TRUE)[[1]][1:2]) %>% sapply(., function(x) paste(x, collapse=" "))
  seq1 <- c(seq(1, length(ttol_names), nrequest), length(ttol_names)+1)
  cat("Searching for exact matches....\n")
  tol_otl <- lapply(1:(length(seq1)-1), function(x) cbind(seq1[x]:(seq1[x+1]-1),tnrs_match_names(ttol_names[seq1[x]:(seq1[x+1]-1)], do_approximate_matching = FALSE)))
  tol_otl <- do.call(rbind, tol_otl)
  colnames(tol_otl) <- c("No", "search_string", "unique_name", "approximate_match", "ott_id", "number_matches", "is_synonym", "is_deprecated")
  not_found <- which(is.na(tol_otl[,2]))
  cat(paste(length(ttol_names) - length(not_found), " exact matches found in OTL database\n", sep=""))
  if(length(not_found) > 0){
    cat(paste("Performing approximate search for ", length(not_found), " unmatched taxa\n", sep=""))
    seq2 <- c(seq(1, length(not_found), nrequest), length(not_found)+1)
    missing <- lapply(1:(length(seq2)-1), function(x) cbind(not_found[seq2[x]:(seq2[x+1]-1)],tnrs_match_names(ttol_names[not_found[seq2[x]:(seq2[x+1]-1)]], do_approximate_matching = TRUE)))
    missing <- do.call(rbind, missing)
    tol_otl[not_found[which(!is.na(missing[,2]))], ] <- missing[which(!is.na(missing[,2])),]
    cat(paste("Found ", length(which(!is.na(missing[,2]))), " approximate matches\n", sep=""))
    cant_find <- not_found[which(is.na(missing[,2]))]
    if(length(cant_find>0) & tree){
      new.ttol <- drop.tip(ttol, cant_find)   
      tol_otl <- tol_otl[-cant_find,]
    } else {
      new.ttol <- gsub(" ", "_", tol_otl$unique_name)
    }
    dropped=ttol_names[cant_find]
  } else {
    dropped=NULL
    new.ttol <- gsub(" ", "_", tol_otl$unique_name)
  }
  return(list(ttol=new.ttol, tol_otl=tol_otl, dropped=dropped))
}


#' Combine mcmc chains
#' 
#' @param chain.list The first chain to be combined
#' @param thin A number or vector specifying the thinning interval to be used. If a single value,
#' then the same proportion will be applied to all chains.
#' @param burnin.prop A number or vector giving the proportion of burnin from each chain to be 
#' discarded. If a single value, then the same proportion will be applied to all chains.
#' 
#' @return A combined bayouMCMC chain
#' 
#' @export
combine.chains <- function(chain.list, thin=1, burnin.prop=0){
  nns <- lapply(chain.list, function(x) names(x))
  if(length(burnin.prop) == 1){
    burnins <- rep(burnin.prop, length(chain.list))
  }
  if(length(thin) == 1){
    thins <- rep(thin, length(chain.list))
  }
  Ls <- sapply(chain.list, function(x) length(x$gen))
  if(!all(sapply(nns, function(x) setequal(nns[[1]], x)))){
    stop ("Not all chains have the same named elements and cannot be combined")
  } else {
    nn <- nns[[1]]
  }
  for(i in 1:length(chain.list)) chain.list[[i]]$gen <- chain.list[[i]]$gen + 0.1*i
  postburns <- lapply(1:length(chain.list), function(x) seq(max(c(floor(burnins[x]*Ls[x]),1)), Ls[x], thins[x]))
  chains <- setNames(vector("list", length(nns[[1]])), nns[[1]])
  attributes(chains) <- attributes(chain.list[[1]])
  for(i in 1:length(nn)){
    chains[[nn[i]]] <- do.call(c, lapply(1:length(chain.list), function(x) chain.list[[x]][[nn[i]]][postburns[[x]]]))
  }
  attributes(chains)$burnin <- 0
  return(chains)
}


#' A function for summarizing the state of a model after a shift
#' 
#' @param chain A bayouMCMC chain
#' @param mcmc A bayou mcmc object
#' @param pp.cutoff The threshold posterior probability for shifts to summarize, if 'branches' 
#' specified than this is ignored.
#' @param branches The specific branches with shifts to summarize, assuming postordered tree
#' 
#' @details shiftSummaries summarizes the immediate parameter values after a shift on a particular
#' branch. Parameters are summarized only for the duration that the particular shift exists. Thus,
#' even global parameters will be different for particular shifts. 
#' 
#' @return A list with elements: 
#' \code{pars} = a bayoupars list giving the location of shifts specified;
#' \code{tree} = The tree; 
#' \code{pred} = Predictor variable matrix; 
#' \code{dat} = A vector of the data; 
#' \code{SE} = A vector of standard errors;
#' \code{PP} = Posterior probabilities of the specified shifts; 
#' \code{model} = A list specifying the model used; 
#' \code{variables} = The variables summarized; 
#' \code{cladesummaries} = A list providing the medians and densities of the distributions of regression 
#' variables for each shift; 
#' \code{descendents} = A list providing the taxa that belong to each regime 
#' \code{regressions} = A matrix providing the regression coefficients for each regime.
#' @export
shiftSummaries <- function(chain, mcmc, pp.cutoff=0.3, branches=NULL, ...){
  cache <- .prepare.ou.univariate(mcmc$tree,mcmc$dat, SE=mcmc$SE, pred=mcmc$pred)
  tree <- cache$phy
  dat <- cache$dat
  pred <- cache$pred
  SE <- cache$SE
  model <- mcmc$model.pars
  
  if(is.null(attributes(chain)$burnin)){
    L <- Lposterior(chain, tree, burnin=0)
  } else {
    L <- Lposterior(chain, tree, burnin=attributes(chain)$burnin)
  }
  if(is.null(branches)){
    branches <- which(L[,1] > pp.cutoff)
    PP <- L[branches,1]
  } else {
    PP <- L[branches,1]
  }
  if(length(branches) == 0){
    stop("No shifts found with posterior probability above cutoff")
  }
  
  if(!is.null(model$call)){
    coefs <- paste("beta_", attr(terms(model$call), "term.labels"), sep="")
  } else {
    coefs <- NULL
  }
  variables <- c('theta', coefs)
  sumpars <- list(k=length(branches), ntheta=length(branches)+1, sb=branches, t2=2:(length(branches)+1), loc=rep(0, length(branches)))
  .summarizeDerivedState <- function(branch, chain, variables){
    if(branch==0){
      values <- lapply(variables, function(x) sapply(chain[[x]], function(y) y[1]))
    } else {
      SB <- unlist(chain$sb)
      gen <- unlist(lapply(1:length(chain$sb), function(x) rep(x, length(chain$sb[[x]]))))
      ind <- which(SB==branch)
      gen <- gen[ind]
      T2 <- unlist(chain$t2)[ind]
      values <- lapply(variables, function(x) sapply(1:length(T2), function(y)if(length(chain[[x]][[gen[y]]]) > 1) {chain[[x]][[gen[y]]][T2[y]]} else chain[[x]][[gen[y]]]))
    }
    medians <- lapply(values, median)
    densities <- lapply(values, density)
    names(medians) <- names(densities) <- variables
    return(list(medians=medians, densities=densities))
  }
  cladesummaries <- lapply(c(0, branches), function(x) .summarizeDerivedState(x, chain, variables))
  regressions <- do.call(rbind, lapply(cladesummaries, function(x) unlist(x$medians)))
  rownames(regressions) = c("root", sumpars$sb)
  tipregs <- .tipregime(sumpars, tree)
  descendents <- lapply(1:(length(sumpars$sb)+1), function(x) names(tipregs[tipregs==x])) 
  out <- list(pars= sumpars, tree=tree, pred=pred, dat=dat, SE=SE, PP=PP, model=model, variables=variables, cladesummaries=cladesummaries, descendents=descendents, regressions=regressions)
  return(out)
}

#' A function to plot a heatmap of reconstructed parameter values on the branches of the tree
#' 
#' @param tree A phylogenetic tree
#' @param chain A bayou MCMC chain
#' @param variable The parameter to reconstruct across the tree
#' @param burnin The initial proportion of burnin samples to discard 
#' @param nn The number of discrete categories to divide the variable into
#' @param pal A color palette function that produces nn colors
#' @param legend_ticks The sequence of values to display a legend for
#' @param legend_settings A list of legend attributes (passed to addColorBar)
#' @param ... Additional options passed to plot.phylo
#' 
#' @details legend_settings is an optional list of any of the following:
#' 
#' legend - a logical indicating whether a legend should be plotted
#' 
#' x - the x location of the legend
#' 
#' y - the y location of the legend
#' 
#' height - the height of the legend
#' 
#' width - the width of the legend
#' 
#' n - the number of gradations in color to plot from the palette
#' 
#' adjx - an x adjustment for placing text next to the legend bar
#' 
#' cex.lab - the size of text labels next to the legend bar
#' 
#' text.col - The color of text labels
#' 
#' locator - if TRUE, then x and y coordinates are ignored and legend is placed
#' interactively.
#' 
#' #' @export
plotBranchHeatMap <- function(tree, chain, variable, burnin=0, nn=NULL, pal=heat.colors, legend_ticks=NULL, legend_settings=list(plot=TRUE), ...){
  dum <- setNames(rep(1, length(tree$tip.label)), tree$tip.label)
  cache <- .prepare.ou.univariate(tree, dum)
  tree <- cache$phy
  seq1 <- floor(max(seq(burnin*length(chain$gen),1), length(chain$gen), 1))
  if(is.null(legend_ticks)){
    legend_ticks <- seq(min(unlist(chain[[variable]][seq1],F,F)), max(unlist(chain[[variable]][seq1],F,F)), length.out=5)
  }
  if(is.null(nn)) nn <- length(seq1) else { seq1 <- floor(seq(max(burnin*length(chain$gen),1), length(chain$gen), length.out=nn))}
  if(length(nn) > length(chain$gen)) stop("Number of samples greater than chain length, lower nn")
  abranches <- lapply(1:nrow(tree$edge), .ancestorBranches, cache=cache)
  allbranches <- suppressWarnings(sapply(1:nrow(tree$edge), function(x) .branchRegime(x, abranches, chain, variable, seq1, summary=TRUE)))
  plot(tree, edge.color=colorRamp(allbranches, pal, 100), ...)
  lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv)
  legend_stuff <- list(x=0.01* lastPP$x.lim[2], 
                       y=0, 
                       height=0.25*diff(lastPP$y.lim), 
                       width=0.01*diff(lastPP$x.lim), 
                       n=100, 
                       trait=allbranches, 
                       ticks=legend_ticks, 
                       adjx=0.01*lastPP$x.lim[2], 
                       cex.lab=0.5, 
                       text.col="black",
                       plot=TRUE,
                       locator=FALSE
  )
  if(length(legend_settings) > 0){
    for(i in 1:length(legend_settings)){
      legend_stuff[[names(legend_settings)[i]]] <- legend_settings[[i]]
    }
  }
  if(legend_stuff$plot) {
    if(legend_stuff$locator){
      lc <- locator(1)
      legend_stuff$x <- lc$x
      legend_stuff$y <- lc$y
      addColorBar(x=legend_stuff$x, y=legend_stuff$y, height=legend_stuff$height, width=legend_stuff$width, pal=pal, n=legend_stuff$n, trait=allbranches, ticks=legend_ticks, adjx=legend_stuff$adjx, cex.lab=legend_stuff$cex.lab, text.col=legend_stuff$text.col)
    } else addColorBar(x=legend_stuff$x, y=legend_stuff$y, height=legend_stuff$height, width=legend_stuff$width, pal=pal, n=legend_stuff$n, trait=allbranches, ticks=legend_ticks, adjx=legend_stuff$adjx, cex.lab=legend_stuff$cex.lab, text.col=legend_stuff$text.col)
  }
}


.ancestorBranches <- function(branch, cache){
  ancbranches <- which(sapply(cache$bdesc, function(x) branch %in% x))
  sort(ancbranches, decreasing=FALSE)
}
.branchRegime <- function(branch, abranches, chain, parameter, seqx, summary=FALSE){
  ancs <- c(branch, abranches[[branch]])
  ancshifts <- lapply(1:length(seqx), function(x) chain$t2[[seqx[x]]][which(chain$sb[[seqx[x]]] == ancs[min(which(ancs %in% chain$sb[[seqx[x]]]))])])
  ancshifts <- sapply(ancshifts, function(x) ifelse(length(x)==0, 1, x))
  ests <- sapply(1:length(ancshifts), function(x) chain[[parameter]][[seqx[x]]][ancshifts[x]])
  res <- cbind(ests)
  if(summary){
    return(apply(res, 2, stats::median))
  } else {
    return(res)
  }
}

