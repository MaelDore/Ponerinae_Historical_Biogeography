#######################################################
# BSM_to_phytools_SM
#######################################################
#' Convert Biogeographic Stochastic Map (BSM) to phytools SIMMAP stochastic map (SM) format
#' 
#' This function converts a Biogeographic Stochastic Map (BSM) to a \code{phytools} "\code{simmap}" 
#' object. This is useful mostly to make use of the more diverse plotting functions for
#' simmaps available in \code{phytools}. 
#'
#' The BioGeoBEARS Biogeographic Stochastic Map (BSM) output is rather complex, as it
#' keeps track of numerous different types of anagenetic and cladogenetic events, the
#' source and destination areas of dispersal events (*), etc. However, BioGeoBEARS only 
#' has plotting functions for the standard, rectangular, right-facing phylogeny. 
#' \code{\link{BSM_to_phytools_SM}} converts a BioGeoBEARS BSM result into a phytools 
#' stochastic mapping object (which is a really a phylo3 tree object, with extra fields
#' to contain the stochastic map.
#'
#' The inputs to \code{BSM_to_phytools_SM} consist of (1) a "res" object (the result of 
#' an ML inference), a \code{clado_events_table} containing
#' the events sampled at each node of the phylogeny, and (optionally, but usually needed)
#' a \code{ana_events_table} containing the anagenetic events on each branch (or branch 
#' segment, in the case of a time-stratified analysis, where the tree has been broken
#' up into subpieces.
#'
#' The resulting object, \code{tr_wSimmap}, is a phytools simmap object, which really 
#' consists of a standard APE phylogeny object, plus the extra fields "map", "mapped.edge",
#' "Q", and "logL". The class of this object, i.e. the result of class(tr_wSimmap), is
#' \code{c("simmap", "phylo")}. 
#'
#' Notes on using the resulting simmap object in phytools:
#'
#' <b>1.</b> The phytools functions, like \code{countSimmap(tr_wSimmap)}, will only count the 
#'    \emph{anagenetic} events, as that is all that is formally recorded in the \code{phytools}
#'    \code{simmap} object. For example, imagine if your model parameters for a DEC+J model were
#'    \emph{d}=0.00224, \emph{e}=0.0, \emph{j}=0.0297. In this case, \code{countSimmap(tr_wSimmap)}
#'    would count transitions like \code{A->AB}, because d is positive. However, the 
#'    reverse transition of \code{AB->A} would require an "\emph{e}" event, but \emph{e} is 
#'    inferred to be 0.0, as is very common in DEC, DEC+J etc. analyses. 
#'
#'    Transitions to 
#'    single-area states in "\emph{j}" events would be common at cladogenesis events (e.g., 
#'    \code{AB->AB,C}), but \code{phytools} doesn't know anything about cladogenesis models, 
#'    as it was written assuming purely anagenetic models.  One could probably write a 
#'    script to infer "<i>j</i>" events off the \code{phytools} \code{simmap} object by 
#'    comparing the last-state-below-a-node with the descendent-pairs-above-a-node, but 
#'    it would probably be easier to just use the BioGeoBEARS BSM outputs 
#'    (\code{clado_events_table} and \code{ana_events_table}), because they record all 
#'    of this information already, and the text files output by the example BSM script
#'    summarize all of the different types of events and their direction.
#'
#' \bold{2.} For similar reasons, the \code{phytools} graphics, while they should branch
#' the branch histories correctly, don't always correctly 
#' draw the colors at the "corners" between a speciation event and the beginning of a 
#' descendant branch. This is for the same reasons as #1: \code{phytools} doesn't 
#' know about cladogenesis models and assumes a purely anagenetic model, where the state
#' just before and just after cladogenesis would always be identical. (Counting the 
#' states at nodes may also be somewhat inaccurate if counted off the \code{phytools} \code{simmap}
#' derived from a BioGeoBEARS BSM, I haven't tested this.)
#' 
#' \bold{(*) Please note carefully:} area-to-area dispersal events are not identical with the 
#' state transitions. For example, a state can be a geographic range with multiple 
#' areas, but under the logic of DEC-type models, a range-expansion event like 
#' ABC->ABCD actually means that a dispersal happened from some specific area (A, B, or C)
#' to the new area, D.  BSMs track this area-to-area sourcing, at least if 
#' \code{\link{simulate_source_areas_ana_clado}} has been run.
#' 
#' @param \code{res} A BioGeoBEARS results object, produced by ML inference via \code{\link{bears_optim_run}}.
#' @param \code{clado_events_table} A cladogenetic events table, from BioGeoBEARS BSM.
#' @param \code{ana_events_table} An anagenetic events table, from BioGeoBEARS BSM. Default is NA, 
#' in which case the function assumes that the BSM had 0 anagenetic range-changing events, and 
#' all range-changing events were cladogenetic. This will only process properly if actually is
#' true that the sampled states at branch bottoms and branch tops in the the \code{clado_events_table}
#' always match (an error check checks for this).
#' @return \code{tr_wSimmap} This is an object of type \code{c("simmap", "phylo"). It contains
#' the following: 
#' \code{tr_wSimmap$tr} is a \code{phylo3} APE tree object
#' \code{tr_wSimmap$maps} = A list for each branch (branch numbers when \code{phylo} object in "cladewise" order),
#' with the state history of each branch (the length of time in each state list of states inhabited, with the 
#' cell names giving the state).
#' \code{tr_wSimmap$mapped.edge} The total amount of time in each state, for each branch (same order as
#' in \code{$maps}).  The row names are the names of the nodes at the bottom and top of each branch/edge.
#' tr_wSimmap$Q The \code{Q} transition matrix (calculated from the ML parameter estimates contained in f;cc).
#' tr_wSimmap$logL The log-likelihood of the data under the ML model (this will be the same across all BSMs).
#' @export
#' @seealso \code{\link{BSMs_to_phytools_SMs}}, \code{\link{simulate_source_areas_ana_clado}}, 
#' \code{phytools::countSimmap},  \code{phytools::make.simmap}
#' @note Go (BioGeo)BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' 
#' @examples
#' 

BSM_to_phytools_SM_custom <- function(res, tree, clado_events_table, ana_events_table = NA)
{
  run_parent_brs_TF = TRUE
  tr <- tree

  # Error check
  if ((class(ana_events_table) == "data.frame") && (dim(ana_events_table)[1] > 0))
  {
    run_parent_brs_TF = TRUE
  } else {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }
  
  if (is.null(ana_events_table) == TRUE)
  {
    ana_events_table = NA
    run_parent_brs_TF = FALSE
  }
  
  # Is it time-stratified?
  stratTF = (length(res$inputs$timeperiods) > 0)
  
  returned_mats = get_Qmat_COOmat_from_BioGeoBEARS_run_object(BioGeoBEARS_run_object=res$inputs, include_null_range=res$inputs$include_null_range)
  returned_mats
  areanames = returned_mats$areanames
  ranges_list = returned_mats$ranges_list
  
  # Fix for issue with using only reduced range list
  reduced_ranges_list = returned_mats$ranges_list
  max_range_size <- max(nchar(reduced_ranges_list))
  ranges_list <- generate_list_ranges(areas_list = areanames, max_range_size = max_range_size, include_null_range = res$inputs$include_null_range)
  
  # Get list of edges:
  trtable = prt(tr, printflag=FALSE)
  trtable
  
  # View(trtable)
  # View(clado_events_table)
  
  # Convert pruningwise edge numbers to cladewise edgenums
  # The clado_events_table (non-stratified) has the edge numbers in
  # pruningwise order
  if (stratTF == FALSE)
  {
    pruningwise_edgenums = clado_events_table$parent_br
    cladewise_edgenums = trtable$parent_br
    translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
    translation_pruning_to_clade_edgenums
    
    # Error trap for when there are no anagenetic events
    if (run_parent_brs_TF == TRUE)
    {
      ana_events_edgenums_indexes_in_clado_events_table = match(x=ana_events_table$parent_br, table=clado_events_table$parent_br)
      ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
    }
    clado_events_table$parent_br = trtable$parent_br
  # } else {
  #   # If same issue with time-stratified clado_events_table, also have duplicates due to branches split across time-strata!
  #   pruningwise_edgenums = unique(clado_events_table$parent_br)
  #   cladewise_edgenums = trtable$parent_br
  #   translation_pruning_to_clade_edgenums = as.data.frame(cbind(pruningwise_edgenums, cladewise_edgenums), stringsAsFactors=FALSE)
  #   translation_pruning_to_clade_edgenums
  #   # Error trap for when there are no anagenetic events
  #   if (run_parent_brs_TF == TRUE)
  #   {
  #     # ana_events_edgenums_indexes_in_clado_events_table = match(x = ana_events_table$parent_br, table = clado_events_table$parent_br)
  #     # ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[ana_events_edgenums_indexes_in_clado_events_table]
  #     ana_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = ana_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
  #   }
  #   clado_events_table$parent_br = translation_pruning_to_clade_edgenums$cladewise_edgenums[match(x = clado_events_table$parent_br, table = translation_pruning_to_clade_edgenums$pruningwise_edgenums)]
  }
   
  
  # Get the edgenums, exclude the "NA" for the root branch
  # Order edgenums from smallest to largest
  nonroot_TF = trtable$node.type != "root"
  edgenums = trtable$parent_br[nonroot_TF]
  edgenums_order = order(edgenums)
  edgenums = edgenums[edgenums_order]
  
  # Get the "ancestor node, descendant node"
  ancnodenums = trtable$ancestor[nonroot_TF][edgenums_order]
  decnodenums = trtable$node[nonroot_TF][edgenums_order]
  
  # rownames for mapped.edge
  rownames_for_mapped_edge = paste0(ancnodenums, ",", decnodenums)
  rownames_for_mapped_edge
  
  # instantiate "maps" for phytools (a list, with array of state residence times
  maps = list()
  
  if (stratTF == TRUE)
  {
    time_tops = sort(unique(clado_events_table$time_top))
    time_bots = sort(unique(clado_events_table$time_bot))
  }
  
  
  # Loop through the edges, record any anagenetic events on the branches
  # fill in "maps"
  for (i in 1:length(edgenums))
  {
    # i <- 410
    edgenum = edgenums[i]
    
    # Trap for if ana_events_table is NA (common, if there are no anagenetic events in
    # the tree at all)
    if ( (length(ana_events_table) == 1) && (is.na(ana_events_table)) )
    {
      edgefound_TF = FALSE
    } else {
      edgefound_TF = ana_events_table$parent_br == edgenum
    }
    
    # If no anagenetic events are found, the whole branchlength is in the 
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) == 0)
    {
      # The states should be the same at the branch bottom and top:
      clado_row_TF = clado_events_table$parent_br == edgenum
      clado_row_TF[is.na(clado_row_TF)] = FALSE
      
      # 			if (stratTF == TRUE)
      # 				{
      # 				match_time_tops_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) >= as.numeric(trtable$time_top)
      # 				match_time_bots_TF = as.numeric(tmptable_rows$abs_event_time[nnr]) < as.numeric(trtable$time_bot)
      # 				}
      
      
      ## NJM 2019-03-12_ fix: doubles can be found in time-strat, FIX
      clado_events_table[clado_row_TF,]
      
      
      # Error check
      if (sum(clado_row_TF) < 1)
      {
        txt = paste0("STOP ERROR #1 in BSM_to_phytools_SM(): 0 rows in clado_events_table match edge number/branch number (parent_br==", edgenum, ").")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      
      if (sum(clado_row_TF) == 1)
      {
        bottom_state_num_1based = clado_events_table$sampled_states_AT_brbots[clado_row_TF]
        top_state_num_1based = clado_events_table$sampled_states_AT_nodes[clado_row_TF]
        
        # Error handle for the root
        if (clado_events_table$node.type[clado_row_TF] == "root")
        {
          bottom_state_num_1based <- top_state_num_1based
        }
      }
      
      if (sum(clado_row_TF) > 1)
      {
        bottom_state_num_1based = unique(clado_events_table$sampled_states_AT_brbots[clado_row_TF])
        top_state_num_1based = unique(clado_events_table$sampled_states_AT_nodes[clado_row_TF])
        
        # Error check: there should be only 1 unique state corresponding to this
        # node (because we are in the section where no anagenetic histories were
        # found on the branch).
        if (length(top_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'top_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")
          
          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])
          
          print("top_state_num_1based")
          print(top_state_num_1based)
          
          stop(txt)
        } # END if (length(top_state_num_1based) != 1)
        
        if (length(bottom_state_num_1based) != 1)
        {
          txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'bottom_state_num_1based' corresponding to the node specified. Printing the matching rows of 'clado_events_table'."
          cat("\n\n")
          cat(txt)
          cat("\n\n")
          
          print("clado_events_table[clado_row_TF,]")
          print(clado_events_table[clado_row_TF,])
          
          print("bottom_state_num_1based")
          print(bottom_state_num_1based)
          
          stop(txt)
        } # END if (length(bottom_state_num_1based) != 1)
      } # END if (sum(clado_row_TF) > 1)
      
      # Error check
      if (bottom_state_num_1based != top_state_num_1based)
      {
        txt = paste0("STOP ERROR #2 in BSM_to_phytools_SM(): the top_state_num_1based (", top_state_num_1based, "), and bottom_state_num_1based (", bottom_state_num_1based, ") have to match at edge number/branch number (parent_br==", edgenum, "), because no anagenetic events were recorded on this branch.")
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        stop(txt)
      }
      
      # No events detected, so put in the states_array just one state
      # But, since there are NO events on this branch, don't add to it if
      # there is already something there (e.g. don't add the same branch for
      # multiple branch segments)
      names_of_states_array = c(ranges_list[bottom_state_num_1based])
      times_in_each_state_array = unique(c(clado_events_table$edge.length[clado_row_TF]))
      
      if (length(times_in_each_state_array) != 1)
      {
        txt = "STOP ERROR in BSM_to_phytools_SM(): more than one 'times_in_each_state_array' corresponding to the node specified. There should only be one time, because no anagenetic events were detected on this branch. Printing the matching rows of 'clado_events_table'."
        cat("\n\n")
        cat(txt)
        cat("\n\n")
        
        print("clado_events_table[clado_row_TF,]")
        print(clado_events_table[clado_row_TF,])
        
        print("times_in_each_state_array")
        print(times_in_each_state_array)
        
        stop(txt)
      } # END if (length(times_in_each_state_array) != 1)
      
      
      names(times_in_each_state_array) = names_of_states_array
      times_in_each_state_array
    } # END if (sum(edgefound_TF) == 0)
    
    # If some anagenetic events are found, the whole branchlength is in the 
    # starting state, as specified in "clado_events_table"
    if (sum(edgefound_TF) > 0)
    {
      # The states will be listed in the ana_events_table
      rows_matching_edgenum_TF = ana_events_table$parent_br == edgenum
      tmp_ana_events_table = ana_events_table[rows_matching_edgenum_TF,]
      
      # Make sure the tmp_ana_events_table is sorted by REVERSE event_time (along branch)
      # (Do REVERSE, because older events have larger event ages)
      tmp_ana_events_table = tmp_ana_events_table[rev(order(tmp_ana_events_table$abs_event_time)),]
      tmp_ana_events_table
      
      numevents = sum(rows_matching_edgenum_TF)
      
      # 1 event, 2 states on branch
      if (numevents == 1)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])
        
        # The rest of the branch is the 2nd state
        further_state_name = c(tmp_ana_events_table$new_rangetxt[1])
        further_state_time = c(tmp_ana_events_table$edge.length[1] - first_state_timelength)
        
        times_in_each_state_array = c(first_state_timelength, further_state_time)
        names_of_states_array = c(first_state_name, further_state_name)
      } # END if (numevents == 1)
      
      # 2+ events, 3+ states on branch
      if (numevents >= 2)
      {
        first_state_name = c(tmp_ana_events_table$current_rangetxt[1])
        #first_state_time = c(tmp_ana_events_table$event_time[1])
        abs_time_bp_at_branch_top = tmp_ana_events_table$time_bp[1]
        abs_time_bp_at_branch_bot = abs_time_bp_at_branch_top + tmp_ana_events_table$edge.length[1]
        first_state_timelength = c(abs_time_bp_at_branch_bot - tmp_ana_events_table$abs_event_time[1])
        
        nonfirst_rows = 2:numevents
        nonlast_rows = 1:(numevents-1)
        
        further_state_names = c(tmp_ana_events_table$new_rangetxt)
        further_state_times = tmp_ana_events_table$abs_event_time[nonlast_rows] - tmp_ana_events_table$abs_event_time[nonfirst_rows]
        
        # How much of the branch is left
        last_time = c(tmp_ana_events_table$edge.length[numevents] - sum(c(first_state_timelength,further_state_times)))
        
        times_in_each_state_array = c(first_state_timelength, further_state_times, last_time)
        names_of_states_array = c(first_state_name, further_state_names)
      } # END if (numevents >= 2)
      
      names(times_in_each_state_array) = names_of_states_array
    } # END if (sum(edgefound_TF) > 0)
    
    # Store
    maps[[i]] = times_in_each_state_array
  } # END for (i in 1:length(edgenums))
  
  maps
  
  
  # Check that the sums of state residence times add up to the branch lengths
  simmap_times <- unlist(lapply(X=maps, FUN=sum))
  tree_times <- tr$edge.length
  all(simmap_times == tree_times)
  cbind(simmap_times, tree_times)
  
  # Check if issue with order of edges
  simmap_times <- simmap_times[order(simmap_times)]
  tree_times <- tree_times[order(tree_times)]
  cbind(simmap_times, tree_times)
  
  # Compute residence times per states
  edge_states_matrix <- plyr::ldply(.data = maps, .fun = base::rbind)
  state_times <- apply(X = edge_states_matrix, MARGIN = 2, FUN = sum, na.rm = T)
  state_times_reordered <- state_times[ranges_list]
  names(state_times_reordered) <- ranges_list
  state_times_reordered[is.na(state_times_reordered)] <- 0
  
  # Add total
  state_times_df <- as.data.frame(t(state_times_reordered))
  state_times_df$total <- sum(state_times_reordered)
  
  # Add percentage of times
  state_times_perc <- state_times_df / state_times_df$total * 100
  state_times_df <- rbind(state_times_df, state_times_perc)
  row.names(state_times_df) <- c("raw", "perc")
  
  # # Remove empty states
  # state_times_df <- state_times_df[, state_times_perc != 0]
  
  # Make the mapped.edge output
  mapped.edge = matrix(data=0.0, nrow=length(edgenums), ncol=length(ranges_list))
  row.names(mapped.edge) = rownames_for_mapped_edge
  colnames(mapped.edge) = ranges_list
  mapped.edge
  
  # For each branch, 
  # 1. Get the list of observed states
  i = 3
  observed_states = sort(unique(names(maps[[i]])))
  observed_states
  
  
  # sapply to get the sum of each
  sapply(X=observed_states, FUN=get_sum_statetime_on_branch, branch_history_map=maps[[i]])
  
  # Fill in the table for each branch
  for (i in 1:nrow(mapped.edge))
  {
    observed_states = sort(unique(names(maps[[i]])))
    total_residence_times = sapply(X=observed_states, FUN=get_sum_statetime_on_branch, branch_history_map=maps[[i]])
    names_observed_states = names(total_residence_times)
    mapped.edge[i,names_observed_states] = unname(total_residence_times)
  }
  mapped.edge
  
  
  # Get the transition matrix and logL
  Q = returned_mats$Qmat
  # row.names(Q) = ranges_list
  # colnames(Q) = ranges_list
  row.names(Q) = reduced_ranges_list
  colnames(Q) = reduced_ranges_list
  logL = res$total_loglikelihood
  
  tr_wSimmap = tr
  tr_wSimmap$maps = maps
  tr_wSimmap$mapped.edge = mapped.edge
  tr_wSimmap$Q = Q
  tr_wSimmap$logL = logL
  class(tr_wSimmap) = c("simmap", "phylo")
  tr_wSimmap
  
  tr_wSimmap$maps
  tr_wSimmap$mapped.edge
  
  return(list(simmap = tr_wSimmap, residence_times = state_times_df))
} 



### Same but for list of multiple BSM maps

BSMs_to_phytools_SMs_custom <- function(res, tree, clado_events_tables, ana_events_tables)
{
  simmaps_list = list()
  
  if (class(clado_events_tables) != "list")
  {
    txt = "WARNING from BSMs_to_phytools_SMs(): Input 'clado_events_tables' was not a list, so we will assume it is instead a single Biogeographical Stochastic Map (BSM) table in data.frame format. Therefore, BSM_to_phytools_SM() is being run instead."
    warning(txt)
    tr_wSimmap = BSM_to_phytools_SM_custom(res = res, tree = tree, clado_events_table = clado_events_tables, ana_events_table = ana_events_tables)$simmap
    return(tr_wSimmap)
  } # END if (class(clado_events_tables) != "list")
  
  for (i in 1:length(clado_events_tables))
  {
    simmaps_list[[i]] = BSM_to_phytools_SM_custom(res = res, tree = tree, clado_events_table = clado_events_tables[[i]], ana_events_table = ana_events_tables[[i]])$simmap
  }
  
  class(simmaps_list) = c("multiSimmap", "multiPhylo")
  return(simmaps_list)
} 
