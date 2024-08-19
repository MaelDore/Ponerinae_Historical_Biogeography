
# Original function: BioGeoBEARS::section_the_tree by Nicholas J. Matzke
# All internal utility functions are from BioGeoBEARS package by Nicholas J. Matzke

# New addition: use a ape phylogenetic tree and tipward and rootward time boundaries as inputs instead of a BioGEOBEARS object


section_the_tree_for_time_window <- function (tree, # ape phylo object
                                              tipward_time, # Tipward time boundary
                                              rootward_time, # Rootward time boundary
                                              make_master_table = FALSE, plot_pieces = TRUE, 
                                              cut_fossils = FALSE, fossils_older_than = 0.001, min_dist_between_node_and_stratum_line = 1e-05, 
                                              remove_root_edge = TRUE, save_phys_before_they_are_chopped = FALSE) 
{
  runjunk = "\n\tmake_master_table=TRUE; plot_pieces=FALSE; cut_fossils=TRUE; fossils_older_than=0.6;\n\t\n\tinputs=BioGeoBEARS_run_object\n\tmake_master_table=TRUE\n\tplot_pieces=TRUE\n\tcut_fossils=TRUE\n\tfossils_older_than=0.1\n\tmin_dist_between_node_and_stratum_line=0.00001\n\t# testing:\n\torig_timeperiods = c(0.5, 1.9, 3.7, 5.1)\n\tsave_phys_before_they_are_chopped = FALSE\n\t"
  if (save_phys_before_they_are_chopped == TRUE) {
    phys_before_they_are_chopped = list()
    pcount = 0
  }
  
  # Store initial tree
  original_tree <- tree
  # Extract root age
  root_age <- max(phytools::nodeHeights(original_tree))
  
  ### Arrange time boundaries in timeperiods format for BioGeoBEARS
  
  if (!(tipward_time < rootward_time))
  {
    stop(paste0("Tipward time must be lower than Rootward time"))
  }

  # Initiate timeperiods vector
  orig_timeperiods <- c(tipward_time, rootward_time)
  
  # Does not need the 0 as lower limit.
  if (tipward_time == 0)
  {
    orig_timeperiods <- orig_timeperiods[-1]
  }
  
  # Need the root age or beyond as upper limit.
  if (rootward_time < root_age) # If boundary is before the root age, add an older age as upper limit
  {
    orig_timeperiods <- c(orig_timeperiods, root_age * 1.1)
  }
  if (rootward_time == root_age) # If boundary is the root age, replace with an older age as upper limit
  {
    orig_timeperiods <- orig_timeperiods[length(orig_timeperiods)] * 1.1
  }
  timeperiods <- orig_timeperiods
  
  
  if (remove_root_edge == TRUE) {
    if ("root.edge" %in% names(original_tree)) {
      txt = paste0("WARNING in chainsaw2: input tree had a 'root.edge', which crashes chainsaw2. Setting original_tree$root.edge=NULL.")
      cat("\n\n")
      cat(txt)
      cat("\n\n")
      warning(txt)
      original_tree$root.edge = NULL
    }
  }
  phy_as_it_is_chopped_down = original_tree
  orig_tr_table = prt(original_tree, printflag = FALSE, get_tipnames = TRUE, 
                      fossils_older_than = fossils_older_than)
  times_older_than_root_node_TF = orig_timeperiods > max(orig_tr_table$node_ht)
  times_younger_than_root_node_TF = orig_timeperiods < max(orig_tr_table$node_ht)
  if (sum(times_older_than_root_node_TF) >= 2) {
    txt = paste0("STOP ERROR in section_the_tree(): the timeperiods file can have ONLY ONE time older than the bottom node in your tree.  Use the function prt() to get a table of node ages for your tree. The oldest node age in your tree is: ", 
                 max(orig_tr_table$node_ht), ". The times in your timeperiods file are: ", 
                 paste(orig_timeperiods, collapse = " ", sep = ""), 
                 ".")
    cat("\n\n")
    cat(txt)
    cat("\n\n")
    stop(txt)
  }
  if (sum(times_younger_than_root_node_TF) == 0) {
    txt = paste0("STOP ERROR in section_the_tree(): the timeperiods file *HAS* to have an oldest time that is older than the bottom node in your tree.  Use the function prt() to get a table of node ages for your tree. The oldest node age in your tree is: ", 
                 max(orig_tr_table$node_ht), ". The times in your timeperiods file are: ", 
                 paste(orig_timeperiods, collapse = " ", sep = ""), 
                 ".")
    cat("\n\n")
    cat(txt)
    cat("\n\n")
    stop(txt)
  }
  for (i in 1:length(orig_timeperiods)) {
    TF = orig_tr_table$time_bp == orig_timeperiods[i]
    if (sum(TF) > 0) {
      errortxt = paste("\n\nERROR in section_the_tree(): your tree has ", 
                       sum(TF), " nodes with date ", orig_timeperiods[i], 
                       ".\nThis is a problem because you have a timeperiod boundary of: ", 
                       orig_timeperiods[i], "\nThe function doesn't know how to section a tree exactly at a node boundary.", 
                       "\nTo fix: change the timeperiod date, or edit the tree so that all nodes are more than ", 
                       min_dist_between_node_and_stratum_line, " time units from a timeperiod boundary\n(specified by 'min_dist_between_node_and_stratum_line', default min_dist_between_node_and_stratum_line=", 
                       min_dist_between_node_and_stratum_line, ").", 
                       "\n\nIf it makes you feel better, there is no way your dating of either phylogenetic or geological events is all that precise anyway.", 
                       sep = "")
      cat(errortxt)
      cat("\n\nNodes with this problem:\n\n")
      print(orig_tr_table[TF, ])
      stop("\n\nStopping on error.\n\n")
    }
    diffs = abs(orig_tr_table$time_bp - orig_timeperiods[i])
    TF = diffs < min_dist_between_node_and_stratum_line
    if (sum(TF) > 0) {
      errortxt = paste("\n\nERROR in section_the_tree(): your tree has ", 
                       sum(TF), " nodes with a date too close to your timeperiod boundary of: ", 
                       orig_timeperiods[i], ".\nThis is a problem because very short branches may cause issues with likelihood calculations, ancestral state estimation, and stochastic mapping.", 
                       "\nSee e.g. the min_branchlength option of calc_loglike_sp().", 
                       "\nTo fix: change the timeperiod date, or edit the tree so that all nodes are more than ", 
                       min_dist_between_node_and_stratum_line, " time units from a timeperiod boundary\n(specified by 'min_dist_between_node_and_stratum_line', default min_dist_between_node_and_stratum_line=", 
                       min_dist_between_node_and_stratum_line, ").", 
                       "\n\nIf it makes you feel better, there is no way your dating of either phylogenetic or geological events is all that precise anyway.", 
                       sep = "")
      cat(errortxt)
      cat("\n\nNodes with this problem:\n\n")
      print(orig_tr_table[TF, ])
      stop("\n\nStopping on error.\n\n")
    }
  }
  tipnums = 1:length(original_tree$tip.label)
  fossils_TF = orig_tr_table$time_bp[tipnums] >= fossils_older_than
  numfossils = sum(fossils_TF)
  fossil_names = original_tree$tip.label[fossils_TF]
  if (numfossils > 0) {
    if (cut_fossils == TRUE) {
      stoptxt = cat("\n\nFATAL ERROR in section_the_tree(): Your tree has ", 
                    numfossils, " fossil tips older than ", fossils_older_than, 
                    " my!\n", "But you have not turned on fossils by setting 'cut_fossils=FALSE' in section_the_tree().\n", 
                    "Fossil tipnames listed below:\n", sep = "")
      cat(stoptxt)
      print(fossil_names)
      cat("\n\nAlso: PLEASE NOTE that several times I have experienced miserable long nights due, apparently, to drop.tip producing weird tree structures, resulting in weird Newick files, without me realizing it.  The solution is usually to open the Newick file in something like FigTree, resort the branches, and save to a new Newick file.\n\n")
      stop(stoptxt)
      junk = "\n\t\t\ttr_nofossils = drop.tip(original_tree, fossil_names)\n\t\t\twrite.tree(tr_nofossils, file=\"venerid_tree_for_biogeog_v1.newick\")\n\t\t\t"
    }
    else {
      warntxt = cat("\n\nWARNING: Your tree has ", numfossils, 
                    " fossil tips older than fossils_older_than=", 
                    fossils_older_than, " my!\n", "If you actually have that many fossil tips, then everything is fine, and you can ignore this warning. If not, make sure that all fossils are older than whatever you set 'fossils_older_than' to be. If you do *not* have any fossils, then you are probably using an undated tree. This is a Very Bad Idea in general, please see 'BioGeoBEARS Mistakes To Avoid' at PhyloWiki.\n", 
                    "(default: fossils_older_than=0.6)\n", "Fossil tipnames listed below:\n", 
                    sep = "")
      warning(warntxt)
      cat(warntxt)
      cat(paste(fossil_names, collapse = "\n", sep = ""))
      cat("\n\n")
      phy_as_it_is_chopped_down = extend_tips_to_ultrametricize(obj = phy_as_it_is_chopped_down, 
                                                                age_of_root = 0, tips_end_at_this_date = NA)
    }
  }
  if (make_master_table == TRUE) {
    master_table = NULL
  }
  if (plot_pieces == TRUE) {
    plot(phy_as_it_is_chopped_down)
    axisPhylo()
  }
  tree_sections_list = NULL
  tnum = 0
  if (length(timeperiods) <= 1) {
    chainsaw_result = list()
    chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
    chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
    tmp_sorted_names_merge = paste(phy_as_it_is_chopped_down$tip.label, 
                                   collapse = ",", sep = "")
    tmp_sorted_names_split = strsplit(x = tmp_sorted_names_merge, 
                                      split = ",")[[1]]
    chainsaw_result$return_pieces_basenames[[1]] = paste(sort(tmp_sorted_names_split), 
                                                         collapse = ",", sep = "")
    attr(chainsaw_result, "class") = "chainsaw_result"
    tree_sections_list[[1]] = chainsaw_result
  }
  else {
    table_colnames = c("node", "node.type", "parent_br", 
                       "edge.length", "ancestor", "daughter_nds", "time_bp", 
                       "fossils", "label")
    SUBtable_colnames = paste("SUB", table_colnames, sep = "")
    if (make_master_table == TRUE) {
      orig_tips_table = orig_tr_table[1:length(original_tree$tip.label), 
      ]
      subtree_table = orig_tips_table
      names(subtree_table) = paste("SUB", names(subtree_table), 
                                   sep = "")
      stratum = 0
      reltimept = 0
      time_bot = 0
      time_top = 0
      piecenum = 0
      piececlass = "orig_tip"
      subtree_table = cbind(stratum, time_top, time_bot, 
                            reltimept, piecenum, piececlass, subtree_table[, 
                                                                           SUBtable_colnames])
      subtree_table$SUBnode.type = "orig_tip"
      tmp_join_table = cbind(orig_tips_table[, table_colnames], 
                             subtree_table)
      tmp_join_table
      master_table = rbind(master_table, tmp_join_table)
    }
    for (i in 1:(length(timeperiods))) {
      stratum = i
      cat("\n", i, "- top: ", orig_timeperiods[i] - timeperiods[i], 
          ", bot: ", orig_timeperiods[i], ", rel_bot: ", 
          timeperiods[i], "\n", sep = "")
      if (i == 1) {
        timepoint = timeperiods[i] - 0
      }
      else {
        timepoint = timeperiods[i]
      }
      timeperiods = timeperiods - timepoint
      timeperiods
      if (save_phys_before_they_are_chopped == TRUE) {
        phys_before_they_are_chopped[[(pcount = pcount + 
                                         1)]] = phy_as_it_is_chopped_down
      }
      if (i < length(timeperiods)) {
        chainsaw_result = chainsaw2(phy_as_it_is_chopped_down, 
                                    timepoint = timepoint, return_pieces = TRUE)
      }
      else {
        chainsaw_result = list()
        chainsaw_result$tree_to_chainsaw = phy_as_it_is_chopped_down
        chainsaw_result$return_pieces_list[[1]] = phy_as_it_is_chopped_down
        tmp_sorted_names_merge = paste(phy_as_it_is_chopped_down$tip.label, 
                                       collapse = ",", sep = "")
        tmp_sorted_names_split = strsplit(x = tmp_sorted_names_merge, 
                                          split = ",")[[1]]
        chainsaw_result$return_pieces_basenames[[1]] = paste(sort(tmp_sorted_names_split), 
                                                             collapse = ",", sep = "")
        attr(chainsaw_result, "class") = "chainsaw_result"
      }
      tree_sections_list[[(tnum = tnum + 1)]] = chainsaw_result
      if (make_master_table == TRUE) {
        tipnames_above_cutpoints = unlist(chainsaw_result$return_pieces_basenames)
        tipnames_above_cutpoints
        pos_of_1st_in_2nd = match(tipnames_above_cutpoints, 
                                  orig_tr_table$tipnames)
        pos_of_1st_in_2nd
        classes_of_pieces = sapply(X = chainsaw_result$return_pieces_list, 
                                   FUN = class)
        classes_of_pieces[classes_of_pieces == "numeric"] = "subbranch"
        classes_of_pieces[classes_of_pieces == "phylo"] = "subtree"
        phy_chopped_down_table = prt(phy_as_it_is_chopped_down, 
                                     printflag = FALSE, get_tipnames = TRUE, fossils_older_than = fossils_older_than)
        for (rownum in 1:nrow(phy_chopped_down_table)) {
          temp_tipnames = phy_chopped_down_table$tipnames[rownum]
          words = strsplit(temp_tipnames, split = ",")[[1]]
          words = sort(words)
          phy_chopped_down_table$tipnames[rownum] = paste(words, 
                                                          collapse = ",", sep = "")
        }
        reltimept = timepoint
        time_bot = orig_timeperiods[i]
        time_top = time_bot - reltimept
        for (p in 1:length(classes_of_pieces)) {
          if (classes_of_pieces[p] == "subtree") {
            tmp_subtree = chainsaw_result$return_pieces_list[[p]]
            subtree_table = prt(tmp_subtree, printflag = FALSE, 
                                get_tipnames = TRUE, fossils_older_than = fossils_older_than)
            for (rownum in 1:nrow(subtree_table)) {
              temp_tipnames = subtree_table$tipnames[rownum]
              words = strsplit(temp_tipnames, split = ",")[[1]]
              words = sort(words)
              subtree_table$tipnames[rownum] = paste(words, 
                                                     collapse = ",", sep = "")
            }
            names(subtree_table) = paste("SUB", names(subtree_table), 
                                         sep = "")
            subtree_table
            tree_piece_nodenums = subtree_table$SUBnode
            tiplabels_for_each_node_in_tree_piece = subtree_table$SUBtipnames
            pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, 
                                      orig_tr_table$tipnames)
            pos_of_1st_in_2nd
            piecenum = p
            piececlass = classes_of_pieces[p]
            subtree_table = cbind(stratum, time_top, 
                                  time_bot, reltimept, piecenum, piececlass, 
                                  subtree_table[, SUBtable_colnames])
            subtree_table
            tmp_join_table = cbind(orig_tr_table[pos_of_1st_in_2nd, 
                                                 table_colnames], subtree_table)
            tmp_join_table
            actual_heights_below_bin_top = tmp_join_table$time_bp - 
              tmp_join_table$time_top
            actual_height_lower_than_bin_top_TF = actual_heights_below_bin_top > 
              1e-10
            subtree_tip_TF = tmp_join_table$SUBnode.type == 
              "tip"
            fossil_in_subtree_TF = (actual_height_lower_than_bin_top_TF + 
                                      subtree_tip_TF) == 2
            tmp_join_table$SUBfossils[fossil_in_subtree_TF] = TRUE
            tmp_join_table$SUBtime_bp[fossil_in_subtree_TF] = actual_heights_below_bin_top[fossil_in_subtree_TF]
            new_subtree_edge_lengths = tmp_join_table$SUBedge.length[fossil_in_subtree_TF]
            tmp_join_table$SUBedge.length[fossil_in_subtree_TF] = tmp_join_table$SUBedge.length[fossil_in_subtree_TF] - 
              actual_heights_below_bin_top[fossil_in_subtree_TF]
            tmp_subtree2 = tmp_subtree
            subtree_tipnums_to_change = tmp_join_table$SUBnode[fossil_in_subtree_TF]
            edge_table_rownums_to_change = match(x = subtree_tipnums_to_change, 
                                                 table = tmp_subtree2$edge[, 2])
            tmp_subtree2$edge.length[edge_table_rownums_to_change] = tmp_join_table$SUBedge.length[fossil_in_subtree_TF]
            chainsaw_result$return_pieces_list[[p]] = tmp_subtree2
            tree_sections_list[[tnum]] = chainsaw_result
            if (is.na(tmp_join_table[1, 1]) == TRUE) {
              stoptxt = "\n\nFATAL ERROR #1 produced in section_the_tree(): NAs in tmp_join_table.\n\n"
              cat(stoptxt)
              print("i")
              print(i)
              print("p")
              print(p)
              print(tmp_join_table)
              stop(stoptxt)
            }
          }
          else {
            tmp_subbranch = chainsaw_result$return_pieces_list[[p]]
            tree_piece_nodenums = 1
            tmp_basenames = chainsaw_result$return_pieces_basenames[[p]]
            tmp_basenames2 = paste(tmp_basenames, collapse = ",", 
                                   sep = "")
            tmp_basenames3 = strsplit(x = tmp_basenames2, 
                                      split = ",")[[1]]
            tiplabels_for_each_node_in_tree_piece = paste(sort(tmp_basenames3), 
                                                          collapse = ",", sep = "")
            pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, 
                                      phy_chopped_down_table$tipnames)
            pos_of_1st_in_2nd
            subtree_table = phy_chopped_down_table
            names(subtree_table) = paste("SUB", names(subtree_table), 
                                         sep = "")
            subtree_table
            piecenum = p
            piececlass = classes_of_pieces[p]
            subtree_table = cbind(stratum, time_top, 
                                  time_bot, reltimept, piecenum, piececlass, 
                                  subtree_table[pos_of_1st_in_2nd, SUBtable_colnames])
            pos_of_1st_in_2nd = match(tiplabels_for_each_node_in_tree_piece, 
                                      orig_tr_table$tipnames)
            pos_of_1st_in_2nd
            tmp_join_table = cbind(orig_tr_table[pos_of_1st_in_2nd, 
                                                 table_colnames], subtree_table)
            tmp_join_table
            if (is.na(tmp_join_table[1, 1]) == TRUE) {
              stoptxt = "\n\nFATAL ERROR #2 produced in section_the_tree(): NAs in tmp_join_table.\n\n"
              cat(stoptxt)
              print("i")
              print(i)
              print("p")
              print(p)
              print(tmp_join_table)
              print(tiplabels_for_each_node_in_tree_piece)
              print(pos_of_1st_in_2nd)
              stop(stoptxt)
            }
          }
          master_table = rbind(master_table, tmp_join_table)
        }
      }
      phy_as_it_is_chopped_down = chainsaw_result$tree_to_chainsaw
      if (plot_pieces == TRUE) {
        plot(phy_as_it_is_chopped_down)
        axisPhylo()
      }
    }
  }
  
  # Initiate output
  output <- list
  
  # Store results
  if (save_phys_before_they_are_chopped == TRUE) {
    output$phys_before_they_are_chopped = phys_before_they_are_chopped
  }
  output$original_tree = original_tree
  output$tree_sections_list = tree_sections_list
  output$master_table = master_table
  
  return(output)
}


# print tree in hierarchical format
#######################################################
# prt
#######################################################
#' Print tree in table format
#' 
#' Learning and using APE's tree structure can be difficult and confusing because much of the information is
#' implicit.  This function prints the entire
#' tree to a table, and makes much of the implicit information explicit.  It is not particularly fast, but
#' it is useful.
#'
#' See \url{http://ape.mpl.ird.fr/ape_development.html} for the official documentation of R tree objects.
#' 
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @param printflag Should the table be printed to screen?  Default TRUE.
#' @param relabel_nodes Manually renumber the internal nodes, if desired. Default FALSE.
#' @param time_bp_digits The number of digits to print in the time_bp (time before present) column. Default=7.
#' @param add_root_edge Should a root edge be added?  Default \code{TRUE}.
#' @param get_tipnames Should the list of tipnames descending from each node be printed as a string in another column?  
#' This is slow-ish, but useful for matching up nodes between differing trees. Default \code{FALSE}.
#' @param fossils_older_than Tips that are older than \code{fossils_older_than} will be marked as \code{TRUE} in a column called \code{fossil}.
#' This is not currently set to 0, because Newick files can have slight precision issues etc. that mean not all tips quite come to zero.  You 
#' can attempt to fix this with \code{\link{average_tr_tips}} (but make sure you do not inappropriately average in fossils!!).
#' @return \code{dtf} A \code{\link[base]{data.frame}} holding the table. (Similar to the printout of a \code{\link[phylobase]{phylo4}} object.)
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{average_tr_tips}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://ape.mpl.ird.fr/ape_development.html}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1
#' 
 
prt <- function(t, printflag=TRUE, relabel_nodes = FALSE, time_bp_digits=7, add_root_edge=TRUE, get_tipnames=FALSE, fossils_older_than=0.6)
{
  # assemble beginning table
  
  # check if internal node labels exist
  if ("node.label" %in% attributes(t)$names == FALSE)
  {
    rootnum = get_nodenum_structural_root(t)
    
    new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
    t$node.label = new_node_labels
  }
  
  # or manually relabel the internal nodes, if desired
  if (relabel_nodes == TRUE)
  {
    rootnum = get_nodenum_structural_root(t)
    
    new_node_labels = paste("inNode", rootnum:(rootnum+t$Nnode-1), sep="")
    t$node.label = new_node_labels
  }
  
  labels = c(t$tip.label, t$node.label)
  ordered_nodenames = get_nodenums(t)
  #nodenums = 1:length(labels)
  node.types1 = rep("tip", length(t$tip.label))
  node.types2 = rep("internal", length(t$node.label))
  node.types2[1] = "root"
  node.types = c(node.types1, node.types2)
  
  # These are the index numbers of the edges below each node
  parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, t$edge[,2])
  #parent_edges = parent_branches
  brlen_to_parent = t$edge.length[parent_branches]
  
  parent_nodes = t$edge[,1][parent_branches]
  daughter_nodes = lapply(ordered_nodenames, get_daughters, t)
  
  # print out the structural root, if desired
  root_nodenum = get_nodenum_structural_root(t)
  tmpstr = paste("prt(t): root=", root_nodenum, "\n", sep="")
  prflag(tmpstr, printflag=printflag)
  
  levels_for_nodes = unlist(lapply(ordered_nodenames, get_level, t))
  #tmplevel = get_level(23, t)
  #print(tmplevel)
  
  
  #height above root
  hts_at_end_of_branches_aka_at_nodes = t$edge.length
  hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
  h = hts_at_end_of_branches_aka_at_nodes
  
  # times before present, below (ultrametric!) tips
  # numbers are positive, i.e. in millions of years before present
  #                       i.e. mybp, Ma
  times_before_present = get_max_height_tree(t) - h
  
  
  # fill in the ages of each node for the edges
  edge_ages = t$edge
  edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
  edge_ages[,2] = h[t$edge[,2]]	# top of branch
  
  
  # fill in the times before present of each node for the edges
  edge_times_bp = t$edge
  edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
  edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
  
  
  # If desired, get the list of all tipnames descended from a node, in alphabetical order
  if (get_tipnames == TRUE)
  {
    # Make the empty list
    list_of_clade_members_lists = rep(list(NA), length(ordered_nodenames))
    
    # Tips have only one descendant
    list_of_clade_members_lists[1:length(t$tip.label)] = t$tip.label
    list_of_clade_members_lists
    
    
    nontip_nodenums = (length(t$tip.label)+1) : length(ordered_nodenames)
    if (length(nontip_nodenums) > 1)
    {
      # More than 1 node
      nontip_nodenames = ordered_nodenames[nontip_nodenums]
      nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
      nontip_cladelists
      
      nontip_cladelists_alphabetical = sapply(X=nontip_cladelists, FUN=sort)
      nontip_cladelists_alphabetical
      
      nontip_cladelists_alphabetical_str = sapply(X=nontip_cladelists_alphabetical, FUN=paste, collapse=",")
      nontip_cladelists_alphabetical_str
      
      # Store the results
      list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
      list_of_clade_members_lists
    } else {
      # Just one node
      nontip_nodenames = ordered_nodenames[nontip_nodenums]
      nontip_cladelists = sapply(X=nontip_nodenames, FUN=get_all_daughter_tips_of_a_node, t=t)
      nontip_cladewords = unlist(sapply(X=nontip_cladelists, FUN=strsplit, split=","))
      
      nontip_cladelists_alphabetical = sort(nontip_cladewords)
      nontip_cladelists_alphabetical
      
      nontip_cladelists_alphabetical_str = paste(nontip_cladelists_alphabetical, collapse=",", sep="")
      nontip_cladelists_alphabetical_str
      
      # Store the results
      list_of_clade_members_lists[nontip_nodenums] = nontip_cladelists_alphabetical_str
      list_of_clade_members_lists			
    }
    
  }
  
  
  # Add fossils TRUE/FALSE column.  You can turn this off with fossils_older_than=NULL.
  fossils = times_before_present > fossils_older_than
  
  # Obviously, internal nodes are irrelevant and should be NA
  tmpnodenums = (length(t$tip.label)+1) : ( length(t$tip.label) + t$Nnode )
  fossils[tmpnodenums] = NA
  
  if (get_tipnames == FALSE)
  {
    # Don't put in the list of clade names
    tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels)
    
    dtf = as.data.frame(tmpdtf, row.names=NULL)
    # nd = node
    
    # edge.length is the same as brlen_2_parent
    names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label")
    
    # convert the cols from class "list" to some natural class
    dtf = unlist_dtf_cols(dtf, printflag=FALSE)
  } else {
    # Put in the list of clade names
    tmpdtf = cbind(1:length(ordered_nodenames), ordered_nodenames, levels_for_nodes, node.types, parent_branches, brlen_to_parent, parent_nodes, daughter_nodes, h, round(times_before_present, digits=time_bp_digits), fossils, labels, list_of_clade_members_lists)
    
    dtf = as.data.frame(tmpdtf, row.names=NULL)
    # nd = node
    
    # edge.length is the same as brlen_2_parent
    names(dtf) = c("node", "ord_ndname", "node_lvl", "node.type", "parent_br", "edge.length", "ancestor", "daughter_nds", "node_ht", "time_bp", "fossils", "label", "tipnames")
    
    # convert the cols from class "list" to some natural class
    dtf = unlist_dtf_cols(dtf, printflag=FALSE)		
  }
  
  
  
  
  
  
  # Add the root edge, if desired
  # (AND, only if t$root.edge exists)
  if ( (add_root_edge == TRUE) && (!is.null(t$root.edge)) )
  {
    root_row_TF = dtf$node.type == "root"
    root_edge_length = t$root.edge
    
    # Stick in this edge length
    dtf$edge.length[root_row_TF] = root_edge_length
    
    # Add the root edge length to all node heights
    dtf$node_ht = dtf$node_ht + root_edge_length
  }
  
  # print if desired
  prflag(dtf, printflag=printflag)
  
  #tree_strings = c()
  #root_str = get_node_info(root_nodenum, t)
  return(dtf)
}


#######################################################
# get_nodenums
#######################################################
#' Get the unique node numbers in a tree
#' 
#' This is a utility function for \code{\link{get_nodenum_structural_root}}.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @return \code{ordered_nodenames} The node numbers, in order.
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenum_structural_root}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah = 1
# this returns the NUMBERS identifying each node

get_nodenums <- function(t)
{
  # get just the unique node numbers from the edge list (left column: start node; right column: end node):
  nodenames = unique(c(t$edge))
  ordered_nodenames = nodenames[order(nodenames)]
  return(ordered_nodenames)
}



#######################################################
# get_indices_where_list1_occurs_in_list2
#######################################################
#' Return (first!) indices in second list matching the first list
#' 
#' This function will return one match (the first) for each item in the list; i.e. the second-list
#' index for each item in the first list.  Only the first hit in the second list is returned.
#' 
#' This is used by \code{\link{prt}}.
#'
#' @param list1 The first list. 
#' @param list2 The second list list.
#' @return \code{match_indices} The match indices.
#' @export
#' @seealso \code{\link{prt}}, \code{\link[base]{LETTERS}}, \code{\link{get_indices_where_list1_occurs_in_list2_noNA}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' list1 = c("N", "I", "C", "K")
#' list2 = LETTERS
#' get_indices_where_list1_occurs_in_list2(list1, list2)

get_indices_where_list1_occurs_in_list2 <- function(list1, list2)
{
  match_indices = match(list1, list2)
  return(match_indices)
}


#######################################################
# get_daughters
#######################################################
#' Get all the direct daughters nodes of a node
#' 
#' @param nodenum The node number to get the daughters of
#' @param t An ape phylo object
#' @return \code{daughter_nodenums} List of the daughter node numbers
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_daughters <- function(nodenum, t)
{
  daughter_edgenums = findall(nodenum, t$edge[,1])
  daughter_nodenums = t$edge[,2][daughter_edgenums]
  return(daughter_nodenums)
}


# Get indices of all matches to a list
#######################################################
# findall
#######################################################
#' Get indices of all matches to a list
#'
#' Just a handy shortcut function
#' 
#' @param what The item to find
#' @param inlist The list to search in 
#' @return \code{matching_indices} List of the matching indices
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

findall <- function(what, inlist)
{
  TFmatches = inlist == what
  indices = 1:length(inlist)
  matching_indices = indices[TFmatches]
  return(matching_indices)
}


#######################################################
# get_nodenum_structural_root
#######################################################
#' Gets the root node 
#' 
#' This function gets the root node by finding the node not in the descendants list (edge[,2]). This
#' may be more reliable than e.g. assuming length(tr$tip.label)+1.
#'
#' @param t A tree object in \code{\link[ape]{phylo}} format.
#' @param print_nodenum Print the node numbers as you go through the list? Default FALSE.
#' @return \code{root_nodenums_list} 
#' @export
#' @seealso \code{\link[ape]{phylo}}, \code{\link{get_nodenums}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' blah=1

get_nodenum_structural_root <- function(t, print_nodenum=FALSE)
{
  #numnodes = length(t$tip.label) + length(t$node.label)
  #ordered_nodes = 1:length(numnodes)
  
  ordered_nodes = get_nodenums(t)
  
  root_nodenums_list = c()
  for (n in 1:length(ordered_nodes))
  {
    tmpnode = ordered_nodes[n]
    if (tmpnode %in% t$edge[,2])
    {
      blah = TRUE
    }
    else
    {
      if (print_nodenum == TRUE)
      {
        cat("get_nodenum_structural_root(): Root nodenum = ", tmpnode, sep="")
      }
      root_nodenums_list = c(root_nodenums_list, tmpnode)
    }
  }
  return(root_nodenums_list)
}

#######################################################
# prflag
#######################################################
#' Utility function to conditionally print intermediate results
#'
#' Just a handy shortcut function, allowing other functions to optionally 
#' print, depending on the value of \code{printflag}.
#' 
#' @param x What to print.
#' @param printflag If TRUE, do the printing
#' @return nothing
#' @export
#' @seealso \code{\link{get_daughters}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples

prflag <- function(x, printflag=TRUE)
{
  # A standard function to print (or not) certain variables,
  #   based on a master printflag
  # This avoids having to comment in/out various code chunks
  #   while debugging.
  if (printflag == TRUE)
  {
    # CAT instead of PRINT if it's a string or numeric
    if (is.character(x))
    {
      cat(x, "\n", sep="")
    }
    if (is.numeric(x))
    {
      cat(x, "\n", sep="")
    } else {
      print(x)
    }
  }
  else
  {
    pass="BLAH"
  }
}


#######################################################
# get_level
#######################################################
#' Get a node's level in the tree
#'
#' Finds how many nodes deep a node is.
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @param tmplevel A starting level (the function is recursive)
#' @return \code{tmplevel} The level of the node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples

get_level <- function(nodenum, t, tmplevel=0)
{
  parent_nodenum = get_parent(nodenum, t)
  if (is.na(parent_nodenum))
  {
    #tmplevel = 0
    return(tmplevel)
  }
  else
  {
    #print(paste("parent_nodenum: ", parent_nodenum, " level: ", tmplevel, sep=""))
    tmplevel = tmplevel + 1
    tmplevel = get_level(parent_nodenum, t, tmplevel)
    return(tmplevel)
  }
  # If an error occurs
  return(NA)
}



#######################################################
# get_parent
#######################################################
#' Get the direct parent node of a node
#' 
#' @param nodenum The node number to get the parent of
#' @param t An ape phylo object
#' @return \code{parent_nodenum}The parent node number
#' @export
#' @seealso \code{\link{findall}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_parent <- function(nodenum, t)
{
  matching_edges = findall(nodenum, t$edge[,2])
  parent_nodenum = t$edge[,1][matching_edges][1]
  return(parent_nodenum)
}


#######################################################
# get_all_node_ages
#######################################################
#' Get the ages of all the nodes in the tree (above the root)
#'
#' A utility function. Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_all_node_ages <- function(obj)
{
  node_ages = dist.nodes(obj)[get_nodenum_structural_root(obj), ]
  return(node_ages)
}


#######################################################
# get_max_height_tree
#######################################################
#' Get the maximum age of all the nodes (above the root)
#'
#' I.e., the distance of the highest node above the root.  A utility function. 
#' Use of \code{\link[ape]{dist.nodes}} may be slow.
#' 
#' @param obj An ape phylo object
#' @return \code{max_height} The age (from the root) of the highest node.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_max_height_tree <- function(obj)
{
  max_height = max(get_node_ages_of_tips(obj))
  return(max_height)
}


#######################################################
# get_node_ages_of_tips
#######################################################
#' Get the ages of each tip above the root
#'
#' A utility function.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The age (from the root) of each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_node_ages_of_tips <- function(obj)
{
  TF_tips = get_TF_tips(obj)
  root_node_num = get_nodenum_structural_root(obj)
  dists_from_root = dist.nodes(obj)[root_node_num, ]
  node_ages_of_tips = dists_from_root[TF_tips]
  return(node_ages_of_tips)
}


#######################################################
# get_TF_tips
#######################################################
#' Get TRUE/FALSE for nodes being tips
#'
#' A utility function that returns \code{TRUE}/\code{FALSE} for whether or not each node is a tip.
#' 
#' @param obj An ape phylo object
#' @return \code{TF_tips} The \code{TRUE}/\code{FALSE} list for each tip.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}, \code{\link{match_list1_in_list2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_TF_tips <- function(obj)
{
  # Get TF for nodes being tips
  
  # BIG CHANGE?
  #TF_tips = match_list1_in_list2(1:length(dists_from_root), obj$tip.label)
  TF_tips = match_list1_in_list2(1:length(obj$edge), 1:length(obj$tip.label))
  #TF_tips = obj$tip.label[TF_tips_indices]
  return(TF_tips)
}


#######################################################
# match_list1_in_list2
#######################################################
#' Return TRUE for list1 items when they occur in list2
#' 
#' Return matching TRUE/FALSE values.  E.g. list1 (e.g. a big list) TRUE if it is found
#' in list2 (e.g. a smaller list)
#'
#' Utility function for %in%, when one's brain gets confused.
#' 
#' @param list1 The list of things you want to check
#' @param list2 The list of things you want to check against
#' @return \code{matchlist} The TRUE/FALSE list for list1
#' @export
#' @seealso \code{\link[base]{match}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

match_list1_in_list2 <- function(list1, list2)
{
  matchlist = list1 %in% list2
  return(matchlist)
}


#######################################################
# get_all_daughter_tips_of_a_node
#######################################################
#' Get all the daughter tips of a node
#' 
#' Like it says. Utility function.
#' 
#' @param nodenum The node to find
#' @param t A \code{\link[ape]{phylo}} tree object.
#' @return \code{temp_tips} The list of daughter tipnodes
#' @export
#' @seealso \code{\link{add_to_downpass_labels}}, \code{\link[ape]{extract.clade}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_all_daughter_tips_of_a_node <- function(nodenum, t)
{
  subtree = extract.clade(t, nodenum)
  temp_tips = subtree$tip.label
  return(temp_tips)
}

#######################################################
# unlist_dtf_cols
#######################################################
#' Unlist the columns in a data.frame
#' 
#' Utility function. What it says.
#' 
#' @param dtf Input \code{\link[base]{data.frame}}
#' @param printflag Print the results if TRUE.
#' @return \code{dtf} The data.frame, hopefully without lists for columns
#' @export
#' @seealso \code{\link[base]{unlist}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' \url{https://code.google.com/p/lagrange/}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

unlist_dtf_cols <- function(dtf, printflag=FALSE)
{
  # Sometimes cbind makes each column a list, this can screw up use/searching of
  #  the column later on.  
  # Unlist each column...
  for (i in 1:ncol(dtf))
  {
    tmpstr = paste("unlisting col: ", names(dtf)[i], "...", sep="")
    prflag(tmpstr, printflag=printflag)		
    
    # catch a possible error from unlisting
    # the number of rows needs to stay the same!!
    tmpcol = unlist(dtf[, i])
    if (length(tmpcol) != length(dtf[, i]))
    {
      tmpstr = paste("...failed! unlist(col) length=", length(tmpcol), "; nrow(dtf) = ", nrow(dtf), sep="")
      prflag(tmpstr, printflag=printflag)
    } 
    else
    {
      dtf[, i] = tmpcol
      tmpstr = paste(" ", " ", sep="")
      prflag(tmpstr, printflag=printflag)
    }
  }
  
  #dtf2 = adf(dtf)
  
  return(dtf)
}


#######################################################
# chainsaw2
#######################################################
#' Saw a tree off at a particular time before present
#' 
#' This function chops a tree like a hedge-trimmer, cutting straight across at a particular timepoint. 
#' The pieces are returned, as is the leftover tree, with branches shortened appropriately.  Pieces
#' that are mini-trees are returned as ape objects, whereas single branches are just lengths.
#'
#' This function is used during stratification, but could have other uses as well.
#' 
#' @param tr An ape phylo object.
#' @param timepoint The time at which the tree should be "chopped".
#' @param return_pieces Default TRUE, which means pieces should be returned
#' @return \code{chainsaw_result} (a list object with the pieces) or \code{tree_to_chainsaw}, just the leftover tree
#' @export
#' @seealso \code{\link{section_the_tree}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

chainsaw2 <- function(tr, timepoint=10, return_pieces=TRUE)
{
  # Take a tree and saw it off evenly across a certain timepoint.
  # This removes any tips above the timepoint, and replaces them 
  # with a single tip representing the lineage crossing
  # the timepoint (with a new tip name).
  
  # Get the tree in a table
  tr_table = prt(tr, printflag=FALSE, get_tipnames=FALSE)
  tr_table
  
  # Find the tips that are less than 10 my old and drop them
  TF_exists_more_recently_than_10mya = tr_table$time_bp < timepoint
  
  # Get the corresponding labels
  labels_for_tips_existing_more_recently_than_10mya = tr_table$label[ TF_exists_more_recently_than_10mya == TRUE ]
  
  ###########################################
  # Draft chainsaw function
  ###########################################
  # loop through the branches that cross 10 mya
  
  # get a list of the edge start/stops in the phylogeny's edges
  edge_times_bp = get_edge_times_before_present(tr)
  
  # which of these branches cross 10 mya (or whatever timepoint)?
  edges_start_earlier_than_10mya = edge_times_bp[, 1] > timepoint
  edges_end_later_than_10mya = edge_times_bp[, 2] <= timepoint
  edges_to_chainsaw = edges_start_earlier_than_10mya + edges_end_later_than_10mya == 2
  
  # then, for each of these edges, figure out how many tips exist descending from it
  # these are the nodes ABOVE the cutoff line
  nodes_to_chainsaw = tr$edge[, 2][edges_to_chainsaw]
  
  
  # Take only internal nodes (? why ?)
  numtips = length(tr$tip.label)
  #nodes_to_chainsaw = nodes_to_chainsaw[nodes_to_chainsaw > numtips]
  
  # create a copy of the tree to chainsaw
  tree_to_chainsaw = tr
  
  if (return_pieces == TRUE)
  {
    return_pieces_list = as.list(rep(NA, length(nodes_to_chainsaw)))
    return_pieces_basenames = as.list(rep(NA, length(nodes_to_chainsaw)))
    chopTable = NULL
    
  }
  
  chainsaw_table = NULL
  
  for (i in 1:length(nodes_to_chainsaw))
  {
    # If this is a tip node on the current tree, shorten the branch rather than cut it off
    if (nodes_to_chainsaw[i] <= numtips)
    {
      # Here, chainsaw is cutting an internal node, so extract the sectioned branch before you cut it down
      # (the cutting happens after the forloop)
      # (This is easy, it is just the length of the timeslab;
      #  which you should UPDATE as you move down the tree
      if (return_pieces == TRUE)
      {
        # Record the length of the branch section, and the name of that tip
        # (which is also the name of that base)
        return_pieces_list[[i]] = timepoint
        tmp_tipname = tr$tip.label[nodes_to_chainsaw[i]]
        return_pieces_basenames[[i]] = tmp_tipname
      }
      # You don't have to do anything else, the chopping of single branches is 
      # covered after the forloop
      #cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")
      
    } else {
      # Here, it's an internal node, so extract the subtree before you drop it
      tmp_subtree = extract.clade(tr, nodes_to_chainsaw[i])
      #plot(tmp_subtree, root.edge=TRUE)
      # Also, record the branchlength below this node
      branchlength_below_subtree_LCA_node = timepoint - get_max_height_tree(tmp_subtree)
      # Add this to the bottom of the subtree
      tmp_subtree$root.edge = branchlength_below_subtree_LCA_node
      #plot(tmp_subtree, root.edge=TRUE)
      
      # Record the piece, if desired
      if (return_pieces == TRUE)
      {
        # Record the length of the branch section, and the name of that tip
        # (which is also the name of that base)
        return_pieces_list[[i]] = tmp_subtree
        
        # Merge THEN split THEN sort!!
        tmp_labels_merge = paste(tmp_subtree$tip.label, collapse=",", sep="")
        tmp_labels_split = strsplit(tmp_labels_merge, split=",")[[1]]
        new_labels = sort(tmp_labels_split)
        basename_after_cutting = paste(new_labels, collapse=",", sep="")
        return_pieces_basenames[[i]] = basename_after_cutting				
      }
      
      #print(tmp_subtree$tip.label)
      
      tmp_number_of_tips = length(tmp_subtree$tip.label)
      #print(tmp_number_of_tips)
      
      # number of tips to drop = (numtips -1)
      numtips_to_drop = tmp_number_of_tips - 1 
      
      # tips_to_drop
      tmp_labels = tmp_subtree$tip.label
      
      labels_to_drop = tmp_labels[1:numtips_to_drop]
      ordered_labels_to_make_into_new_name = sort(tmp_labels)
      name_new_tip = paste(ordered_labels_to_make_into_new_name, collapse=",", sep="")
      
      # new label
      label_kept_num = length(tmp_labels)
      label_kept = tmp_labels[label_kept_num]
      #new_label = paste("CA_", label_kept, "+", numtips_to_drop, "_tips", sep="")
      new_label = name_new_tip
      tree_to_chainsaw$tip.label[tree_to_chainsaw$tip.label == label_kept] = new_label
      
      # chop off e.g. 2 of the 3 tips
      tree_to_chainsaw = drop.tip(tree_to_chainsaw, labels_to_drop)
      #cat("\ni=", i, "	ntips=", length(tree_to_chainsaw$tip.label), sep="")
    } # end else
  } # end for loop
  #plot(tree_to_chainsaw)
  #axisPhylo()
  
  tree_to_chainsaw_table = prt(tree_to_chainsaw, printflag=FALSE)
  
  tree_to_chainsaw_table_tips_TF_time_bp_LT_10my = tree_to_chainsaw_table$time_bp < timepoint
  
  
  tmp_edge_lengths =  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  
  times_bp_for_edges_to_chainsaw = tree_to_chainsaw_table$time_bp[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my]
  
  adjustment = times_bp_for_edges_to_chainsaw - timepoint
  
  revised_tmp_edge_lengths = tmp_edge_lengths + adjustment
  
  tree_to_chainsaw_table$edge.length[tree_to_chainsaw_table_tips_TF_time_bp_LT_10my] = revised_tmp_edge_lengths
  
  # revised
  ordered_nodenames = get_nodenums(tree_to_chainsaw)
  parent_branches = get_indices_where_list1_occurs_in_list2(ordered_nodenames, tree_to_chainsaw$edge[,2])
  
  NA_false = is.not.na(tree_to_chainsaw_table$edge.length)
  
  tree_to_chainsaw$edge.length[parent_branches[NA_false]] = tree_to_chainsaw_table$edge.length[NA_false]
  
  if (return_pieces == TRUE)
  {
    chainsaw_result = NULL
    chainsaw_result$tree_to_chainsaw = tree_to_chainsaw
    chainsaw_result$return_pieces_list = return_pieces_list
    chainsaw_result$return_pieces_basenames = return_pieces_basenames
    class(chainsaw_result) = "chainsaw_result"
    return(chainsaw_result)
  } else {
    return(tree_to_chainsaw)
  }
}


#######################################################
# get_edge_times_before_present
#######################################################
#' Get the times of the top and bottom of each edge
#'
#' A utility function. 
#' 
#' @param t An ape phylo object
#' @return \code{edge_times_bp} A 2-column matrix with the age (from the present) of the top and bottom of each edge.
#' @export
#' @seealso \code{\link{prt}}, \code{\link{chainsaw2}}
#' @note Go BEARS!
#' @author Nicholas J. Matzke \email{matzke@@berkeley.edu} 
#' @references
#' \url{http://phylo.wikidot.com/matzke-2013-international-biogeography-society-poster}
#' @bibliography /Dropbox/_njm/__packages/BioGeoBEARS_setup/BioGeoBEARS_refs.bib
#'   @cite Matzke_2012_IBS
#' @examples
#' test=1

get_edge_times_before_present <- function(t)
{
  #height above root
  hts_at_end_of_branches_aka_at_nodes = t$edge.length
  hts_at_end_of_branches_aka_at_nodes = get_all_node_ages(t)
  h = hts_at_end_of_branches_aka_at_nodes
  
  # times before present, below (ultrametric!) tips
  # numbers are positive, i.e. in millions of years before present
  #                       i.e. mybp, Ma
  times_before_present = get_max_height_tree(t) - h
  
  
  # fill in the ages of each node for the edges
  edge_ages = t$edge
  edge_ages[,1] = h[t$edge[,1]]	# bottom of branch
  edge_ages[,2] = h[t$edge[,2]]	# top of branch
  
  # fill in the times before present of each node for the edges
  edge_times_bp = t$edge
  edge_times_bp[,1] = times_before_present[t$edge[,1]]	# bottom of branch
  edge_times_bp[,2] = times_before_present[t$edge[,2]]	# top of branch
  
  return(edge_times_bp)
}


