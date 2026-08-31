
#' Count independent foreground lineages
#'
#' @param tree A `phylo` tree with unique tips.
#' @param trait A data frame containing `species` and one or more binary traits.
#' @return A named list of independent foreground-lineage counts.
#' @examples
#' tree <- ape::read.tree(text="((A:1,B:1):1,C:1);")
#' trait <- data.frame(species=c("A", "B", "C"), foreground=c(1, 1, 0))
#' count_foreground_lineage(tree, trait)
#' @export
count_foreground_lineage = function(tree, trait) {
    num_fg_lineage = list()
    .validate_phylo_input(tree, context='tree', unique_tips=TRUE)
    if (!("species" %in% colnames(trait))) {
        stop('trait must contain a "species" column.')
    }
    trait_cols = setdiff(colnames(trait), "species")
    if (!length(trait_cols)) {
        stop('trait must contain at least one foreground trait column.')
    }
    species_names = as.character(trait[['species']])
    if (any(is.na(species_names) | trimws(species_names) == '')) {
        stop('trait$species must contain only non-missing, non-empty names.')
    }
    if (anyDuplicated(species_names)) {
        stop(
            'trait$species contains duplicated names: ',
            paste(unique(species_names[duplicated(species_names)]), collapse=', ')
        )
    }
    unknown_species = setdiff(stats::na.omit(species_names), tree[['tip.label']])
    if (length(unknown_species)) {
        stop(
            'trait$species contains names not present in tree tip labels: ',
            paste(unique(unknown_species), collapse=', ')
        )
    }
    missing_species = setdiff(tree[['tip.label']], species_names)
    if (length(missing_species)) {
        stop(
            'trait is missing rows for tree tip label(s): ',
            paste(missing_species, collapse=', ')
        )
    }
    index = .build_phy_index(tree, context='tree')
    node_nums = seq_len(index[['max_node']])
    num_tip = index[['num_tip']]
    children_by_parent = index[['children']]
    parent_by_child = index[['parent']]
    root_num = index[['root']]
    traversal_order = index[['preorder']]
    for (trait_col in trait_cols) {
        trait_values = trait[[trait_col]]
        if (is.factor(trait_values)) {
            trait_values = as.character(trait_values)
        }
        if (is.logical(trait_values)) {
            trait_values = as.numeric(trait_values)
        } else {
            suppressWarnings(trait_values <- as.numeric(trait_values))
        }
        is_binary = all(is.na(trait_values) | (trait_values %in% c(0, 1)))
        if (!is_binary) {
            warning(
                trait_col,
                ': count_foreground_lineage() supports 0/1 traits. Returning NA.',
                call.=FALSE
            )
            num_fg_lineage[[trait_col]] = NA
            next
        }
        if (anyNA(trait_values)) {
            warning(
                trait_col,
                ': missing foreground states make the lineage count undefined. Returning NA.'
            )
            num_fg_lineage[[trait_col]] = NA
            next
        }
        fg_spp = species_names[(!is.na(trait_values)) & (trait_values == 1)]
        is_fg_only_clade = logical(length(node_nums))
        is_fg_only_clade[seq_len(num_tip)] = tree[['tip.label']] %in% fg_spp
        for (node_num in rev(traversal_order)) {
            if (node_num <= num_tip) {
                next
            }
            children = children_by_parent[[node_num]]
            is_fg_only_clade[[node_num]] = length(children) > 0L &&
                all(is_fg_only_clade[children])
        }
        is_fg_stem = logical(length(node_nums))
        is_fg_stem[[root_num]] = is_fg_only_clade[[root_num]]
        nonroot_nodes = setdiff(node_nums, root_num)
        nonroot_parents = unname(parent_by_child[nonroot_nodes])
        if (anyNA(nonroot_parents)) {
            stop('Invalid tree topology in count_foreground_lineage().')
        }
        is_fg_stem[nonroot_nodes] = is_fg_only_clade[nonroot_nodes] &
            !is_fg_only_clade[nonroot_parents]
        num_fg_lineage[[trait_col]] = sum(is_fg_stem)
    }
    return(num_fg_lineage)
}
