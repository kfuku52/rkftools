# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

remove_invariant_traits = function(trait_table, small_dif=0.001) {
    is_small_dif = c()
    for (i in 1:length(trait_table)) {
        trait_min = min(trait_table[,i])
        trait_max = max(trait_table[,i])
        is_small_dif = c(is_small_dif, trait_max - trait_min < small_dif)
    }
    removed_traits = original_traits[is_small_dif]
    if (length(removed_traits)) {
        cat("Trait removed due to small difference (<", small_dif, '):',  removed_traits, '\n')
    } else {
        cat('All traits passed small difference check.\n')
    }
    trait_table = trait_table[,!is_small_dif]
    out = list(trait_table=trait_table, removed_traits=removed_traits)
    return(out)
}

get_high_similarity_clades = function(tree, trait_table, method, threshold, verbose=FALSE, num_test=0) {
    num_all_leaves = length(tree$tip.label)
    collapse_node_nums = c()
    root_num = get_root_num(tree)
    subroot_nums = get_children_num(tree, root_num)
    next_node_nums = subroot_nums
    for (nnn in next_node_nums) {
        if (length(get_tip_labels(tree, nnn))==1) {
            next_node_nums = next_node_nums[next_node_nums!=nnn]
        }
    }
    num_processed_node = 0
    while (length(next_node_nums) > 0) {
        current_node_num = next_node_nums[1]
        tip_labels = get_tip_labels(tree, current_node_num)
        num_leaves = length(tip_labels)
        leaf_combinations = expand.grid(c1=1:num_leaves, c2=1:num_leaves)
        leaf_combinations = subset(leaf_combinations, c1 < c2)
        rownames(leaf_combinations) = NULL
        num_combinations = nrow(leaf_combinations)
        if (verbose) {
            cat('number of next_node_nums =', length(next_node_nums), 'number of leaf combinations =', num_combinations, '\n')
        }
        min_similarity = 1
        do_collapse = TRUE
        for (i in 1:num_combinations) {
            current_similarity = 1
            leaf_c1 = tip_labels[leaf_combinations[i,'c1']]
            leaf_c2 = tip_labels[leaf_combinations[i,'c2']]
            trait_c1 = unlist(trait_table[leaf_c1,])
            trait_c2 = unlist(trait_table[leaf_c2,])
            if (method=='complementarity') {
                current_similarity = 1 - calc_complementarity(trait_c1, trait_c2, method='weighted')
            } else {
                current_similarity = suppressWarnings(cor(trait_c1, trait_c2, method=method))
                current_similarity = ifelse(is.na(current_similarity), 1, current_similarity)
            }
            min_similarity = min(min_similarity, current_similarity)
            if (min_similarity < threshold) {
                do_collapse = FALSE
                if (verbose) {
                    cat('A low similarity (', min_similarity, ') found at', i, 'th leaf combinations (', num_combinations, ')\n')
                }
                break
            }
        }
        if (length(next_node_nums)==1) {
            next_node_nums = c()
        } else {
            next_node_nums = next_node_nums[2:length(next_node_nums)]
        }
        if (do_collapse) {
            cat('node_num =', current_node_num, 'size =', num_leaves, 'Min similarity =', min_similarity, '\n')
            collapse_node_nums = c(collapse_node_nums, current_node_num)
       } else {
            children_nums = get_children_num(tree, current_node_num)
            children_nums = na.omit(children_nums)
            children_nums = children_nums[children_nums > num_all_leaves]
            next_node_nums = c(next_node_nums, children_nums)
            if (verbose) {
                cat('node_num =', current_node_num, 'size =', num_leaves, 'Min similarity =', min_similarity, '\n')
            }
        }
        num_processed_node = num_processed_node + 1
        if (num_processed_node%%100==0) {
            cat('processed', num_processed_node, 'nodes\n')
        }
        if ((num_test!=0)&(num_processed_node==num_test)) {
            cat('Reaching num_test, Exiting at the ', num_test, ' th processed node.\n')
            break
        }
    }
    return(collapse_node_nums)

}

collapse_clades = function(tree, trait_table, collapse_node_nums) {
    cat('number of leaves before collapse =', length(tree$tip.label), '\n')
    trait_table = trait_table[tree$tip.label,]
    retain_tips = list()
    drop_tips = list()
    collapse_leaf_names = list()
    for (cnn in collapse_node_nums) {
        collapse_leaf_names[[as.character(cnn)]] = get_tip_labels(tree, cnn)
        retain_tips[[cnn]] = collapse_leaf_names[[as.character(cnn)]][1]
        drop_tips[[cnn]] = collapse_leaf_names[[as.character(cnn)]][2:length(collapse_leaf_names[[as.character(cnn)]])]
    }
    out_tree = tree
    out_trait = trait_table
    for (cnn in collapse_node_nums) {
        out_tree = ape::drop.tip(out_tree, tip=drop_tips[[cnn]], trim.internal=TRUE)
        out_tree$tip.label[out_tree$tip.label==retain_tips[[cnn]]] = cnn
        subtree = ape::extract.clade(tree, cnn)
        subtree_traits = trait_table[subtree$tip.label,]
        subtree_root_states = c()
        for (i in 1:ncol(trait_table)) {
            ace_result = ape::ace(subtree_traits[,i], subtree, type='continuous', method='REML', CI=FALSE, model='BM')
            subtree_root_state = ace_result$ace[names(ace_result$ace)==get_root_num(subtree)]
            subtree_root_states = c(subtree_root_states, subtree_root_state)
        }
        df_subtree_root = data.frame(subtree_root_states)
        out_trait = out_trait[!(rownames(out_trait) %in% subtree$tip.label),]
        out_trait = rbind(out_trait, subtree_root_states)
        rownames(out_trait)[nrow(out_trait)] = cnn
    }
    out_trait = out_trait[out_tree$tip.label,]
    cat('number of leaves after collapse =', length(out_tree$tip.label), '\n')
    out = list(tree=out_tree, trait=out_trait, collapse_leaf_names=collapse_leaf_names)
    return(out)
}

map_node_num = function(tree_original, tree_collapsed, collapse_leaf_names, verbose=FALSE) {
    # tree_collapsed should be a collapsed tree
    stopifnot(length(tree_original$tip.label)>=length(tree_collapsed$tip.label))
    df = data.frame(matrix(NA, max(tree_original$edge[,1]), 2))
    colnames(df) = c('tree_original', 'tree_collapsed')
    tree_original_node_nums = 1:max(tree_original$edge[,1])
    tree_collapsed_node_nums = 1:max(tree_collapsed$edge[,1])
    while (length(tree_original_node_nums)>0) {
        nn1 = tree_original_node_nums[1]
        nn1_leaf_names = sort(get_tip_labels(tree_original, nn1))
        for (i in 1:length(tree_collapsed_node_nums)) {
            nn2 = tree_collapsed_node_nums[i]
            nn2_leaf_names = sort(get_tip_labels(tree_collapsed, nn2))
            collapsed_node_nums = nn2_leaf_names[as.character(nn2_leaf_names) %in% names(collapse_leaf_names)]
            any_collapse = (length(collapsed_node_nums)>0)
            is_collapse_nn2 = (length(nn2_leaf_names)==1)&(any_collapse)
            if (any_collapse) {
                nn2_leaf_names = nn2_leaf_names[!nn2_leaf_names %in% collapsed_node_nums]
                for (cnn in collapsed_node_nums) {
                    nn2_leaf_names = c(nn2_leaf_names, collapse_leaf_names[[as.character(cnn)]])
                }
                nn2_leaf_names = sort(nn2_leaf_names)
            }
            if (is_collapse_nn2) {
                flag = all(nn1_leaf_names %in% nn2_leaf_names)
            } else if (length(nn1_leaf_names)==length(nn2_leaf_names)) {
                flag = identical(nn1_leaf_names, nn2_leaf_names)
            } else {
                flag = FALSE
            }
            if (verbose) {
                cat('nn1:', nn1, 'nn2:', nn2,
                    'collapse match:', all(nn1_leaf_names %in% nn2_leaf_names),
                    'noncollapse match:', all(nn1_leaf_names==nn2_leaf_names),
                    'is_collapse_nn2:', is_collapse_nn2, '\n')
            }
            if (flag) {
                df[nn1,] = c(nn1,nn2)
                tree_original_node_nums = tree_original_node_nums[-1]
                if (!is_collapse_nn2) {
                    tree_collapsed_node_nums = tree_collapsed_node_nums[-i]
                    tree_collapsed_node_nums = tree_collapsed_node_nums[tree_collapsed_node_nums!=nn2]
                }
                break
            }
        }
        if (!flag) {
            cat('Cannot find the original tree node:', nn1, '\n')
            cat('tree_original_node_nums:', tree_original_node_nums, '\n')
            cat('tree_collapsed_node_nums:', tree_collapsed_node_nums, '\n')
            break
        }
    }
    num_match_before = nrow(df)
    df = df[(!is.na(df[['tree_original']])),]
    df = df[(!is.na(df[['tree_collapsed']])),]
    num_match_after = nrow(df)
    if (num_match_before!=num_match_after) {
        num_na = num_match_before-num_match_after
        cat('Warning:', num_na, 'NA are produced in map_node_num().\n')
    }
    return(df)
}

get_tree_table = function(pcm_out, mode) {
    if (mode=='l1ou') {
        leaves = pcm_out$tree$tip.label
        leaf2sp = function(leaf_name) {
            leaf_split = strsplit(leaf_name, "_")
            sp = paste(leaf_split[[1]][1], leaf_split[[1]][2], sep="_")
            return(sp)
        }
        spp = unique(sapply(leaves, leaf2sp))
        if (is.null(names(pcm_out$shift.configuration))) {
            shift_conf = c("0", as.character(1:length(pcm_out$shift.configuration)))
        } else {
            shift_conf = c("0", attr(pcm_out$shift.configuration, "name"))
        }
        num_regime = length(unique(shift_conf))
        num_conv_regime = length(unique(shift_conf[duplicated(shift_conf)]))
        tree_table = data.frame(
            num_shift = pcm_out$nShifts,
            num_regime = num_regime,
            num_conv_regime = num_conv_regime,
            num_uniq_regime = num_regime - num_conv_regime,
            num_species = length(spp),
            num_leaf = length(leaves),
            model_score = pcm_out$score,
            stringsAsFactors = FALSE
        )
    } else if (mode=='PhylogeneticEM') {
        pp = PhylogeneticEM::params_process(pcm_out)
        df = data.frame(matrix(NA,1,4))
        colnames(df) = c('num_shift','log_likelihood','num_species','num_leaf')
        df$num_shift = ncol(pp[['shifts']][['values']])
        df$log_likelihood = attr(pp, 'log_likelihood')
        df$num_species = length(unique(suppressWarnings(leaf2species(leaf_names=pcm_out[['phylo']][['tip.label']]))))
        df$num_leaf = length(pcm_out[['phylo']][['tip.label']])
        tree_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
        return(tree_table)
}

get_regime_table = function(pcm_out, mode) {
    if (mode=='l1ou') {
        column_names = c("regime", "node_name", "param", colnames(pcm_out$Y))
        regime_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
        colnames(regime_table) = column_names
        shift_conf = sort(pcm_out$shift.configuration, decreasing=TRUE)
        if (is.null(names(shift_conf))) {
            regimes = 1:length(shift_conf)
        } else {
            regimes = as.integer(names(shift_conf))
        }
        branch_index = shift_conf
        node_names = c()
        for (ind in branch_index) {
            node_names = c(node_names, get_node_name_by_num(phy=adj_data$tree, node_num=adj_data$tree$edge[ind,2]))
        }

        if (pcm_out$nShifts > 0) {
            tmp_table = cbind(regimes, node_names, "shift_value", pcm_out$shift.values)
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)

            tmp_table = cbind(regimes, node_names, "shift_mean", pcm_out$shift.means)
            colnames(tmp_table) = column_names
            regime_table = rbind(regime_table, tmp_table)
        }

        tmp_table = cbind(NA, NA, "alpha", matrix(pcm_out$alpha, 1, length(pcm_out$alpha)))
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)

        tmp_table = cbind(NA, NA, "sigma2", matrix(pcm_out$sigma2, 1, length(pcm_out$sigma2)))
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)

        tmp_table = cbind(NA, NA, "intercept", matrix(pcm_out$intercept, 1, length(pcm_out$intercept)))
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)

        tmp_table = cbind(NA, NA, "log_likelihood", matrix(pcm_out$logLik, 1, length(pcm_out$logLik)))
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)
    } else if (mode=='PhylogeneticEM') {
        pp = PhylogeneticEM::params_process(pcm_out)
        traits = rownames(pcm_out[['Y_data']])
        num_trait = length(traits)
        shift_node_indices = pp[['shifts']][['edges']]
        shift_node_nums = pcm_out[['phylo']]$edge[shift_node_indices,2]
        num_shift = length(shift_node_nums)
        df = data.frame(matrix(NA, num_shift+3, num_trait+3))
        colnames(df) = c('regime','node_name','param',traits)
        row=1
        if (num_shift) {
            df[row:(row+num_shift-1),'regime'] = 1:num_shift
            df[row:(row+num_shift-1),'node_name'] = get_node_name_by_num(pcm_out[['phylo']], shift_node_nums)
            df[row:(row+num_shift-1),'param'] = 'shift_value'
            df[row:(row+num_shift-1),traits] = t(pp[['shifts']][['values']])
        }
        row = row+num_shift - 1
        row = row+1; df[row,] = c(NA, NA, 'alpha', diag(as.matrix(pp$selection.strength)))
        row = row+1; df[row,] = c(NA, NA, 'sigma2', diag(as.matrix(pp$variance)))
        row = row+1; df[row,] = c(NA, NA, 'intercept', pp$optimal.value)
        regime_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    rownames(regime_table) = NULL
    return(regime_table)
}

get_leaf_regimes = function(pcm_out, mode) {
    if (mode=='l1ou') {
        tree = pcm_out$tree
        shift_nodes = c()
        leaf_regimes = data.frame(regime=0, leaf=rownames(pcm_out$Y))
        shift_conf = sort(pcm_out$shift.configuration, decreasing=TRUE)
        if ((length(shift_conf)>0)&(is.null(names(shift_conf)))) {
            names(shift_conf) = 1:length(shift_conf)
        }
        for (node_index in shift_conf) {
            node_name = get_node_name_by_num(phy=tree, node_num=tree$edge[node_index,2])
            regime = ifelse(is.null(names(shift_conf)[shift_conf==node_index]), 1, names(shift_conf)[shift_conf==node_index])
            for (subtree in subtrees(tree)) {
                if (subtree$node.label[1]==node_name) {
                    table_leaves = as.character(leaf_regimes$leaf)
                    subtree_leaves = as.character(subtree$tip.label)
                    leaf_regimes[table_leaves %in% subtree_leaves, "regime"] = regime
                    break
                }
            }
            for (leaf in tree$tip.label) {
                if (leaf==node_name) {
                    leaf_regimes[leaf_regimes$leaf==leaf, "regime"] = regime
                    break
                }
            }
        }
    } else if (mode=='PhylogeneticEM') {
        tree = pcm_out[['phylo']]
        pp = PhylogeneticEM::params_process(pcm_out)
        shift_nodes = pp[['shifts']][['edges']]
        leaf_regimes = data.frame(regime=0, leaf=tree$tip.label)
        shift_conf = shift_nodes
        if (!is.null(shift_conf)) {
            if (is.null(names(shift_conf))) {
                names(shift_conf) = 1:length(shift_conf)
            }
        }
        shift_conf = shift_conf[order(tree$edge[shift_conf,1], decreasing=FALSE)]
        for (node_index in shift_conf) {
            node_name = get_node_name_by_num(phy=tree, node_num=tree$edge[node_index,2])
            regime = ifelse(is.null(names(shift_conf)[shift_conf==node_index]), 1, names(shift_conf)[shift_conf==node_index])
            for (subtree in subtrees(tree)) {
                if (subtree$node.label[1]==node_name) {
                    table_leaves = as.character(leaf_regimes$leaf)
                    subtree_leaves = as.character(subtree$tip.label)
                    leaf_regimes[table_leaves %in% subtree_leaves, "regime"] = regime
                    break
                }
            }
            for (leaf in tree$tip.label) {
                if (leaf==node_name) {
                    leaf_regimes[leaf_regimes$leaf==leaf, "regime"] = regime
                    break
                }
            }
        }
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    return(leaf_regimes)
}

get_leaf_table = function(pcm_out, mode) {
    if (mode=='l1ou') {
        column_names = c("regime", "leaf", "param", colnames(pcm_out$Y))
        leaf_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
        colnames(leaf_table) = column_names
        leaf_regimes = get_leaf_regimes(pcm_out, mode)

        tmp_table = cbind(leaf_regimes, "Y", pcm_out$Y)
        colnames(tmp_table) = column_names
        leaf_table = rbind(leaf_table, tmp_table)

        tmp_table = cbind(leaf_regimes, "optima", pcm_out$optima)
        colnames(tmp_table) = column_names
        leaf_table = rbind(leaf_table, tmp_table)

        tmp_table = cbind(leaf_regimes, "mu", pcm_out$mu)
        colnames(tmp_table) = column_names
        leaf_table = rbind(leaf_table, tmp_table)

        tmp_table = cbind(leaf_regimes, "residuals", pcm_out$residuals)
        colnames(tmp_table) = column_names
        leaf_table = rbind(leaf_table, tmp_table)
    } else if (mode=='PhylogeneticEM') {
        num_leaf = length(pcm_out[['phylo']][['tip.label']])
        traits = rownames(pcm_out[['Y_data']])
        num_trait = length(traits)
        params = c('imputed','expectations')
        num_param = length(params)
        df_leaf_regimes = get_leaf_regimes(pcm_out, mode='PhylogeneticEM')
        regimes = df_leaf_regimes$regime
        df = data.frame(matrix(NA, num_leaf*num_param, num_trait+3))
        colnames(df) = c('regime','leaf','param', rownames(pcm_out[['Y_data']]))
        ind_start = 1
        for (param in params) {
            ind_end = ind_start + num_leaf - 1
            df[ind_start:ind_end,'regime'] = regimes
            df[ind_start:ind_end,'leaf'] = pcm_out[['phylo']][['tip.label']]
            df[ind_start:ind_end,'param'] = param
            df[ind_start:ind_end,traits] = t(PhylogeneticEM::imputed_traits(pcm_out, trait=1:num_trait, where='tips', what=param))
            ind_start = ind_end + 1
        }
        leaf_table = df
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    rownames(leaf_table) = NULL
    return(leaf_table)
}

get_bootstrap_table = function(pcm_out, bootstrap_result) {
    if (mode=='l1ou') {
        node_names = c(pcm_out$tree$tip.label, pcm_out$tree$node.label)
        nleaf = length(pcm_out$tree$tip.label)
        column_names = c("node_name", "bootstrap_support")
        bp_table = data.frame(matrix(0,length(node_names),length(column_names)))
        colnames(bp_table) = column_names
        bp_table["node_name"] = node_names
        bp_table["bootstrap_support"] = c(bootstrap_result$detection.rate[1:nleaf], NA, bootstrap_result$detection.rate[(nleaf+1):length(bootstrap_result$detection.rate)])
    } else {
        cat('mode "', mode, '" is not supported.\n')
    }
    return(bp_table)
}

tree_table_collapse2original = function(tree_table, tree_original) {
    num_leaf = length(tree_original$tip.label)
    num_species = length(unique(leaf2species(tree_original$tip.label)))
    tree_table[,'num_leaf'] = num_leaf
    tree_table[,'num_species'] = num_species
    return(tree_table)
}

get_deepest_node_num = function(tree, node_nums) {
    depth_values = ape::node.depth(tree)
    target_depth_values = depth_values[node_nums]
    deepest_node_num = node_nums[target_depth_values==max(target_depth_values)]
    return(deepest_node_num)
}

regime_table_collapse2original = function(regime_table, tree_original, tree_collapsed, node_num_mapping) {
    node_names_collapsed = unlist(unique(na.omit(regime_table['node_name'])))
    for (node_name_collapsed in node_names_collapsed) {
        if (node_name_collapsed %in% c(tree_collapsed$tip.label, tree_collapsed$node.label)) {
            node_num_collapsed = get_node_num_by_name(tree_collapsed, node_name_collapsed)
            node_num_originals = node_num_mapping[(node_num_mapping['tree_collapsed']==node_num_collapsed),'tree_original']
            node_num_original = get_deepest_node_num(tree_original, node_num_originals)
            node_name_original = get_node_name_by_num(tree_original, node_num_original)
            is_target = (regime_table['node_name']==node_name_collapsed)
            is_target[is.na(is_target)] = FALSE
            regime_table[is_target,'node_name'] = node_name_original
        }
    }
    rownames(df) = NULL
    return(regime_table)
}

leaf_table_collapse2original = function(leaf_table, tree_original, tree_collapsed, node_num_mapping) {
    num_leaf_original = length(tree_original$tip.label)
    node_names_collapsed = unlist(unique(na.omit(leaf_table['leaf'])))
    params = unlist(unique(na.omit(leaf_table['param'])))
    df = data.frame()
    for (param in params) {
        for (node_name_collapsed in node_names_collapsed) {
            if (node_name_collapsed %in% c(tree_collapsed$tip.label, tree_collapsed$node.label)) {
                node_num_collapsed = get_node_num_by_name(tree_collapsed, node_name_collapsed)
                node_num_originals = node_num_mapping[(node_num_mapping['tree_collapsed']==node_num_collapsed),'tree_original']
                for (node_num_original in node_num_originals) {
                    if (node_num_original <= num_leaf_original) {
                        node_name_original = get_node_name_by_num(tree_original, node_num_original)
                        conditions = (leaf_table['param']==param)&(leaf_table['leaf']==node_name_collapsed)
                        row = leaf_table[conditions,]
                        row[,'leaf'] = node_name_original
                        df = rbind(df, row)
                    }
                }
            }
        }
    }
    rownames(df) = NULL
    return(df)
}

restore_imputed_leaves = function(leaf_table, original_trait_table) {
    traits = colnames(original_trait_table)
    leaf_names = leaf_table$leaf
    for (leaf_name in leaf_names) {
        if (leaf_name %in% rownames(original_trait_table)) {
            conditions = ((leaf_table$leaf==leaf_name)&(leaf_table$param=='imputed'))
            leaf_table[conditions, traits] = original_trait_table[leaf_name, traits]
        }
    }
    return(leaf_table)
}