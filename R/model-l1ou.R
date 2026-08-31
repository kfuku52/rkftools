# Adapters for l1ou fit objects. No l1ou installation is needed to read a fit.

.l1ou_regime_table = function(pcm_out) {
    mode = 'l1ou'
    if (is.null(colnames(pcm_out$Y))) {
        cols = 'trait1'
    } else {
        cols = colnames(pcm_out$Y)
    }
    column_names = c("regime", "node_name", "param", cols)
    regime_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
    colnames(regime_table) = column_names
    shift_conf_raw = pcm_out$shift.configuration
    shift_conf = tryCatch(
        .normalize_integerish(
            shift_conf_raw,
            'shift.configuration',
            allow_empty=TRUE
        ),
        error=function(e) integer(0)
    )
    if (length(shift_conf) != length(shift_conf_raw)) {
        stop('shift.configuration must contain integer edge indices in get_regime_table(mode="l1ou").')
    }
    if (is.null(names(shift_conf_raw))) {
        regimes = as.character(seq_along(shift_conf_raw))
    } else {
        regimes = as.character(names(shift_conf_raw))
    }
    invalid_shift_idx = shift_conf[(shift_conf < 1L) | (shift_conf > nrow(pcm_out[['tree']][['edge']]))]
    if (length(invalid_shift_idx)) {
        stop(
            'shift.configuration contains invalid edge index/indices in get_regime_table(mode="l1ou"): ',
            paste(unique(invalid_shift_idx), collapse=', ')
        )
    }
    branch_index = shift_conf
    expected_shift_count = length(branch_index)
    n_shifts_field = tryCatch(
        .normalize_integerish(pcm_out$nShifts, 'pcm_out$nShifts', min_value=0L),
        error=function(e) integer(0)
    )
    if (length(n_shifts_field) != 1L) {
        stop('pcm_out$nShifts must be a single non-negative integer in get_regime_table(mode="l1ou").')
    }
    if (n_shifts_field != expected_shift_count) {
        stop(
            'Mismatch between pcm_out$nShifts (',
            n_shifts_field,
            ') and length(pcm_out$shift.configuration) (',
            expected_shift_count,
            ') in get_regime_table(mode="l1ou").'
        )
    }
    tree = fill_node_labels(pcm_out[['tree']])
    node_names = vapply(branch_index, function(ind) {
        if (is.na(ind) || ind < 1 || ind > nrow(tree$edge)) {
            return(NA_character_)
        }
        node_values = get_node_name_by_num(phy=tree, node_num=tree$edge[ind,2])
        if (length(node_values) != 1) {
            return(NA_character_)
        }
        as.character(node_values)
    }, character(1))

    if (expected_shift_count > 0) {
        shift_values = .as_l1ou_shift_matrix(
            values=pcm_out$shift.values,
            param_name='shift.values',
            expected_rows=expected_shift_count,
            expected_cols=ncol(pcm_out$Y)
        )
        shift_means = .as_l1ou_shift_matrix(
            values=pcm_out$shift.means,
            param_name='shift.means',
            expected_rows=expected_shift_count,
            expected_cols=ncol(pcm_out$Y)
        )
        shift_values = .align_l1ou_trait_columns(
            shift_values, 'shift.values', cols
        )
        shift_means = .align_l1ou_trait_columns(
            shift_means, 'shift.means', cols
        )
        tmp_table = data.frame(
            regime=regimes,
            node_name=node_names,
            param='shift_value',
            shift_values,
            stringsAsFactors=FALSE,
            check.names=FALSE
        )
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)

        tmp_table = data.frame(
            regime=regimes,
            node_name=node_names,
            param='shift_mean',
            shift_means,
            stringsAsFactors=FALSE,
            check.names=FALSE
        )
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)
    }

    param_values = list(
        alpha=pcm_out$alpha,
        sigma2=pcm_out$sigma2,
        intercept=pcm_out$intercept,
        log_likelihood=pcm_out$logLik
    )
    for (param_name in names(param_values)) {
        values = .as_l1ou_param_vector(
            values=param_values[[param_name]],
            param_name=param_name,
            expected_cols=ncol(pcm_out$Y)
        )
        value_names = names(values)
        if (!is.null(value_names)) {
            if (anyDuplicated(value_names) || !setequal(value_names, cols)) {
                stop(
                    'Names of pcm_out$', param_name,
                    ' must match pcm_out$Y columns.'
                )
            }
            values = values[cols]
        }
        value_table = as.data.frame(
            matrix(as.numeric(values), 1, ncol(pcm_out$Y)),
            stringsAsFactors=FALSE
        )
        colnames(value_table) = cols
        tmp_table = data.frame(
            regime=NA_character_,
            node_name=NA_character_,
            param=param_name,
            value_table,
            stringsAsFactors=FALSE,
            check.names=FALSE
        )
        colnames(tmp_table) = column_names
        regime_table = rbind(regime_table, tmp_table)
    }
    rownames(regime_table) = NULL
    regime_table
}

.l1ou_leaf_table = function(pcm_out) {
    mode = 'l1ou'
    trait_cols = colnames(pcm_out$Y)
    if (is.null(trait_cols)) {
        trait_cols = paste0("trait", seq_len(ncol(pcm_out$Y)))
    }
    column_names = c("regime", 'node_name', "param", trait_cols)
    leaf_table = data.frame(matrix(0,0,ncol(pcm_out$Y)+3))
    colnames(leaf_table) = column_names
    leaf_regimes = get_leaf_regimes(pcm_out, mode)
    param2table = list(
        Y=pcm_out$Y,
        optima=pcm_out$optima,
        mu=pcm_out$mu,
        residuals=pcm_out$residuals
    )
    for (param_name in names(param2table)) {
        param_table_raw = param2table[[param_name]]
        param_col_names = colnames(param_table_raw)
        param_table = as.data.frame(param_table_raw, stringsAsFactors=FALSE)
        if (!is.null(param_col_names)) {
            colnames(param_table) = param_col_names
        } else {
            colnames(param_table) = NULL
        }
        param_row_names = rownames(param_table)
        if (!is.null(param_row_names)) {
            if (anyDuplicated(param_row_names)) {
                duplicated_rows = unique(param_row_names[duplicated(param_row_names)])
                stop(
                    'Duplicate row name(s) in pcm_out$',
                    param_name,
                    ': ',
                    paste(duplicated_rows, collapse=', ')
                )
            }
            leaf_labels = as.character(leaf_regimes[['label']])
            missing_rows = setdiff(leaf_labels, param_row_names)
            if (length(missing_rows)) {
                stop(
                    'pcm_out$',
                    param_name,
                    ' is missing row(s) for tree tip label(s): ',
                    paste(missing_rows, collapse=', ')
                )
            }
            param_table = param_table[leaf_labels,,drop=FALSE]
        } else if (nrow(param_table) != nrow(leaf_regimes)) {
            stop(
                'pcm_out$',
                param_name,
                ' must have ',
                nrow(leaf_regimes),
                ' row(s) to match tree tips.'
            )
        }
        param_table = .align_l1ou_trait_columns(param_table, param_name, trait_cols)
        tmp_table = cbind(leaf_regimes, param_name, param_table)
        colnames(tmp_table) = column_names
        leaf_table = rbind(leaf_table, tmp_table)
    }
    rownames(leaf_table) = NULL
    leaf_table
}

.as_l1ou_shift_matrix = function(values, param_name, expected_rows, expected_cols) {
    mat = as.matrix(values)
    if (expected_rows == 0L) {
        return(matrix(numeric(0), nrow=0, ncol=expected_cols))
    }
    if (is.null(dim(mat)) || length(dim(mat)) != 2L) {
        stop(
            'pcm_out$',
            param_name,
            ' must be a matrix/data.frame with ',
            expected_rows,
            ' row(s) and ',
            expected_cols,
            ' column(s) in get_regime_table(mode="l1ou").'
        )
    }
    if (nrow(mat) != expected_rows || ncol(mat) != expected_cols) {
        stop(
            'pcm_out$',
            param_name,
            ' has invalid dimensions in get_regime_table(mode="l1ou"): expected ',
            expected_rows,
            'x',
            expected_cols,
            ', got ',
            nrow(mat),
            'x',
            ncol(mat),
            '.'
        )
    }
    mat
}


.as_l1ou_param_vector = function(values, param_name, expected_cols) {
    original_names = names(values)
    vec = as.vector(values)
    if (!is.null(original_names) && length(original_names) == length(vec)) {
        names(vec) = original_names
    }
    if (expected_cols == 0L) {
        return(numeric(0))
    }
    if (length(vec) == expected_cols) {
        return(vec)
    }
    if (length(vec) == 1L) {
        return(unname(rep(vec, expected_cols)))
    }
    stop(
        'pcm_out$',
        param_name,
        ' must have length 1 or ',
        expected_cols,
        ' in get_regime_table(mode="l1ou"), got ',
        length(vec),
        '.'
    )
}


.align_l1ou_trait_columns = function(param_table, param_name, trait_cols) {
    if (ncol(param_table) != length(trait_cols)) {
        stop(
            'pcm_out$',
            param_name,
            ' must have ',
            length(trait_cols),
            ' column(s) to match pcm_out$Y.'
        )
    }

    param_cols = colnames(param_table)
    if (is.null(param_cols)) {
        colnames(param_table) = trait_cols
        return(param_table)
    }
    if (anyDuplicated(param_cols)) {
        duplicated_cols = unique(param_cols[duplicated(param_cols)])
        stop(
            'Duplicate column name(s) in pcm_out$',
            param_name,
            ': ',
            paste(duplicated_cols, collapse=', ')
        )
    }

    missing_cols = setdiff(trait_cols, param_cols)
    extra_cols = setdiff(param_cols, trait_cols)
    if (length(missing_cols) || length(extra_cols)) {
        stop(
            'Column names of pcm_out$',
            param_name,
            ' must match pcm_out$Y. Missing: ',
            ifelse(length(missing_cols), paste(missing_cols, collapse=', '), 'none'),
            '. Extra: ',
            ifelse(length(extra_cols), paste(extra_cols, collapse=', '), 'none'),
            '.'
        )
    }

    param_table[, trait_cols, drop=FALSE]
}
