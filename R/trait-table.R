
.validate_numeric_trait_table = function(trait_table, context='trait_table') {
    if (is.null(dim(trait_table)) || length(dim(trait_table)) != 2L) {
        stop(context, ' must be a matrix or data.frame.')
    }
    if (ncol(trait_table) == 0L) {
        stop(context, ' must contain at least one trait column.')
    }
    valid_columns = vapply(seq_len(ncol(trait_table)), function(i) {
        values = trait_table[,i]
        is.numeric(values) || is.integer(values) || is.logical(values)
    }, logical(1))
    if (any(!valid_columns)) {
        stop(context, ' must contain only numeric or logical trait columns.')
    }
    finite_columns = vapply(seq_len(ncol(trait_table)), function(i) {
        values = as.numeric(trait_table[,i])
        all(is.finite(values[!is.na(values)]))
    }, logical(1))
    if (any(!finite_columns)) {
        trait_names = colnames(trait_table)
        if (is.null(trait_names)) {
            trait_names = as.character(seq_len(ncol(trait_table)))
        }
        stop(
            context, ' contains non-finite trait values in column(s): ',
            paste(trait_names[!finite_columns], collapse=', ')
        )
    }
    invisible(TRUE)
}


#' Remove invariant trait columns
#'
#' @param trait_table A numeric matrix or data frame.
#' @param small_dif Finite non-negative range below which a trait is removed.
#' @param verbose Whether to report removed traits.
#' @return A list with the filtered table and removed trait names.
#' @examples
#' remove_invariant_traits(data.frame(a=c(1, 1), b=c(1, 2)))
#' @export
remove_invariant_traits = function(trait_table, small_dif=0.001, verbose=FALSE) {
    .validate_numeric_trait_table(trait_table)
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    small_dif = .normalize_finite_numeric_scalar(
        small_dif, 'small_dif', min_value=0
    )
    out_trait_table = trait_table
    num_traits = ncol(trait_table)
    trait_names = colnames(trait_table)
    if (is.null(trait_names)) {
        trait_names = as.character(seq_len(num_traits))
    }

    is_small_dif = logical(num_traits)
    for (i in seq_len(num_traits)) {
        trait_values = trait_table[,i]
        trait_values = as.numeric(trait_values)
        trait_values = trait_values[!is.na(trait_values)]
        if (length(trait_values) == 0) {
            is_small_dif[i] = TRUE
            next
        }
        trait_min = min(trait_values)
        trait_max = max(trait_values)
        is_small_dif[i] = (trait_max - trait_min < small_dif)
    }

    removed_traits = trait_names[is_small_dif]
    if (verbose && length(removed_traits)) {
        message(
            'Trait removed due to small difference (< ', small_dif,
            '): ', paste(removed_traits, collapse=', ')
        )
    } else if (verbose) {
        message('All traits passed small difference check.')
    }
    out_trait_table = out_trait_table[,!is_small_dif, drop=FALSE]
    out = list(trait_table=out_trait_table, removed_traits=removed_traits)
    return(out)
}


.replicate_base_name = function(column_name, replicate_sep) {
    split_name = strsplit(column_name, replicate_sep, fixed=TRUE)[[1]]
    if (length(split_name) <= 1L) {
        return(column_name)
    }
    paste(split_name[-length(split_name)], collapse=replicate_sep)
}


#' Merge replicate trait columns
#'
#' @param trait_table A named numeric matrix or data frame.
#' @param replicate_sep Literal separator before replicate suffixes.
#' @param verbose Whether to report detected replicate groups.
#' @return If `replicate_sep` is empty or no replicate groups are found, the
#'   original `trait_table`, retaining its matrix or data-frame class. Otherwise,
#'   a data frame with replicate groups averaged by row, ignoring missing values.
#'   An entirely missing group remains `NA`.
#' @examples
#' x = matrix(1:4, 2, dimnames=list(c("A", "B"), c("x", "y")))
#' class(merge_replicates(x, "_")) # No groups: matrix
#' colnames(x) = c("x_1", "x_2")
#' class(merge_replicates(x, "_")) # Averaged groups: data.frame
#' @export
merge_replicates = function(trait_table, replicate_sep, verbose=FALSE) {
    .validate_numeric_trait_table(trait_table)
    verbose = .normalize_single_logical_arg(verbose, 'verbose')
    replicate_sep = .normalize_single_string_arg(
        value=replicate_sep,
        arg_name='replicate_sep',
        allow_empty=TRUE
    )
    if (replicate_sep=='') {
        return(trait_table)
    }
    if (is.null(colnames(trait_table)) || any(is.na(colnames(trait_table)) | colnames(trait_table) == '')) {
        stop('trait_table must have non-empty column names in merge_replicates().')
    }
    without_reps = vapply(
        colnames(trait_table),
        .replicate_base_name,
        character(1),
        replicate_sep=replicate_sep
    )
    if (length(unique(without_reps))==ncol(trait_table)) {
        if (verbose) {
            message('No replicate was found with replicate_sep="', replicate_sep, '".')
        }
        return(trait_table)
    } else {
        if (verbose) {
            message('Replicates were found with replicate_sep="', replicate_sep,
                '". Mean values will be used.')
        }
    }
    new_cols = unique(without_reps)
    out = data.frame(matrix(ncol=length(new_cols), nrow=nrow(trait_table)))
    colnames(out) = new_cols
    rownames(out) = rownames(trait_table)
    for (new_col in colnames(out)) {
        is_col = without_reps == new_col
        values = rowMeans(trait_table[,is_col,drop=FALSE], na.rm=TRUE)
        values[is.nan(values)] = NA_real_
        out[,new_col] = values
    }
    return(out)
}


#' Sort an expression table by tree-tip order
#'
#' @param exp A data frame containing one row per tree tip.
#' @param tree A `phylo` tree.
#' @param col Name of the key column.
#' @return `exp` reordered to `tree$tip.label`.
#' @export
sort_exp = function(exp, tree, col='gene_id') {
    out_exp = exp
    if (!(col %in% colnames(out_exp))) {
        stop('Column "', col, '" is not present in exp.')
    }
    key_values = as.character(out_exp[,col])
    if (any(is.na(key_values) | trimws(key_values) == '')) {
        stop('exp contains missing/blank values in key column "', col, '".')
    }
    if (anyDuplicated(key_values)) {
        duplicated_keys = unique(key_values[duplicated(key_values)])
        stop(
            'exp contains duplicated values in key column "', col, '": ',
            paste(duplicated_keys, collapse=', ')
        )
    }
    missing_tips = setdiff(tree[['tip.label']], key_values)
    if (length(missing_tips)) {
        stop(
            'exp is missing rows for tree tip label(s): ',
            paste(missing_tips, collapse=', ')
        )
    }
    rownames(out_exp) = key_values
    out_exp = out_exp[tree[['tip.label']],,drop=FALSE]
    out_exp[,col] = tree[['tip.label']]
    rownames(out_exp) = NULL
    return(out_exp)
}


#' Get base names of expression replicate columns
#'
#' @param trait_table A table with non-empty column names.
#' @param replicate_sep Literal separator before replicate suffixes.
#' @return Unique base expression names excluding `"gene"`.
#' @export
get_expression_bases = function(trait_table, replicate_sep) {
    replicate_sep = .normalize_single_string_arg(
        value=replicate_sep,
        arg_name='replicate_sep',
        allow_empty=TRUE
    )
    if (is.null(colnames(trait_table)) || any(is.na(colnames(trait_table)) | colnames(trait_table) == '')) {
        stop('trait_table must have non-empty column names in get_expression_bases().')
    }
    if (replicate_sep=='') {
        out = unique(colnames(trait_table))
        out = out[out!='gene']
        return(out)
    }
    without_reps = vapply(
        colnames(trait_table),
        .replicate_base_name,
        character(1),
        replicate_sep=replicate_sep
    )
    out = unique(without_reps)
    out = out[out!='gene']
    return(out)
}
