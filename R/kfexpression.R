# Expression-profile utilities.

#' Calculate expression-profile complementarity
#'
#' Compares two non-negative numeric profiles using weighted or independent
#' relative differences.
#'
#' @param array1,array2 Equal-length finite, non-negative numeric vectors.
#' @param method Either `"weighted"` or `"independent"`.
#' @return A numeric complementarity score.
#' @examples
#' calc_complementarity(c(1, 2, 3), c(1, 2, 0))
#' @export
calc_complementarity = function(array1, array2, method='weighted') {
    # only works for positive values
    if (length(array1) != length(array2)) {
        stop('array1 and array2 must have equal length in calc_complementarity().')
    }
    if (!is.numeric(array1) || !is.numeric(array2)) {
        stop('array1 and array2 must be numeric in calc_complementarity().')
    }
    if (anyNA(array1) || anyNA(array2) || any(!is.finite(array1)) || any(!is.finite(array2))) {
        stop('array1 and array2 must contain only finite, non-missing values.')
    }
    if (any(array1 < 0) || any(array2 < 0)) {
        stop('array1 and array2 must contain only non-negative values.')
    }
    method_name = match.arg(method, c('weighted', 'independent'))
    num_item = length(array1)
    if (num_item == 0) {
        return(0)
    }

    abs_diff = abs(array1 - array2)
    is_different = abs_diff != 0
    if (!any(is_different)) {
        return(0)
    }

    max_value = pmax(array1, array2)
    contribution = abs_diff[is_different] / max_value[is_different]
    if (method_name == 'weighted') {
        total_value = sum(array1, array2)
        weights = (array1[is_different] + array2[is_different]) / total_value
    } else {
        weights = rep(1 / num_item, sum(is_different))
    }
    complementarity = sum(contribution * weights)
    return(complementarity)
}
