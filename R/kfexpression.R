# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

calc_complementarity = function(array1, array2, method='weighted') {
    # only works for positive values
    stopifnot(length(array1)==length(array2))
    stopifnot(all(array1>=0))
    stopifnot(all(array2>=0))
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
