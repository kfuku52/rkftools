# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

calc_complementarity = function(array1, array2) {
    stopifnot(length(array1)==length(array2))
    num_item = length(array1)
    sum_dif = 0
    for (i in 1:num_item) {
        if (array1[i]!=array2[i]) {
            sum_dif = sum_dif + (abs(array1[i]-array2[i]) / max(array1[i], array2[i]))
        }
    }
    normalized_dif = sum_dif / num_item
    return(normalized_dif)
}

