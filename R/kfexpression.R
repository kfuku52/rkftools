# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

calc_complementarity = function(array1, array2, method='weighted') {
    # only works for positive values
    stopifnot(length(array1)==length(array2))
    stopifnot(all(array1>=0))
    stopifnot(all(array2>=0))
    num_item = length(array1)
    total_value = sum(array1, array2)
    complementarity = 0
    for (i in 1:num_item) {
        if (array1[i]!=array2[i]) {
            max_value = max(array1[i], array2[i])
            if (method=='weighted') {
                weight = (array1[i]+array2[i]) / total_value
                complementarity = complementarity + (abs(array1[i]-array2[i]) / max_value * weight)
            } else if (method=='independent') {
                complementarity = complementarity + ((abs(array1[i]-array2[i]) / max_value) / num_item)

            }
        }
    }
    return(complementarity)
}