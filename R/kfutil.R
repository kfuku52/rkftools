# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/20/18

get_parsed_args = function(args, print=TRUE) {
    split = strsplit(sub("^--", "", args), "=")
    parsed = list()
    for (i in 1:length(split)) {
        param = split[[i]][1]
        value = split[[i]][2]
        if (!is.na(suppressWarnings(as.numeric(value)))) {
            value = as.numeric(value)
        }
        parsed[[param]] = value
    }
    if (print) {
        for (name in names(parsed)) {
            cat(name, '=', parsed[[name]], '\n')
        }
        cat('\n')
    }
    return(parsed)
}
