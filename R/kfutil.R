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

is.blank = function(x, false.triggers=FALSE){
    # https://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
    if(is.function(x)) return(FALSE) # Some of the tests below trigger
                                     # warnings when used on functions
    return(
        is.null(x) ||                # Actually this line is unnecessary since
        length(x) == 0 ||            # length(NULL) = 0, but I like to be clear
        all(is.na(x)) ||
        all(x=="") ||
        (false.triggers && all(!x))
    )
}