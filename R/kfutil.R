# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/20/18

.coerce_parsed_arg_value = function(value) {
    if (length(value) != 1) {
        return(value)
    }
    numeric_value = suppressWarnings(as.numeric(value))
    if (!is.na(numeric_value)) {
        return(numeric_value)
    }
    return(value)
}

.parse_long_arg = function(arg) {
    if (length(arg) != 1 || is.na(arg)) {
        stop('Invalid long argument: argument must be a single non-missing string.')
    }
    arg = as.character(arg)
    normalized = sub("^--", "", arg)
    eq_pos = regexpr("=", normalized, fixed=TRUE)[1]
    if (eq_pos == -1) {
        if (normalized == "") {
            stop('Invalid long argument: parameter name is empty in "', arg, '".')
        }
        return(list(param=normalized, value=TRUE))
    }
    param = substr(normalized, 1, eq_pos - 1)
    if (param == "") {
        stop('Invalid long argument: parameter name is empty in "', arg, '".')
    }
    value = substr(normalized, eq_pos + 1, nchar(normalized))
    value = .coerce_parsed_arg_value(value)
    return(list(param=param, value=value))
}

.normalize_single_string_arg = function(value, arg_name, allow_empty=TRUE) {
    if (length(value) != 1 || is.na(value)) {
        stop(arg_name, ' must be a single non-missing string.')
    }
    value = as.character(value)
    if (!allow_empty && trimws(value) == '') {
        stop(arg_name, ' must be a single non-empty string.')
    }
    value
}

.normalize_single_logical_arg = function(value, arg_name) {
    if (length(value) != 1 || is.na(value)) {
        stop(arg_name, ' must be a single non-missing logical value.')
    }
    logical_value = suppressWarnings(as.logical(value))
    if (is.na(logical_value)) {
        stop(arg_name, ' must be a single non-missing logical value.')
    }
    logical_value
}

get_parsed_args = function(args, print=TRUE) {
    print = .normalize_single_logical_arg(
        value=print,
        arg_name='print'
    )
    parsed = list()
    if (length(args) == 0) {
        if (print) {
            cat('\n')
        }
        return(parsed)
    }
    for (arg in args) {
        parsed_item = .parse_long_arg(arg)
        parsed[[parsed_item[['param']]]] = parsed_item[['value']]
    }
    if (print) {
        for (name in names(parsed)) {
            cat(name, '=', parsed[[name]], '\n')
        }
        cat('\n')
    }
    return(parsed)
}

.get_optional_pkg_fun = function(pkg_name, fun_name) {
    if (!requireNamespace(pkg_name, quietly=TRUE)) {
        stop("'", pkg_name, "' package not found, please install it.")
    }
    if (!exists(fun_name, envir=asNamespace(pkg_name), mode="function", inherits=FALSE)) {
        stop("Function '", fun_name, "' was not found in package '", pkg_name, "'.")
    }
    get(fun_name, envir=asNamespace(pkg_name), mode="function")
}

.resolve_parallel_cores = function(requested=NULL, max_tasks=Inf, auto_when_missing=FALSE) {
    detected_cores = suppressWarnings(as.integer(parallel::detectCores(logical=FALSE)))
    if (is.na(detected_cores) || detected_cores < 1) {
        detected_cores = 1L
    }
    auto_cores = if (detected_cores < 2) 1L else detected_cores - 1L

    has_requested = !(is.null(requested) || length(requested) == 0)
    requested_cores = suppressWarnings(as.integer(requested))
    if (!length(requested_cores)) {
        requested_cores = NA_integer_
    } else {
        requested_cores = requested_cores[[1]]
    }
    if (!has_requested && auto_when_missing) {
        num_parallel = auto_cores
    } else if (!has_requested || is.na(requested_cores) || requested_cores < 1) {
        num_parallel = 1L
    } else {
        num_parallel = requested_cores
    }

    option_cap = suppressWarnings(as.integer(getOption("rkftools.max_cores", NA)))
    if (!length(option_cap)) {
        option_cap = NA_integer_
    } else {
        option_cap = option_cap[[1]]
    }
    if (!is.na(option_cap) && option_cap >= 1) {
        num_parallel = min(num_parallel, option_cap)
    }

    task_cap = suppressWarnings(as.integer(max_tasks))
    if (!length(task_cap)) {
        task_cap = NA_integer_
    } else {
        task_cap = task_cap[[1]]
    }
    if (!is.na(task_cap) && task_cap >= 1) {
        num_parallel = min(num_parallel, task_cap)
    }

    check_limit_cores = Sys.getenv("_R_CHECK_LIMIT_CORES_", unset="")
    is_check_limited = nzchar(check_limit_cores) &&
        !(tolower(check_limit_cores) %in% c("0", "false", "no"))
    if (is_check_limited) {
        num_parallel = min(num_parallel, 2L)
    }

    if (is.na(num_parallel) || num_parallel < 1) {
        num_parallel = 1L
    }
    as.integer(num_parallel)
}

is.blank = function(x, false.triggers=FALSE){
    # https://stackoverflow.com/questions/19655579/a-function-that-returns-true-on-na-null-nan-in-r
    false.triggers = .normalize_single_logical_arg(
        value=false.triggers,
        arg_name='false.triggers'
    )
    if(is.function(x)) return(FALSE) # Some of the tests below trigger
                                     # warnings when used on functions
    if (is.null(x) || length(x) == 0) {
        return(TRUE)
    }
    if (all(is.na(x))) {
        return(TRUE)
    }
    if (is.character(x)) {
        x_trim = trimws(x)
        is_empty_or_na = is.na(x_trim) | (x_trim == "")
        if (all(is_empty_or_na)) {
            return(TRUE)
        }
    }
    if (false.triggers) {
        logical_x = suppressWarnings(as.logical(x))
        if (!any(is.na(logical_x)) && all(!logical_x)) {
            return(TRUE)
        }
    }
    return(FALSE)
}
