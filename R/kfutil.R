# Shared argument, validation, optional-package, and parallelism helpers.

.coerce_parsed_arg_value = function(value) {
    if (length(value) != 1) {
        return(value)
    }
    if (grepl('^[+-]?0[0-9]+$', value)) {
        return(value)
    }
    is_numeric_literal = grepl(
        '^[+-]?((0|[1-9][0-9]*)(\\.[0-9]*)?|\\.[0-9]+)([eE][+-]?[0-9]+)?$',
        value
    )
    if (!is_numeric_literal) {
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
    if (!startsWith(arg, '--')) {
        stop('Invalid long argument: expected a "--name" or "--name=value" argument, got "', arg, '".')
    }
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
    if (!is.character(value) || length(value) != 1 || is.na(value)) {
        stop(arg_name, ' must be a single non-missing string.')
    }
    if (!allow_empty && trimws(value) == '') {
        stop(arg_name, ' must be a single non-empty string.')
    }
    value
}

.normalize_single_logical_arg = function(value, arg_name) {
    if (!is.logical(value) || length(value) != 1 || is.na(value)) {
        stop(arg_name, ' must be a single non-missing logical value.')
    }
    value
}

.normalize_finite_numeric_scalar = function(
    value,
    arg_name,
    min_value=-Inf,
    max_value=Inf
) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
            !is.finite(value) || value < min_value || value > max_value) {
        range_text = if (is.finite(min_value) || is.finite(max_value)) {
            paste0(' in [', min_value, ', ', max_value, ']')
        } else {
            ''
        }
        stop(arg_name, ' must be a single finite numeric value', range_text, '.')
    }
    as.numeric(value)
}

.normalize_integerish = function(
    value,
    arg_name,
    allow_na=FALSE,
    allow_empty=FALSE,
    min_value=-.Machine$integer.max,
    max_value=.Machine$integer.max
) {
    if (!is.atomic(value) || is.factor(value)) {
        stop(arg_name, ' must contain integer values.')
    }
    if (!length(value)) {
        if (allow_empty) {
            return(integer(0))
        }
        stop(arg_name, ' must contain at least one integer value.')
    }
    numeric_value = suppressWarnings(as.numeric(value))
    conversion_failed = is.na(numeric_value) & !is.na(value)
    non_integer = !is.na(numeric_value) & (
        !is.finite(numeric_value) |
        numeric_value != trunc(numeric_value) |
        numeric_value < min_value |
        numeric_value > max_value
    )
    if (any(conversion_failed | non_integer)) {
        invalid = unique(as.character(value[conversion_failed | non_integer]))
        stop(
            arg_name,
            ' must contain integer values in [', min_value, ', ', max_value,
            ']. Invalid value(s): ', paste(invalid, collapse=', '), '.'
        )
    }
    if (!allow_na && anyNA(numeric_value)) {
        stop(arg_name, ' must not contain missing values.')
    }
    as.integer(numeric_value)
}

.normalize_choice_arg = function(value, arg_name, choices) {
    value = .normalize_single_string_arg(value, arg_name, allow_empty=FALSE)
    if (!(value %in% choices)) {
        stop(
            arg_name, ' must be one of: ',
            paste(sprintf('"%s"', choices), collapse=', '),
            '.'
        )
    }
    value
}

#' Parse command-line long arguments
#'
#' @param args Character values in `--name` or `--name=value` form.
#' @param print Whether to print parsed values. Credential-like values are
#'   redacted when printed.
#' @return A named list of parsed values.
#' @examples
#' get_parsed_args(c("--threads=2", "--dry-run"))
#' @export
get_parsed_args = function(args, print=FALSE) {
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
        if (parsed_item[['param']] %in% names(parsed)) {
            stop('Duplicate long argument: --', parsed_item[['param']])
        }
        parsed[[parsed_item[['param']]]] = parsed_item[['value']]
    }
    if (print) {
        for (name in names(parsed)) {
            display_value = parsed[[name]]
            if (grepl('(token|secret|password|passwd|credential|api[-_]?key|private[-_]?key)', name, ignore.case=TRUE)) {
                display_value = '<redacted>'
            }
            cat(name, '=', display_value, '\n')
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
    auto_cores = if (detected_cores < 2) 1L else min(detected_cores - 1L, 8L)

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

#' Test whether a value is blank
#'
#' @param x An object to inspect.
#' @param false.triggers Whether logical false values should count as blank.
#' @return A single logical value.
#' @export
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
