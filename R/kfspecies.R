# Shared species-label parsers.

.normalize_species_parser_arg = function(value, arg_name='species_parser') {
    value = .normalize_single_string_arg(
        value=value,
        arg_name=arg_name,
        allow_empty=FALSE
    )
    supported_parsers = c('legacy', 'taxonomic')
    if (!(value %in% supported_parsers)) {
        stop(
            arg_name,
            ' must be one of: ',
            paste(sprintf('"%s"', supported_parsers), collapse=', '),
            '.'
        )
    }
    value
}

.taxonomic_genus_only_placeholders = c('sp', 'spp')
.taxonomic_single_token_qualifiers = c('cf', 'aff', 'nr', 'x')
.taxonomic_paired_token_qualifiers = c(
    'subsp', 'ssp', 'subspecies',
    'var', 'variety',
    'forma', 'form', 'f',
    'strain', 'substrain',
    'serovar', 'serotype', 'serogroup',
    'pathovar', 'pv',
    'biovar', 'biotype', 'chemovar', 'morphovar',
    'cultivar', 'cv',
    'isolate',
    'group', 'subgroup', 'complex', 'clade', 'lineage',
    'section', 'series', 'ecotype', 'breed'
)

.taxonomic_species_parser_token_count = function(split_label) {
    species_len = 2L
    if (length(split_label) < species_len) {
        return(length(split_label))
    }

    if (tolower(split_label[2]) %in% .taxonomic_genus_only_placeholders) {
        if (length(split_label) >= 3L) {
            return(3L)
        }
        return(2L)
    }
    if (tolower(split_label[2]) %in% .taxonomic_single_token_qualifiers) {
        if (length(split_label) >= 3L) {
            return(3L)
        }
        return(2L)
    }

    while (length(split_label) > species_len) {
        next_token = split_label[species_len + 1L]
        next_token_norm = tolower(next_token)
        if (next_token_norm %in% .taxonomic_single_token_qualifiers) {
            species_len = species_len + 1L
            next
        }
        if (next_token_norm %in% .taxonomic_paired_token_qualifiers) {
            species_len = species_len + 1L
            if (length(split_label) > species_len) {
                species_len = species_len + 1L
            }
            next
        }
        break
    }
    species_len
}

.extract_species_tokens_from_label = function(label, species_parser='legacy', sep='_', require_gene=FALSE) {
    if (length(label) != 1L || is.na(label)) {
        return(list(ok=FALSE, species_tokens=character(0), has_gene=FALSE))
    }

    split_label = strsplit(as.character(label), sep, fixed=TRUE)[[1]]
    if (species_parser == 'legacy') {
        species_len = 2L
    } else {
        species_len = .taxonomic_species_parser_token_count(split_label)
    }
    if (length(split_label) < species_len || species_len < 2L) {
        return(list(ok=FALSE, species_tokens=character(0), has_gene=FALSE))
    }
    species_tokens = split_label[seq_len(species_len)]
    if (any(is.na(species_tokens) | trimws(species_tokens) == '')) {
        return(list(ok=FALSE, species_tokens=character(0), has_gene=FALSE))
    }

    has_gene = length(split_label) > species_len
    if (require_gene && !has_gene) {
        return(list(ok=FALSE, species_tokens=character(0), has_gene=FALSE))
    }

    list(
        ok=TRUE,
        species_tokens=species_tokens,
        has_gene=has_gene
    )
}

.parse_species_labels = function(
    labels,
    species_parser='legacy',
    sep='_',
    output_sep=' ',
    require_gene=FALSE,
    fallback_label=FALSE
) {
    species_labels = rep(NA_character_, length(labels))
    parsed_ok = rep(FALSE, length(labels))
    for (i in seq_along(labels)) {
        parsed = .extract_species_tokens_from_label(
            label=labels[[i]],
            species_parser=species_parser,
            sep=sep,
            require_gene=require_gene
        )
        parsed_ok[[i]] = parsed[['ok']]
        if (parsed[['ok']]) {
            species_labels[[i]] = paste(parsed[['species_tokens']], collapse=output_sep)
        } else if (fallback_label && length(labels[[i]]) == 1L && !is.na(labels[[i]])) {
            species_labels[[i]] = as.character(labels[[i]])
        }
    }
    list(species_labels=species_labels, parsed_ok=parsed_ok)
}

#' Parse species names from labels
#'
#' @param a One or more labels.
#' @param species_parser Species-label convention.
#' @param sep Literal input separator.
#' @return Species names separated by spaces.
#' @examples
#' get_species_name("Homo_sapiens_gene1")
#' get_species_name("Genus_cf_species_gene1", species_parser="taxonomic")
#' @export
get_species_name = function(a, species_parser='legacy', sep='_') {
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    parsed = .parse_species_labels(
        labels=a,
        species_parser=species_parser,
        sep=sep,
        output_sep=' ',
        require_gene=FALSE,
        fallback_label=TRUE
    )
    return(parsed[['species_labels']])
}


#' Parse species names from tree tips
#'
#' @param phy A `phylo` tree.
#' @param sep Literal input and output separator.
#' @param species_parser Species-label convention.
#' @return A character vector aligned with tree tips.
#' @export
get_species_names = function(phy, sep='_', species_parser='legacy') {
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    parsed = .parse_species_labels(
        labels=phy[['tip.label']],
        species_parser=species_parser,
        sep=sep,
        output_sep=sep,
        require_gene=FALSE,
        fallback_label=FALSE
    )
    species_names = parsed[['species_labels']]
    if (any(!parsed[['parsed_ok']])) {
        bad_labels = phy[['tip.label']][!parsed[['parsed_ok']]]
        warning(
            'Leaf name(s) could not be interpreted with species_parser="',
            species_parser, '": ', paste(bad_labels, collapse=', '),
            call.=FALSE
        )
    }
    return(species_names)
}


#' Convert gene-bearing leaf labels to species names
#'
#' @param leaf_names Gene-bearing leaf labels.
#' @param use_underbar Whether output species names retain underscores.
#' @param species_parser Species-label convention.
#' @param sep Literal input separator.
#' @return Parsed species names, with `NA` for malformed labels.
#' @export
leaf2species = function(leaf_names, use_underbar=FALSE, species_parser='legacy', sep='_') {
    use_underbar = .normalize_single_logical_arg(
        value=use_underbar,
        arg_name='use_underbar'
    )
    species_parser = .normalize_species_parser_arg(
        value=species_parser,
        arg_name='species_parser'
    )
    sep = .normalize_single_string_arg(
        value=sep,
        arg_name='sep',
        allow_empty=FALSE
    )
    parsed = .parse_species_labels(
        labels=leaf_names,
        species_parser=species_parser,
        sep=sep,
        output_sep=if (use_underbar) '_' else ' ',
        require_gene=TRUE,
        fallback_label=FALSE
    )
    species_names = parsed[['species_labels']]
    if (any(!parsed[['parsed_ok']])) {
        bad_labels = leaf_names[!parsed[['parsed_ok']]]
        warning(
            'Leaf name(s) could not be interpreted as species-bearing labels with species_parser="',
            species_parser, '": ', paste(bad_labels, collapse=', '),
            call.=FALSE
        )
    }
    return(species_names)
}
