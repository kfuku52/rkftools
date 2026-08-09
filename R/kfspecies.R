# Title     : TODO
# Objective : TODO
# Created by: Codex
# Created on: 2026-03-21

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
