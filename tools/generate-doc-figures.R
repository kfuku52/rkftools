# Run from the repository root, normally with make figures. The published
# examples are the source of truth for all tree, table, and score values.
source("tools/load-local.R")
invisible(load_local_package())

run_usage_examples <- function(path) {
    examples <- new.env(parent=globalenv())
    inside <- FALSE
    code <- character()
    for (line in readLines(path, warn=FALSE)) {
        if (!inside && identical(trimws(line), "```r")) {
            inside <- TRUE
            code <- character()
        } else if (inside && identical(trimws(line), "```")) {
            eval(parse(text=code), envir=examples)
            inside <- FALSE
        } else if (inside) {
            code <- c(code, line)
        }
    }
    if (inside) stop("Unclosed R code block in ", path)
    examples
}

examples <- run_usage_examples("docs/usage.md")
ink <- "#23313F"
muted <- "#5F7080"
border <- "#CEDAE5"
blue <- "#4779A8"
green <- "#14866D"
purple <- "#7957A5"

render_png <- function(path, width, height, draw) {
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    grDevices::png(path, width=width, height=height, res=160, pointsize=14,
        bg="white")
    on.exit(grDevices::dev.off())
    draw()
    message("Generated ", path)
}

label <- function(text, x, y, size=13, color=ink, bold=FALSE, mono=FALSE,
                  just="left") {
    grid::grid.text(text, x, y, just=just, gp=grid::gpar(col=color,
        fontsize=size, fontface=if (bold) "bold" else "plain",
        fontfamily=if (mono) "mono" else "sans"))
}

panel <- function(x, y, width, height, fill="white") {
    grid::grid.roundrect(x, y, width, height, just=c("left", "bottom"),
        r=grid::unit(0.08, "inches"), gp=grid::gpar(col=border, fill=fill))
}

draw_table <- function(data, x, y, width, height, widths=rep(1, ncol(data))) {
    values <- rbind(names(data), vapply(data, function(column) {
        format(column, trim=TRUE, digits=6)
    }, character(nrow(data))))
    widths <- width * widths / sum(widths)
    left <- x + c(0, head(cumsum(widths), -1L))
    row_height <- height / nrow(values)
    for (row in seq_len(nrow(values))) for (column in seq_len(ncol(values))) {
        middle_x <- left[[column]] + widths[[column]] / 2
        middle_y <- y + height - (row - 0.5) * row_height
        fill <- if (row == 1L) "#E8EFF5" else if (row %% 2L) "#F6F9FB" else "white"
        grid::grid.rect(middle_x, middle_y, widths[[column]], row_height,
            gp=grid::gpar(col=border, fill=fill, lwd=0.5))
        label(values[row,column], middle_x, middle_y, size=11,
            bold=row == 1L, just="center")
    }
}

conversion_arrow <- function(x, y, text, color, label_x, label_y) {
    grid::grid.lines(x, y, arrow=grid::arrow(length=grid::unit(0.1, "inches"),
        type="closed"), gp=grid::gpar(col=color, fill=color, lwd=2))
    label(text, label_x, label_y, size=11, color=color, bold=TRUE, just="center")
}

render_png("man/figures/table2phylo_roundtrip.png", 2400, 1500, function() {
    tree <- examples$tree
    grid::grid.newpage()
    label("One tree, three representations", 0.5, 0.955, size=23, bold=TRUE,
        just="center")
    label("Computed from the conversion example in docs/usage.md", 0.5, 0.915,
        color=muted, just="center")

    panel(0.06, 0.79, 0.88, 0.10, fill="#F6F9FB")
    label("Newick text", 0.08, 0.862, size=15, bold=TRUE)
    label(examples$newick, 0.08, 0.82, size=12, mono=TRUE)
    conversion_arrow(c(0.23, 0.23), c(0.785, 0.677), "ape::read.tree()",
        green, 0.135, 0.735)
    conversion_arrow(c(0.40, 0.40), c(0.677, 0.785), "ape::write.tree()",
        purple, 0.495, 0.735)

    panel(0.04, 0.10, 0.43, 0.57)
    label("ape::phylo object", 0.06, 0.628, size=18, bold=TRUE)
    label("$edge and $edge.length", 0.06, 0.578, mono=TRUE)
    edges <- data.frame(parent=tree$edge[,1], child=tree$edge[,2],
        edge.length=tree$edge.length)
    draw_table(edges, 0.06, 0.30, 0.39, 0.25, widths=c(1, 1, 1.5))
    label(paste0("$tip.label: ", paste(tree$tip.label, collapse=", ")),
        0.06, 0.26, size=12, mono=TRUE)
    label(paste0("$node.label: ", paste(tree$node.label, collapse=", ")),
        0.06, 0.22, size=11, mono=TRUE)
    label(paste0("$Nnode: ", tree$Nnode), 0.06, 0.178, size=12, mono=TRUE)
    label(paste0("$root.edge: ", if (is.null(tree$root.edge)) "NULL" else tree$root.edge),
        0.06, 0.138, size=12, mono=TRUE)

    panel(0.60, 0.10, 0.36, 0.57)
    label("Branch table (round trip)", 0.618, 0.628, size=18, bold=TRUE)
    draw_table(examples$roundtrip_table, 0.618, 0.25, 0.324, 0.33,
        widths=c(1.1, 1, 1, 1.7, 0.8))
    label("parent and sister refer to branch_id values", 0.618, 0.20,
        size=10, color=muted)
    label("node_name stores the tip or internal-node label", 0.618, 0.165,
        size=10, color=muted)
    conversion_arrow(c(0.595, 0.475), c(0.44, 0.44), "table2phylo()",
        green, 0.535, 0.47)
    conversion_arrow(c(0.475, 0.595), c(0.32, 0.32), "phylo2table()",
        purple, 0.535, 0.29)
    label("Newick and branch tables use ape::phylo as their bridge.", 0.5, 0.05,
        size=11, color=muted, just="center")
})

render_png("man/figures/root_position_species_overlap.png", 2400, 1400, function() {
    plot_root_position_species_overlap(examples$gene_tree)
})
