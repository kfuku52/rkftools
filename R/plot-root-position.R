#' Example gene tree for species-overlap scoring
#'
#' A synthetic rooted gene tree with two gene copies from each of four
#' vertebrate species. Tip labels follow the `Genus_species_gene` convention.
#'
#' @format An object of class `phylo` with eight tips and branch lengths.
"species_overlap_example_tree"

#' Plot species-overlap scores for candidate roots
#'
#' Draws the tree topology alongside a bar chart of candidate-root scores on
#' the active graphics device. Edge numbers refer to rows of `phy$edge` in
#' both panels. All candidates with the minimum score are highlighted.
#' The plot uses the full device and sets a two-panel layout.
#'
#' @inheritParams get_root_position_dependent_species_overlap_scores
#' @return Invisibly, the numeric vector of scores in the input edge order.
#' @examples
#' plot_root_position_species_overlap(ape::unroot(species_overlap_example_tree))
#' @export
plot_root_position_species_overlap = function(phy, species_parser='legacy', sep='_') {
    scores = get_root_position_dependent_species_overlap_scores(
        phy, species_parser=species_parser, sep=sep
    )
    minimum = min(scores)
    best = scores == minimum
    ink = '#23313F'
    muted = '#5F7080'
    blue = '#4779A8'
    green = '#14866D'
    highlight = ifelse(best, green, ink)
    rooted = ape::is.rooted(phy)

    previous_par = graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(previous_par))
    graphics::layout(matrix(1:2, nrow=1), widths=c(1.35, 1))
    graphics::par(oma=c(2, 0, 4, 0), mar=c(3, 1, 4, 1), col=ink,
        col.axis=ink, col.lab=ink, col.main=ink, family='sans')
    ape::plot.phylo(phy, type=if (rooted) 'phylogram' else 'unrooted',
        use.edge.length=FALSE, lab4ut='horizontal', cex=0.85,
        edge.color=highlight, edge.width=ifelse(best, 5, 2),
        tip.color=ink, font=3, underscore=TRUE)

    # ape may reorder edges for plotting. The child node uniquely identifies
    # each edge, so translate the plotted rows back to the input edge order.
    plotted = get('last_plot.phylo', envir=ape::.PlotPhyloEnv)
    edge_index = match(plotted$edge[,2], phy$edge[,2])
    ape::edgelabels(text=edge_index, frame='circle', cex=0.75, font=2,
        bg=ifelse(best[edge_index], green, 'white'),
        col=ifelse(best[edge_index], 'white', ink))
    graphics::title(if (rooted) 'Rooted candidate edges' else 'Unrooted candidate edges',
        line=2, cex.main=1.1)
    graphics::mtext('Edge indices in phy$edge, topology only', line=0.6,
        cex=0.75, col=muted)
    graphics::legend('bottom', legend='Minimum-score candidate', pch=21,
        pt.bg=green, col=green, bty='n', cex=0.8, inset=-0.08, xpd=TRUE)

    graphics::par(mar=c(5, 4, 4, 2))
    label_padding = 0.05 * max(1, scores)
    positions = graphics::barplot(scores, names.arg=seq_along(scores),
        col=ifelse(best, green, blue), border=NA,
        ylim=c(0, max(scores) + 4 * label_padding), las=1, cex.names=0.85,
        xlab='Edge index in phy$edge', ylab='Species-overlap score')
    graphics::abline(h=minimum, lty=2, col=green, lwd=1.5)
    graphics::text(positions, scores + label_padding, labels=scores,
        cex=0.85, col=highlight)
    graphics::title('Root-position scores', line=2, cex.main=1.1)
    graphics::mtext('Lower scores indicate less species overlap', line=0.6,
        cex=0.75, col=muted)
    graphics::mtext('Species overlap across candidate roots', outer=TRUE,
        line=2, cex=1.5, font=2, col=ink)
    graphics::mtext(sprintf('Minimum score %s across %s of %s candidate edges.',
        minimum, sum(best), length(scores)), side=1, outer=TRUE,
        line=0.5, cex=0.85, col=muted)

    invisible(scores)
}
