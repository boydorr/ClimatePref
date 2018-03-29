function getinternalnodes(t::AbstractTree)
    return collect(nodenamefilter(x->!isleaf(x)& !isroot(x), t))
end

function drop_tip!(t::AbstractTree, tip::Vector{String})
    tree_names = getleafnames(t)
    cut_names = setdiff(tree_names, tip)
    for i in cut_names[2:end]
        deletenode!(t, i)
    end
    while length(setdiff(collect(nodenamefilter(isleaf, t)), tip)) > 0
        nodes = setdiff(collect(nodenamefilter(isleaf, t)), tip)
        map(x -> deletenode!(t, x), nodes)
    end
    while sum(map(x-> length(getchildren(t, x)).< 2, getinternalnodes(t))) > 0
        inner_nodes = getinternalnodes(t)
        remove_nodes = find(map(x-> length(getchildren(t, x)).< 2, inner_nodes))
        for i in remove_nodes
            parent = getparent(t, inner_nodes[i])
            parentbranch = getinbound(getnode(t, inner_nodes[i]))

            child = getchildren(t, inner_nodes[i])[1]
            childbranch = getoutbounds(getnode(t, inner_nodes[i]))[1]

            len = distance(t, parent, child)

function drop_tip(t::AbstractTree, tip::Vector{String}, trim.internal = TRUE, subtree = FALSE, root.edge = 0,
    rooted = is.rooted(phy), collapse.singles = TRUE, interactive = FALSE)
    Ntip <- length(getleafnames(t))
    tree_names = getleafnames(t)
    tip <- map(x -> find(x .== tree_names)[1], tip)
    out.of.range <- tip > Ntip
    if (any(out.of.range)) {
        warning("some tip numbers were larger than the number of tips: they were ignored")
        tip <- tip[!out.of.range]
    }

    wbl <- !is.null(phy$edge.length)
    if (length(tip) == Ntip - 1 && trim.internal) {
        i <- which(phy$edge[, 2] == (1:Ntip)[-tip])
        res <- list(edge = matrix(2:1, 1, 2), tip.label = phy$tip.label[phy$edge[i,
            2]], Nnode = 1L)
        class(res) <- "phylo"
        if (wbl)
            res$edge.length <- phy$edge.length[i]
        if (!is.null(phy$node.label))
            res$node.label <- phy$node.label[phy$edge[i, 1] -
                Ntip]
        return(res)
    }
    if (!rooted && subtree) {
        phy <- root(phy, (1:Ntip)[-tip][1])
        root.edge <- 0
    }
    phy <- reorder(phy)
    NEWROOT <- ROOT <- Ntip + 1
    Nnode <- phy$Nnode
    Nedge <- dim(phy$edge)[1]
    if (subtree) {
        trim.internal <- TRUE
        tr <- reorder(phy, "postorder")
        N <- .C(node_depth, as.integer(Ntip), as.integer(tr$edge[,
            1]), as.integer(tr$edge[, 2]), as.integer(Nedge),
            double(Ntip + Nnode), 1L)[[5]]
    }
    edge1 <- phy$edge[, 1]
    edge2 <- phy$edge[, 2]
    keep <- !logical(Nedge)
    keep[match(tip, edge2)] <- FALSE
    if (trim.internal) {
        ints <- edge2 > Ntip
        repeat {
            sel <- !(edge2 %in% edge1[keep]) & ints & keep
            if (!sum(sel))
                break
            keep[sel] <- FALSE
        }
        if (subtree) {
            subt <- edge1 %in% edge1[keep] & edge1 %in% edge1[!keep]
            keep[subt] <- TRUE
        }
        if (root.edge && wbl) {
            degree <- tabulate(edge1[keep])
            if (degree[ROOT] == 1) {
                j <- integer(0)
                repeat {
                  i <- which(edge1 == NEWROOT & keep)
                  j <- c(i, j)
                  NEWROOT <- edge2[i]
                  degree <- tabulate(edge1[keep])
                  if (degree[NEWROOT] > 1)
                    break
                }
                keep[j] <- FALSE
                if (length(j) > root.edge)
                  j <- 1:root.edge
                NewRootEdge <- sum(phy$edge.length[j])
                if (length(j) < root.edge && !is.null(phy$root.edge))
                  NewRootEdge <- NewRootEdge + phy$root.edge
                phy$root.edge <- NewRootEdge
            }
        }
    }
    if (!root.edge)
        phy$root.edge <- NULL
    phy$edge <- phy$edge[keep, ]
    if (wbl)
        phy$edge.length <- phy$edge.length[keep]
    TERMS <- !(phy$edge[, 2] %in% phy$edge[, 1])
    oldNo.ofNewTips <- phy$edge[TERMS, 2]
    if (subtree) {
        i <- which(tip %in% oldNo.ofNewTips)
        if (length(i)) {
            phy$tip.label[tip[i]] <- "[1_tip]"
            tip <- tip[-i]
        }
    }
    n <- length(oldNo.ofNewTips)
    phy$edge[TERMS, 2] <- rank(phy$edge[TERMS, 2])
    phy$tip.label <- phy$tip.label[-tip]
    if (subtree || !trim.internal) {
        node2tip <- oldNo.ofNewTips[oldNo.ofNewTips > Ntip]
        new.tip.label <- if (subtree) {
            paste("[", N[node2tip], "_tips]", sep = "")
        }
        else {
            if (is.null(phy$node.label))
                rep("NA", length(node2tip))
            else phy$node.label[node2tip - Ntip]
        }
        phy$tip.label <- c(phy$tip.label, new.tip.label)
    }
    phy$Nnode <- dim(phy$edge)[1] - n + 1L
    newNb <- integer(Ntip + Nnode)
    newNb[NEWROOT] <- n + 1L
    sndcol <- phy$edge[, 2] > n
    newNb[sort(phy$edge[sndcol, 2])] <- (n + 2):(n + phy$Nnode)
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    storage.mode(phy$edge) <- "integer"
    if (!is.null(phy$node.label))
        phy$node.label <- phy$node.label[which(newNb > 0) - Ntip]
    if (collapse.singles)
        phy <- collapse.singles(phy)
    phy
}
