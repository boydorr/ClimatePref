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

            deletebranch!(t, parentbranch)
            deletebranch!(t, childbranch)
            delete!(getnodes(t), inner_nodes[i])
            delete!(t.noderecords, inner_nodes[i])

            addbranch!(t, parent, child, len)
        end
    end
end
