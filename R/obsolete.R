# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

get_tip_labels_obsolete = function(phy, node_num, out=NULL, mode='repetitive') {
    if (mode=='recursive') {
        if (is.null(out)) {
            out = list()
            out[['nn']] = c(node_num)
            out[['label']] = c()
            node_num = NULL
        }
        if (length(out[['nn']])==0) {
            return(out[['label']])
        }
        if (out[['nn']][1] <= length(phy$tip.label)) {
            out[['label']] = c(out[['label']], phy$tip.label[out[['nn']][1]])
        } else {
            out[['nn']] = c(out[['nn']], get_children_num(phy, out[['nn']][1]))
        }
        out[['nn']] = out[['nn']][out[['nn']]!=out[['nn']][1]]
        Recall(phy=phy, node_num=NULL, out=out)
    } else if (mode=='repetitive') {
        node_nums = as.integer(node_num)
        leaf_nums = 1:length(phy$tip.label)
        while (! all(node_nums %in% leaf_nums)) {
            children_nums = c()
            for (nn in node_nums) {
                if (nn > length(phy$tip.label)) {
                    cnums = get_children_num(phy, nn)
                    children_nums = c(children_nums, cnums)
                    node_nums = node_nums[node_nums!=nn]
                }
            }
            node_nums = c(node_nums, children_nums)
        }
        tip_labels = get_node_name_by_num(phy, node_nums)
        return(tip_labels)
    }
}