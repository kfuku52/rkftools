# Title     : TODO
# Objective : TODO
# Created by: kf
# Created on: 5/19/18

get_tip_labels_obsolete = function(phy, node_num, out=NULL, mode='repetitive') {
    if (mode=='recursive') {
        state = out
        start_node_num = node_num
        if (is.null(state)) {
            state = list()
            state[['nn']] = c(start_node_num)
            state[['label']] = c()
        }
        if (length(state[['nn']])==0) {
            return(state[['label']])
        }
        current_node_num = state[['nn']][1]
        if (current_node_num <= length(phy[['tip.label']])) {
            state[['label']] = c(state[['label']], phy[['tip.label']][current_node_num])
        } else {
            state[['nn']] = c(state[['nn']], get_children_num(phy, current_node_num))
        }
        state[['nn']] = state[['nn']][state[['nn']]!=current_node_num]
        Recall(phy=phy, node_num=NULL, out=state)
    } else if (mode=='repetitive') {
        node_nums = as.integer(node_num)
        leaf_nums = seq_along(phy[['tip.label']])
        while (! all(node_nums %in% leaf_nums)) {
            children_nums = c()
            for (nn in node_nums) {
                if (nn > length(phy[['tip.label']])) {
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
