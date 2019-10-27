import numpy as np


def order_corner(node):
    corn = node.corners_original
    return [corn[0][1], corn[1][1], corn[0][0], corn[1][0]]

def filter_wcorners(node, corners):
    """
    Search node based on position 
    """
    filtered = {}
    for k in node:
        one_branch = node[k]
        branch_corner = one_branch.corners_original
        if (corners[1][0] < branch_corner[0][0]) or (corners[0][0] > branch_corner[1][0]) or (corners[1][1] < branch_corner[0][1]) or (corners[0][1] > branch_corner[1][1]):
            pass
        else:
            filtered[k] = one_branch
            
    return filtered

def filter_tree_wcorners(trees, ref_tree):
    """
    Search tree based on position
    """
    corners = ref_tree.root_node.corners_original
    filtered = {}
    for k in trees:
        one_tree = trees[k]
        one_corner = one_tree.root_node.corners_original
        if (corners[1][0] < one_corner[0][0]) or (corners[0][0] > one_corner[1][0]) or (corners[1][1] < one_corner[0][1]) or (corners[0][1] > one_corner[1][1]):
            pass
        else:
            filtered[k] = one_tree
    return filtered


def config_pixeloc(node):
    """
    Config "valid" pixel locations of a given node mask with corner info
    """
    yinx, xinx = np.where(node.mask[::-1,:]) # valid indices and flip to origin='lower'
    corners = node.corners_original
    #print (corners)
    xmin = corners[0][1]
    ymax = corners[1][0]
    #print (xmin, ymax)
    xloc = xinx + xmin
    yloc = ymax - yinx
    #pixloc = np.vstack((xloc, yloc))
    return xloc, yloc



def compute_overlap(test, ref):
    """
    Compute the test's overlapping fraction based on ref
    """
    xref, yref = config_pixeloc(ref)
    ref_coord = np.array(zip(xref, yref))
    ref_set = set(map(tuple,ref_coord))
    
    xtest, ytest = config_pixeloc(test)
    test_coord = np.array(zip(xtest, ytest))
    test_set = set(map(tuple,test_coord))
    
    overlap_set = test_set.intersection(ref_set)
    #print len(ref_set)
    overlap_frac = float(len(overlap_set)) / float(len(ref_set))
    return overlap_frac

def compute_matchability(ref_tree, new_tree):
    """
    Check overlapping fraction of found tree
    """
    oldv = ref_tree.getTreeVelocityRange()
    newv = new_tree.getTreeVelocityRange()
    overlpv = list(set(oldv).intersection(set(newv)))
    
    if len(oldv) > len(newv):
        matchability = np.zeros(len(oldv))
        pvs = oldv
        print ("old one is longer?")
    else:
        matchability = np.zeros(len(newv))
        pvs = newv
        
    for i in overlpv:
        new_ind = np.where(newv==i)[0][0]
        old_ind = np.where(oldv==i)[0][0]
        new_node = new_tree.getNode(new_ind)
        old_node = ref_tree.getNode(old_ind)
        overlap = compute_overlap(new_node, old_node)
        matchability[new_ind] = overlap
    return pvs, matchability

def possible_merge(trees, ovlp_thresh=0.85):
    """
    Go through trees and find other trees that are overlaped with velocity and position
    """
    merg_tree = {}
    for k in trees:
        one_tree = trees[k]
        onevs = set(one_tree.getTreeVelocityRange())
        xoverlap = filter_tree_wcorners(trees, one_tree)
        candidates = []
        for j in xoverlap:
            two_tree = xoverlap[j]
            ovlp_frac = compute_overlap(two_tree.root_node, one_tree.root_node)
            twovs = set(two_tree.getTreeVelocityRange())
            inter = onevs.intersection(twovs)
            if (len(inter) > 0) and (ovlp_frac > ovlp_thresh) and (k is not j):
                candidates.append(j)
            
        if len(candidates) > 0:
            merg_tree[k] = candidates
    
    return merg_tree



## merge_flag = possible_merge(data)
