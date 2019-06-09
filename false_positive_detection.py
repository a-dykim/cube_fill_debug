"""
Find the false positive cases 
"""

import numpy as np
import pickle 
from astropy.io import fits
from cube_fil_finder.util import moments
from astropy import wcs
import pickle 


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
    print len(ref_set)
    overlap_frac = float(len(overlap_set)) / float(len(ref_set))
    return overlap_frac
    


def tree_candidate(tree, all_nodes):
    """
    Based on the tree's velocity & position, search through node dictionary to find possible tree candidate at first and last nodes
    
    Returns to candidate diction
    """
    candidate_ovlp = {}
    
    first = tree.getNode(0) ### first node
    firstv = int(first.v_slice_index[0])
    first_candidate = filter_wcorners(all_nodes[firstv-1], first.corners_original)

    for k in first_candidate:
        test = first_candidate[k]
        testovlp = compute_overlap(test, first)
        candidate_ovlp[k] = [firstv-1,testovlp]
    
    last = tree.getNode(-1) ### last node
    lastv = int(last.v_slice_index[0])
    last_candidate = filter_wcorners(all_nodes[lastv+1], last.corners_original)
    
    for k in last_candidate:
        test = last_candidate[k]
        testovlp = compute_overlap(test, last)
        candidate_ovlp[k] = [lastv+1,testovlp]
    
    return candidate_ovlp
    

def false_positive(tree, all_nodes):
	"""
	Find the false positive cases
	"""
	false_positive = {}
	for k in tree:
	    one_tree = tree[k]
	    candidates = tree_candidate(one_tree, all_nodes)
	    for n in candidates:
	        if candidates[n][1] >= 0.75:
	            false_positive[n] = [k, candidates[n][0], candidates[n][1]]

	return false_positive


































    

    