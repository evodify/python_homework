#!/usr/bin/env python2

from ete2 import Tree

############################# options #############################



############################# functions #############################

def assign_node_names(tree):
	'''Assignes names to all internal nodes without names.
	The string "inter_" + number are used for names.'''
	myNodeName = 1
	for node in tree.traverse():
		if not node.name:
			node.add_features(name='inter_'+str(myNodeName))
		  	myNodeName += 1

def count_subtree_internal_nodes(tree, group):
	''' Counts intenal nodes from leaves of group
	to the most recent common ancestor (MRCA) of the group'''
	ancestor = tree.get_common_ancestor(group).name # get MRCA name
	nodes = [ancestor] # make if the fist node
	for group_name in group:
		leave = tree.search_nodes(name=group_name)[0] # extract a leave from a tree
		for a in leave.iter_ancestors():
			if a.name != ancestor:
				nodes.append(a.name)
			else:
				break # break if the MRCA is reached
	print set(nodes) # for debuging
	return int(len(set(nodes))) # count unique node names

def count_all_internal_nodes(tree):
	''' Counts all intenal nodes in a tree'''
	nodeL = 0
	for node in tree.traverse():
	 	if not node.is_leaf():
	 		nodeL +=1
	return nodeL

def gsi(tree, group):
	'''
	GS = n/sum(du-2)
	where d is the degree of node u of U total nodes uniting a group
	(estimated coalescent events) through a most recent common an
	cestor
	n is the minimum number of nodes (coalescent events)
	required to unite a group of size n + 1 through a most recent
	GSI = (obsGS - minGS)/(maxGS-minGS)
	'''
	maxGS = 1.0
	# count internal nodes
	dGr = count_subtree_internal_nodes(tree, group)
	n = len(group) - 1
	GS = float(n)/float(dGr)
	minGS = float(n)/float(count_all_internal_nodes(tree))
	GSI = (GS - minGS) / (maxGS - minGS)
	return GSI

############################# analysis #############################


t = Tree("tree.nwk")
# print t.get_ascii(show_internal=True) # for debuging

gr1 = ["a1", "a2", "a3", "a4"]
gr2 = ["b1", "b2", "b3", "b4"]

# add names to internal nodes for convinicence of codding and debuging.
assign_node_names(t)
print t.get_ascii(show_internal=True) # for debuging 

# amono = t.check_monophyly(values=["a1", "a2", "a3", "a4"], target_attr="name")
# bmono = t.check_monophyly(values=["b1", "b2", "b3", "b4"], target_attr="name")
# print bmono

# # count internal nodes
# dGr1 = count_subtree_internal_nodes(t, gr1) # for debuging 
# dGr2 = count_subtree_internal_nodes(t, gr2) # for debuging 
# # print "Number of internal nodes = ", dGr1, "&", dGr2 # for debuging 

# calculate GSI
print
print 'GSI of a:\n', gsi(t, gr1)
print
print 'GSI of b:\n', gsi(t, gr2)

