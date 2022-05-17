#!/usr/bin/env python
#Change log
#2.6 Added attributes to edges July 2016 
#2.6.2 Added argument to take the solution file in .csv
#Error log
import sys
import os
import argparse
import csv #for reading a comma seperated value (CSV) table 
from xml.dom import minidom #for parsing XML Document Object Models 
import semanticnet as sn #this is the library for creating the JSON graphs 
from math import log10 

parser = argparse.ArgumentParser()
parser.add_argument("VCML_File", help="Enter a .VCML file to convert")
parser.add_argument("Solution_file", help="Enter a .csv file with time dependent solutions")
parser.add_argument("outfile",help="Enter an output file name with .json extension")
args = parser.parse_args()

#g=sn.Graph() Changed to DiGraph which I think stands for Directed Graph!!
g=sn.DiGraph()

xmldoc=minidom.parse(args.VCML_File) 

mathmodel1=xmldoc.getElementsByTagName("Model")[0]

#Get a list of all the compounds 
localCompounds= mathmodel1.getElementsByTagName("LocalizedCompound")

for compound in localCompounds:
    compoundName = compound.getAttribute("Name")
    compoundID = compound.getAttribute("KeyValue")
    compoundNode = g.add_node({"label":compoundName, "key val":compoundID,"type":"Compound"})
    g.set_node_attribute(compoundNode, "og:space:color", [0.640, 0.313, 0.746, 1.0])
    g.set_node_attribute(compoundNode, "og:space:icon", "shapes/hexagon")
    g.set_node_attribute(compoundNode, "og:space:activity", "0")
    print("Read compound: ", compoundName)

#Now to get reactions which will be reaction nodes. Then we will link reaction and compound nodes.  It might also be good to look at reactions as nodes and species as edges  
reactions = mathmodel1.getElementsByTagName("SimpleReaction")
for reaction in reactions:
    rname = reaction.getAttribute("Name")
    rKeyVal= reaction.getAttribute("KeyValue")
    
    # save the list of reaction nodes locally for making edges 
    reactionNode = g.add_node({"label":rname, "key val":rKeyVal,"type":"Reaction"})
    g.set_node_attribute(reactionNode, "og:space:color", [1.0, 0.0, 0.0, 1.0])
    g.set_node_attribute(reactionNode, "og:space:icon", "shapes/star")
    g.set_node_attribute(reactionNode, "og:space:activity", "0")
    print ("Read reaction: ",rname)

# Lastly load up 'Modifiers' which I think are just enzymes  
#modifiers = mathmodel1.getElementsByTagName("Modifier")
#for modifier in modifiers:
#    modname = modifier.getAttribute("Name"

#Add edges: Go through reactions to add reactants and products as edges
#Run between reaction node and corrosponding compound nodes   
g.cache_nodes_by("label")
for reaction in reactions:
    reactants = reaction.getElementsByTagName("Reactant")  
    #for each reaction, run through the compounds in that reaction 
    reactionName = reaction.getAttribute("Name")

    
    reactionNode= g.get_nodes_by_attr("label",reactionName)
    #reactionNode ends up as a one element list where the element is a dict data type
    
    reactionNode_ID=reactionNode[0]['id']
    for compound in reactants:
        compoundName=compound.getAttribute("LocalizedCompoundRef")  
        compoundNode=g.get_nodes_by_attr("label",compoundName)
        compoundNode_ID=compoundNode[0]['id'] 
        print ("In reaction: {} added reactant: {}".format(reactionName, compoundName))
        #assuming a net forward reaction use the edge command below.  If it is going the otherway switch reaction and compound ID's 
        reactantEdge=g.add_edge(compoundNode_ID, reactionNode_ID)
	    #Set reactant edge properties 
        g.set_edge_attribute(reactantEdge, "space:activity", 0.4) #Make this proportional to concentration of species 
        g.set_edge_attribute(reactantEdge, "space:width", 0.1)

    #Now add connetions between reaction nodes and product nodes 
    products = reaction.getElementsByTagName("Product")
    for product in products:
        productName = product.getAttribute("LocalizedCompoundRef")
        productNode = g.get_nodes_by_attr("label",productName)
        productNode_ID = productNode[0]['id'] 
        productEdge=g.add_edge(reactionNode_ID,productNode_ID)
        print ("In reaction: {} added product: {}".format(reactionName,productName))
        #Set reactant edge properties 
        g.set_edge_attribute(productEdge, "space:activity", 0.4) #Make this proportional to concentration of species 
        g.set_edge_attribute(productEdge, "space:width", 0.1)

    #Connect Modifiers (enzymes) to reactions 
    modifiers = reaction.getElementsByTagName("Modifier")
    print ("modifiers are:", modifiers)
    for modifier in modifiers:
        modifierName = modifier.getAttribute("LocalizedCompoundRef")
        modifierNode = g.get_nodes_by_attr("label",modifierName)
        modifierNode_ID = modifierNode[0]['id'] 
        modifierEdge=g.add_edge(modifierNode_ID,reactionNode_ID)
        print ("In reaction: {} added modifier: {}".format(reactionName,modifierName))

	    #Set modifier edge properties 
        g.set_edge_attribute(modifierEdge, "space:icon", "styles/dashed")
        g.set_edge_attribute(modifierEdge, "space:activity", 0.4) #Make this proportional to concentration of species 
        g.set_edge_attribute(modifierEdge, "space:width", 0.1)
 
#Read in a solution file from a CSV datafile 
solutionFile=open(args.Solution_file)
solfileReader=csv.reader(solutionFile)
#1. Read in solution file: Columns are time, reaction rates and species 
#2. Find the max concentration of each species and normalize 
#3 For each time write the current concentrations for each species as space:activity, node:activity or edge:width 
solData=list(solfileReader)  #Read in the solution data using a file reader type, ends up as a list of strings
floatData = []  #floatData is a list 
specNames=solData[0][1:]  #Get the first row excluding the "time" heading
for row in solData[1:]: #Skip the first row (list) with the headings  
    floatData.append([float (i) for i in row[1:]]) 

#Fix problem with names in the solution file not matching the VCML 
newspecNames=[] #Create an empty list
print("The species names are:",specNames)
for specName in specNames:
    spaceLocation = specName.find(' ')
    if spaceLocation!=-1: 
        specName=specName[(spaceLocation+1):] 
    newspecNames.append(specName)
    print("the current name is:",specName)
specNames=newspecNames 
print("The species names after the loop are:",specNames)

#Find the max value in each column and normalize the values 
maxVals=[] 
minVals=[]
normVals=[]
colNo=0  
for col in specNames: 
     maxVals.append(max(l[colNo] for l in floatData))
     minVals.append(min(l[colNo] for l in floatData))
     # Deal with the case where all the data in a column are zeros
     if maxVals[colNo]==0:                          
         normVals.append(list(1 for l in floatData))  #Set all values to 1 so we can see the edges
     else: #Normalize the column data 
         if maxVals[colNo]!=minVals[colNo]:
             normVals.append(list((l[colNo]-minVals[colNo])/maxVals[colNo] for l in floatData)) 
         else:  #normalize without subtracting minVals which just makes everything in the column to 1
             normVals.append(list(1 for l in floatData))         
     colNo+=1 
print("The Max vals are: ",maxVals)
print("the number of columns is:",colNo)



specNo=int(0) 
timeline_counter=0 
'''
#Now add property nodes to json file to reflect time change
#Nodes should still be sorted by label 

specNames = [w.replace('J_', '') for w in specNames]  #Need this to fix a little issue with VCell where reaction nodes have a different name in the solution file than in the vcml 
for nodes in specNames:
    graphNode= g.get_nodes_by_attr("label",nodes)
    print(nodes,graphNode) 
    graphNode_ID=graphNode[0]['id']
    print(type(graphNode_ID)) 
    for concentration in normVals[0:][specNo]: #Which column we are using here deosnt really matter but skip the time
	 g.add_event(timeline_counter, "graph:set_node_attribute",
         {
            'id': graphNode_ID,
            'name': 'og:space:activity',
            'type': 'float',
            'value': str(concentration)
         } 
         )
         print('Node:', nodes,'Concentration:',concentration, 'specNo:',specNo)
         timeline_counter+=100
    timeline_counter=0
    specNo+=1
'''
#Change width of edges from species to indicate concentration and activity of reaction nodes to indicate flux 

#set all the concentrations to '1' at the first time point so we can see all of the edges in the graph
#Otherwise the first view of the graph shows the model's initial conditions 
for i in range(colNo): 
     normVals[i][0]= 1 
     print("i is: ", i)

print("normVals is ", normVals)

g.cache_edges_by("src")
print(specNames) 
for nodes in specNames: 
    if 'J_' in nodes: 	#Then this is a reaction node. Can do this as all reactions are labeled J_r# in solution files 
        nodes=nodes.replace('J_','') #Need this to fix a little issue with VCell where reaction nodes have a different name in the solution file than in the vcml 
        graphNode= g.get_nodes_by_attr("label",nodes)
        print(nodes,graphNode) 
        graphNode_ID=graphNode[0]['id']
        print("In IF statement",type(graphNode_ID)) 
        for concentration in normVals[0:][specNo]: #Which column we are using here deosnt really matter but skip the time
	    g.add_event(timeline_counter, "graph:set_node_attribute",
            {
                'id': graphNode_ID,
                'name': 'og:space:activity',
                'type': 'float',
                'value': str(concentration)
            } 
            )
            timeline_counter+=100
        timeline_counter=0
    else:  #This is a species node so change all edges coming out of it to a width proportional to concentration +
        graphNode= g.get_nodes_by_attr("label",nodes)
        print("In else statement",nodes,graphNode) 
        graphNode_ID=graphNode[0]['id']
        edgeList=g.get_edges_by_attr('src',graphNode_ID)
        print ("edge list is",edgeList) 
        for edge in edgeList:  #each node could have several edges, need to deal with src and dst nodes sperately because of semanticnet limitations
            print ("edge is ",edge) 
            edge_ID=edge['id']
            for concentration in normVals[0:][specNo]: #Which column we are using here deosnt really matter but skip the time
	        g.add_event(timeline_counter, "graph:set_link_attribute",
                {
                'id': edge_ID,
                'name': 'og:space:width',
                'type': 'float',
                'value': str(concentration)
                } 
                ) 
                timeline_counter+=100 
            timeline_counter=0  
        g.cache_edges_by("dst")   
        edgeList=g.get_edges_by_attr('dst',graphNode_ID)
        print ("edge list2 is",edgeList,"graphNode_ID:",graphNode_ID) 
        for edge in edgeList:  #each node could have several edges, need to deal with src and dst nodes sperately because of semanticnet limitations 
            edge_ID=edge['id']    
            for concentration in normVals[0:][specNo]: #Which column we are using here deosnt really matter but skip the time
	        g.add_event(timeline_counter, "graph:set_link_attribute",
                {
                'id': edge_ID,
                'name': 'og:space:width',
                'type': 'float',
                'value': str(concentration)
                } 
                )  
                # print("In destination loop adding concentration")    
                timeline_counter+=100 
            timeline_counter=0     
    specNo+=1	
#Need to seperate reaction nodes and species nodes solutions.  Can do this as all reactions are labeled J_r# in solution files 
 

g.save_json(args.outfile)


