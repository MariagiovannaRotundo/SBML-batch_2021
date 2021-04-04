import roadrunner
import libsbml
import os

import ModifiersRole
from SBML_batch import BaseFunctions


#to create Petri nets of a set of models
#in input it takes reference to roadrunner, path of directory that contains the files, path
#of directory that will contain the new files, number of test to classify a modifier, value of increase of quantities
def createPetriNets(rr, path, pathResult, nTest, valueIncrease):
    files=os.listdir(path)
    for file in files:
        if file.endswith('.xml'):
            print(file)
            #to create a Petri net
            createOnePetriNet(rr, path, file, pathResult, nTest, valueIncrease)


#to create the Petri Net of model 
#in input it takes reference to roadrunner, path of directory that contains the file, file's name, path
#of directory that will contain the new file, number of test to classify a modifier, value of increase of quantities
#It returns 1 if an error occurs, 0 otherwise
def createOnePetriNet(rr, path, file, pathResult, nTest, valueIncrease):

    filefile = open("modifiers_discovered.txt", "a")
    filefilelist = []

    if BaseFunctions.createDirectory(pathResult)==1:
        return 1
    
    #it clears the currently loaded model and all associated memory
    rr.clearModel()

    path=path+"\\"

    try:
        #to load the file with roadrunner
        rr.load(path+file)
    except:
        print("An error occurs during file load")
        return 1
    
    #to get the model using libSBML
    model = libsbml.SBMLReader().readSBMLFromString(rr.getSBML()).getModel()

    nameFile=file[0:(len(file)-4)]
    try:
        f=open(pathResult+"\\"+nameFile+".gv", "w")
    except IOError:
        print("An error occurs during opening file")
        return 1

    #string that will contains petri net of model
    graph='digraph "'+str(file)+'"{\n'

    #to create a node for each specie
    for i in range(model.getNumSpecies()):
        graph=graph+str(model.getSpecies(i).getId())+" [shape=circle];\n"
    
    #to examine one reaction at a time
    for i in range(model.getNumReactions()):

        #to create a box (node) for each reaction
        graph=graph+str(model.getListOfReactions().get(i).getId())+" [shape=box];\n"

        #to get reaction's reactants
        for reactant in model.getReaction(i).getListOfReactants():
            graph=graph+str(reactant.getSpecies())+" -> "+str(model.getReaction(i).getId())+" [arrowhead=vee, label="+str(reactant.getStoichiometry())+"];\n"
        
        #to get reaction's products
        for product in model.getListOfReactions().get(i).getListOfProducts():
            graph=graph+str(model.getReaction(i).getId())+" -> "+str(product.getSpecies())+" [arrowhead=vee, label="+str(product.getStoichiometry())+"];\n"

        #lista di tutti i modificatori
        listModifiers = [modifier.getSpecies() for modifier in model.getReaction(i).getListOfModifiers()]

        #to get reaction's modifiers and establish if they are promoter or inhibitor
        for j in range(model.getReaction(i).getNumModifiers()):
            
            #se il modificatore è anche un reagente non lo considero
            if model.getReaction(i).getModifier(j).getSpecies() in [reactant.getSpecies() for reactant in model.getReaction(i).getListOfReactants()]:
                continue

            #se il modificatore compare più volte lo considero solo l'ultima volta
            if model.getReaction(i).getModifier(j).getSpecies() in listModifiers[j+1:]:
                continue
            
            response = findRole(rr, model, model.getReaction(i).getModifier(j).getSpecies(), i, nTest, valueIncrease)
            
            #if it's a promoter
            if response==1:
                graph=graph+str(model.getReaction(i).getModifier(j).getSpecies())+" -> "+str(model.getReaction(i).getId())+" [arrowhead=dot];\n"
            #if it's a inhibitor
            elif response==0:
                graph=graph+str(model.getReaction(i).getModifier(j).getSpecies())+" -> "+str(model.getReaction(i).getId())+" [arrowhead=tee];\n"    
            


        listReactants = [reactant.getSpecies() for reactant in model.getReaction(i).getListOfReactants()]
        listProducts = [product.getSpecies() for product in model.getReaction(i).getListOfProducts()]
        
        alreadySeen = []

        #vedo se ci sono specie nè reagenti nè modificatori nella legge cinetica
        for k in range(model.getReaction(i).getKineticLaw().getMath().getListOfNodes().getSize()):
            #se il nodo è una specie
            nameOfNode = model.getReaction(i).getKineticLaw().getMath().getListOfNodes().get(k).getName()
            if nameOfNode in [specie.getId() for specie in model.getListOfSpecies()]:
                #se non è nè nella lista dei reagenti nè dei modificatore trattalo come un modificatore
                if nameOfNode in alreadySeen or nameOfNode in listReactants or  nameOfNode in listModifiers or \
                    (nameOfNode in listProducts and model.getReaction(i).getReversible()):
                    continue

                alreadySeen.append(nameOfNode)
                
                response = findRole(rr, model, nameOfNode, i, nTest, valueIncrease)
                #if it's a promoter
                if response==1:
                    graph=graph+str(model.getReaction(i).getKineticLaw().getMath().getListOfNodes().get(k).getName())+" -> "+str(model.getReaction(i).getId())+" [arrowhead=dot];\n"
                #if it's a inhibitor
                elif response==0:
                    graph=graph+str(model.getReaction(i).getKineticLaw().getMath().getListOfNodes().get(k).getName())+" -> "+str(model.getReaction(i).getId())+" [arrowhead=tee];\n"    
                if file not in filefilelist:
                    filefilelist.append(file)
            
        #if there is the inverse reaction
        if model.getReaction(i).getReversible():
            #to create a box (node) for the inverse reaction
            graph=graph+str(model.getReaction(i).getId())+"_Reverse [shape=box];\n"

            #to get inverse reaction's reactants
            for reactant in model.getListOfReactions().get(i).getListOfReactants():
                graph=graph+str(model.getListOfReactions().get(i).getId())+"_Reverse -> "+str(reactant.getSpecies())+" [arrowhead=vee, label="+str(reactant.getStoichiometry())+"];\n"

            #to get inverse reaction's products
            for product in model.getListOfReactions().get(i).getListOfProducts():
                graph=graph+str(product.getSpecies())+" -> "+str(model.getListOfReactions().get(i).getId())+"_Reverse [arrowhead=vee, label="+str(product.getStoichiometry())+"];\n"

    graph=graph+"}"
    
    try:
        f.write(graph)
    finally:
        try:
            f.close()
        except IOError:
            return 1


    for elem in filefilelist:
        filefile.write(elem+"\n")

    filefile.close()

    return 0





"""
in input it takes reference to roadrunner, reference to model ceate using libsbml, number of specie modifiers,
number of reaction, number of test to do, value of increase of quantities.
output: 
if it's an inhibitor it returns 0,
if promoter 1,
if it fails to establish the role -1,
if the species it's not in law -2
"""
def findRole(rr, model, specieId, reaction, nTest, valueIncrease):

    
    i=0
    while i<model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize() and \
        specieId!=model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().get(i).getName():
        i=i+1

    #check if the species it's not in law
    if i==model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize():
        return -2




    #controllo se c'è un assignment rule per la specie modificatore
    assignmentRule = model.getAssignmentRule(specieId)
    if assignmentRule is not None:
        return ModifiersRole.withAssignmentRule(rr, model, specieId, reaction, nTest, valueIncrease, assignmentRule)

    return ModifiersRole.withFloatingAndBoundary(rr, model, specieId, reaction, nTest, valueIncrease)










