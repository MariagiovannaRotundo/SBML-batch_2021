import roadrunner
import libsbml



def withAssignmentRule(rr, model, specie, reaction, nTest, valueIncrease, assignmentRule):

    formula = assignmentRule.getFormula()
    #print(formula)

    elements = assignmentRule.getMath().getListOfNodes()

    allSpeciesInRule = [model.getSpecies(elements.get(i).getName()) for i in range(elements.getSize()) if elements.get(i).isName() and model.getSpecies(elements.get(i).getName()) is not None]

 
    #use a time which is not 0
    rr.model.setTime(0.1)

    originalQuantities = rr.model.getFloatingSpeciesAmounts()
    originalConcentrationsBoundary = rr.model.getBoundarySpeciesConcentrations()

    #vedo per quali specie nell'assignment rule non c'è un'assignment rule
    for oneSpecie in allSpeciesInRule:

        #non considerare se la specie appare anche nella kinetic low perchè influenza i risultati
        i=0
        while i<model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize() and \
            oneSpecie.getId()!=model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().get(i).getName():
            i=i+1
        #check if the species it's not in law
        if i == model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize():

            #se non c'è la regola
            if model.getAssignmentRule(oneSpecie.getId()) is None:

                #vedo se la specie è boundary
                if oneSpecie.isSetBoundaryCondition() and oneSpecie.getBoundaryCondition():
                    role = analysisBoundary(rr, model, reaction, nTest, valueIncrease, oneSpecie)
                else:
                    role = analysisFloating(rr, model, reaction, nTest, valueIncrease, oneSpecie)

                if role == 0 or role == 1:
                    return role

    
    #se arriva qui non ha capito il ruolo del modificatore: provo ad incrementare il valore
    #delle altre specie cercando  di non avere specie a 0 nella regola così che non 
    #mi annullino la specie che sto incrementando ( es. 0/specie): 
    #incremento tutte le specie di un po' 

    
  
    #to try the role changing the quantities of molecules
    #before there is a less change then a major change 
    for j in range(nTest):
        quantities=originalQuantities.copy()
        
        for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
            if model.getAssignmentRule(specieId) is None:
                
                quantities[i]=quantities[i]+(0.1+j)
                rr.model.setFloatingSpeciesAmounts([i],[quantities[i]])
                
                #quantities are not changed if the change produce a negative rate
                if rr.model.getReactionRates()[reaction]<0:
                    quantities[i]=quantities[i]-(0.1+j)
                    rr.model.setFloatingSpeciesAmounts([i],[quantities[i]])


        concentrations = originalConcentrationsBoundary.copy()

        for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
            if model.getAssignmentRule(specieId) is None:
                
                concentrations[i]=concentrations[i]+(0.1+j)
                rr.model.setBoundarySpeciesConcentrations([i],[concentrations[i]])
                
                #concentrations are not changed if the change produce a negative rate
                if rr.model.getReactionRates()[reaction]<0:
                    concentrations[i]=concentrations[i]-(0.1+j)
                    rr.model.setBoundarySpeciesConcentrations([i],[concentrations[i]])

                
        if rr.model.getReactionRates()[reaction]<0:
            for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setFloatingSpeciesAmounts([i], [originalQuantities[i]])
            for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setBoundarySpeciesConcentrations([i], [originalConcentrationsBoundary[i]])
            continue

        #analyze the role of modifier and reset the original quantities of species

        #vedo per quali specie nell'assignment rule non c'è un'assignment rule
        for oneSpecie in allSpeciesInRule:

            #non considerare se la specie appare anche nella kinetic low perchè influenza i risultati
            i=0
            while i<model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize() and \
                oneSpecie.getId()!=model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().get(i).getName():
                i=i+1
            #check if the species it's not in law
            if i == model.getReaction(reaction).getKineticLaw().getMath().getListOfNodes().getSize():

                #se non c'è la regola
                if model.getAssignmentRule(oneSpecie.getId()) is None:

                    #vedo se la specie è boundary
                    if oneSpecie.isSetBoundaryCondition() and oneSpecie.getBoundaryCondition():
                        role = analysisBoundary(rr, model, reaction, nTest, valueIncrease, oneSpecie)
                    else:
                        role = analysisFloating(rr, model, reaction, nTest, valueIncrease, oneSpecie)

                    if role == 0 or role == 1:
                        return role

    return -1




#this meThod is called if there are not Assignment rules in the whole model!!!
def withFloatingAndBoundary(rr, model, specie, reaction, nTest, valueIncrease):

    #print("\nFORMULA!!!")
    #print(model.getReaction(reaction).getKineticLaw().getFormula())
    #print("MODIFICATORE: ", end="")
    #print(model.getReaction(reaction).getModifier(specie).getSpecies())
    #print()

    #original quantities of species
    originalBoundaryQuantities = rr.model.getBoundarySpeciesConcentrations()
    originalFloatingQuantities = rr.model.getFloatingSpeciesAmounts()
    
    #use a time which is not 0
    rr.model.setTime(0.1)

    #remove parameters equal to 0
    parametersIds = rr.model.getGlobalParameterIds()
    parametersValue = rr.model.getGlobalParameterValues()
    for (i, value) in enumerate(parametersValue):
        if value == 0:
            if model.getAssignmentRule(parametersIds[i]) is None:  
                rr.model.setGlobalParameterValues([i], [0.1])
    
    
    #analyze the role of modifier
    modifier = model.getSpecies(specie)
    
    if modifier.isSetBoundaryCondition() and modifier.getBoundaryCondition():
        role = analysisBoundary(rr, model, reaction, nTest, valueIncrease, modifier)

        for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
            if model.getAssignmentRule(specieId) is None:   
                rr.model.setBoundarySpeciesConcentrations([i], [originalBoundaryQuantities[i]])
        
    else:
        role = analysisFloating(rr, model, reaction, nTest, valueIncrease, modifier)
        
        for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
            if model.getAssignmentRule(specieId) is None:   
                rr.model.setFloatingSpeciesAmounts([i], [originalFloatingQuantities[i]])


    # if the role of modifier is known, it returns the corresponding code
    if role==0 or role==1:

        for (i, value) in enumerate(parametersValue):
            if model.getAssignmentRule(parametersIds[i]) is None:  
                rr.model.setGlobalParameterValues([i], [parametersValue[i]])

        return role
        


    for j in range(nTest):
        quantities=originalFloatingQuantities.copy()
        
        for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
            if model.getAssignmentRule(specieId) is None:
                
                quantities[i]=quantities[i]+(0.1+j)
                rr.model.setFloatingSpeciesAmounts([i],[quantities[i]])

                #quantities are not changed if the change produce a negative rate
                if rr.model.getReactionRates()[reaction]<0:
                    quantities[i]=quantities[i]-(0.1+j)
                    rr.model.setFloatingSpeciesAmounts([i],[quantities[i]])

                    
                    
        #print("specie ") 
        #print(rr.model.getFloatingSpeciesIds())
        #print(rr.model.getFloatingSpeciesAmounts())
        

        concentrations = originalBoundaryQuantities.copy()

        for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
            if model.getAssignmentRule(specieId) is None:

                concentrations[i]=concentrations[i]+(0.1+j)
                rr.model.setBoundarySpeciesConcentrations([i],[concentrations[i]])
                
                #concentrations are not changed if the change produce a negative rate
                if rr.model.getReactionRates()[reaction]<0:
                    concentrations[i]=concentrations[i]-(0.1+j)
                    rr.model.setBoundarySpeciesConcentrations([i],[concentrations[i]])

        

        if rr.model.getReactionRates()[reaction]<0:
            for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setFloatingSpeciesAmounts([i], [originalFloatingQuantities[i]])
            for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setBoundarySpeciesConcentrations([i], [originalBoundaryQuantities[i]])
            continue



        if modifier.isSetBoundaryCondition() and modifier.getBoundaryCondition():
            role = analysisBoundary(rr, model, reaction, nTest, valueIncrease, modifier)
            for (i,specieId) in enumerate(rr.model.getBoundarySpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setBoundarySpeciesConcentrations([i], [originalBoundaryQuantities[i]])
        else:
            role = analysisFloating(rr, model, reaction, nTest, valueIncrease, modifier)
            for (i,specieId) in enumerate(rr.model.getFloatingSpeciesIds()):
                if model.getAssignmentRule(specieId) is None:   
                    rr.model.setFloatingSpeciesAmounts([i], [originalFloatingQuantities[i]])


        # if the role of modifier is known, it returns the corresponding code
        if role==0 or role==1:

            for (i, value) in enumerate(parametersValue):
                if model.getAssignmentRule(parametersIds[i]) is None:  
                    rr.model.setGlobalParameterValues([i], [parametersValue[i]])

            return role

    for (i, value) in enumerate(parametersValue):
        if model.getAssignmentRule(parametersIds[i]) is None:  
            rr.model.setGlobalParameterValues([i], [parametersValue[i]])
    #if the role of modifier is not known
    return -1







"""
*********************************

METODI PER L'ANALISI DEL RUOLO!!

*********************************
"""



#to establish the role of a modifier having a fixed configuration 
#in input it takes reference to roadrunner, reference to model ceate using libsbml, number of specie modifiers,
#number of reaction, number of test to do, vector with quantities of species, value of increase of quantities.
#it returns 1 if it's a promoter, 0 if it's an inhibitor, -1 if don't known
def analysisFloating(rr, model, reaction, nTest, valueIncrease, specie):

    

    #index of specie to examine
    index=rr.model.getFloatingSpeciesIds().index(specie.getId())
    
    quantities = rr.model.getFloatingSpeciesAmounts([index])[0]

    if quantities == 0:
        rr.model.setFloatingSpeciesAmounts([index],[0.1])

    #print("specie ") 
    #print(rr.model.getFloatingSpeciesIds())
    #print(rr.model.getFloatingSpeciesAmounts())


    #original rates of sistem
    originalRate=rr.model.getReactionRates()[reaction]

    countProm=0
    countInhib=0

    #execute the specified number of tests 
    for i in range(1, nTest+1):
        #increase the quantity of specie to examine

        rr.model.setFloatingSpeciesAmounts([index],[quantities+quantities*(i*valueIncrease)/100])
        

        #if rates increase
        if rr.model.getReactionRates()[reaction]>originalRate:
            countProm=countProm+1
        #if rates decrease
        elif rr.model.getReactionRates()[reaction]<originalRate:
            countInhib=countInhib+1

        #print("RATE "+str(rr.model.getReactionRates()[reaction]))

        originalRate=rr.model.getReactionRates()[reaction]

    
    #print("Inibitore : "+ str(countInhib==nTest))
    #print("Promotore : "+ str(countProm==nTest))


    rr.model.setFloatingSpeciesAmounts([index],[quantities])

    #if rates decrease in every test it's a inhibitor
    if countInhib==nTest:
        return 0
    #if rates increase in every test it's a promoter
    elif countProm==nTest:
        return 1


    #using a major increase of quantity
    if quantities == 0:
        rr.model.setFloatingSpeciesAmounts([index],[5])
    
    #print("specie ") 
    #print(rr.model.getFloatingSpeciesIds())
    #print(rr.model.getFloatingSpeciesAmounts())


    originalRate=rr.model.getReactionRates()[reaction]
    countProm=0
    countInhib=0

    for i in range(1, nTest+1):
        rr.model.setFloatingSpeciesAmounts([index],[quantities*i*valueIncrease])
        
        if rr.model.getReactionRates()[reaction]>originalRate:
            countProm=countProm+1
        elif rr.model.getReactionRates()[reaction]<originalRate:
            countInhib=countInhib+1

        #print("RATE "+str(rr.model.getReactionRates()[reaction]))

        originalRate=rr.model.getReactionRates()[reaction]

    
    #print("Inibitore : "+ str(countInhib==nTest))
    #print("Promotore : "+ str(countProm==nTest))

    rr.model.setFloatingSpeciesAmounts([index],[quantities])


    if countInhib==nTest:
        return 0
    elif countProm==nTest:
        return 1
    
    #if the role is not known return -1
    return -1



def analysisBoundary(rr, model, reaction, nTest, valueIncrease, specie):

    #index of specie to examine
    index=rr.model.getBoundarySpeciesIds().index(specie.getId())
    
    quantities = rr.model.getBoundarySpeciesConcentrations([index])[0]

    if quantities == 0:
        rr.model.setBoundarySpeciesConcentrations([index],[0.1])

    #original rates of sistem
    originalRate=rr.model.getReactionRates()[reaction]

    countProm=0
    countInhib=0

    #execute the specified number of tests 
    for i in range(1, nTest+1):
        #increase the quantity of specie to examine

        rr.model.setBoundarySpeciesConcentrations([index],[quantities+quantities*(i*valueIncrease)/100])
        

        #if rates increase
        if rr.model.getReactionRates()[reaction]>originalRate:
            countProm=countProm+1
        #if rates decrease
        elif rr.model.getReactionRates()[reaction]<originalRate:
            countInhib=countInhib+1

        originalRate=rr.model.getReactionRates()[reaction]

    #print("RATE "+str(rr.model.getReactionRates()[reaction]))
    #print("Inibitore : "+ str(countInhib==nTest))
    #print("Promotore : "+ str(countProm==nTest))


    rr.model.setBoundarySpeciesConcentrations([index],[quantities])

    #if rates decrease in every test it's a inhibitor
    if countInhib==nTest:
        return 0
    #if rates increase in every test it's a promoter
    elif countProm==nTest:
        return 1


    #using a major increase of quantity
    if quantities == 0:
        rr.model.setBoundarySpeciesConcentrations([index],[2])
    
    originalRate=rr.model.getReactionRates()[reaction]
    countProm=0
    countInhib=0

    for i in range(1, nTest+1):
        rr.model.setBoundarySpeciesConcentrations([index],[quantities*i*valueIncrease])
        
        if rr.model.getReactionRates()[reaction]>originalRate:
            countProm=countProm+1
        elif rr.model.getReactionRates()[reaction]<originalRate:
            countInhib=countInhib+1

        originalRate=rr.model.getReactionRates()[reaction]

    #print("RATE "+str(rr.model.getReactionRates()[reaction]))
    #print("Inibitore : "+ str(countInhib==nTest))
    #print("Promotore : "+ str(countProm==nTest))

    rr.model.setBoundarySpeciesConcentrations([index],[quantities])


    if countInhib==nTest:
        return 0
    elif countProm==nTest:
        return 1
    
    #if the role is not known return -1
    return -1