#!/usr/bin/python3
# Prisoner's dilemma simulator program
# By Prof. M. Colvin, UC Merced
# This software is in the Public Domain for all uses

# Bring in useful modules

import random
import math
import copy
from ipd_class import *

### Create Roulette selection function ###
# Adapted from stackoverflow.com/questions/298301/roulette-wheel-selection-algorithm
def roulette_select(fitnesses, num):
    total_fitness = float(sum(fitnesses))
    rel_fitness = [f/total_fitness for f in fitnesses]
    # Generate probability intervals for each individual
    probs = [sum(rel_fitness[:i+1]) for i in range(len(rel_fitness))]
    # Draw new population
    new_population = [0]*len(fitnesses)
    # Calculate even dist
    evendist = [ (f+1)/num for f in range(num) ]
    #print probs
    #print evendist
    
    epsilon=0.000001
    for n in range(num):
        #r = random.random()
        r=evendist[n]
        for i in range(len(fitnesses)):
            if r <= probs[i]+epsilon:
                new_population[i]+=1
                break
    return new_population

### Set up payoff matrix ###
payoff={("Defect","Defect"):(1,1),
        ("Cooperate","Defect"):(0,3),
        ("Defect","Cooperate"):(3,0),
        ("Cooperate","Cooperate"):(2,2)}

### Open file for population output ###
outfile=open("ipd_pop.txt","w")

### Simulation Settings ###
ntransactions=100000          # Number of two-actor transactions
ngenerations=40               # Number of generations

### Strategies available ###
strats=[defect,cooperate,tit_for_tat,coin_flip,biased_random,susp_tit_for_tat,grudger]

### Set list for the number of each strategy ###
Nactor_list=[5, 15, 20, 10, 10, 10, 10]
Nactors=sum(Nactor_list)

### Print these to output file
for i in strats:
    memo="%s "%i.__name__
    outfile.write(memo)
outfile.write("\n")

# Create random number class
ran=ranclass()

for gen in range(ngenerations):
    if Nactors<=1:
        print("At least two players must be defined")
        exit()
    
    actors=[]
    iactor=0

    # Construct list of lists to store interaction matrix
    transactions=[]
    teamname={}
    for i in range(len(Nactor_list)):
        transactions.append([])
        teamname[strats[i].__name__]=i
        for j in range(len(Nactor_list)):
            transactions[i].append(0)

    # Construct a set of actor objects with different strategy objects
    for i in range(len(Nactor_list)):
        for j in range(Nactor_list[i]):
            actors.append(actor(iactor,strats[i].__name__,strats[i](Nactors,iactor)))
            iactor+=1

    #for act in actors:
    #    print act.return_name()
    
    ### Run simulation ###
    for i in range(ntransactions):
        #Pick two actors at random
        pick2=pick1=ran.ranint(Nactors)
        while (pick2==pick1):
            pick2=ran.ranint(Nactors)

        #Update transaction matrix
        transactions[teamname[actors[pick1].return_name()]]\
            [teamname[actors[pick2].return_name()]]+=1

        #Calculate responses
        response1=actors[pick1].strategy.response(pick2)
        response2=actors[pick2].strategy.response(pick1)

        #Give payoff
        payoff1,payoff2=payoff[(response1,response2)]
        actors[pick1].pay(payoff1)
        actors[pick2].pay(payoff2)

        #Inform actors of the other actors' responses
        actors[pick1].strategy.inform(pick2,response2)
        actors[pick2].strategy.inform(pick1,response1)

    #print "Number of transactions",ntransactions

    #Create a set of dictionaries to hold team scores and team sizes        
    team_scores_totals={}
    team_scores={}
    team_size={}
    team_ntrans={}
    for i in range(Nactors):
        if actors[i].return_name() not in team_size:
            team_scores_totals[actors[i].return_name()]=actors[i].return_score()/actors[i].return_ntrans()
            team_size[actors[i].return_name()]=1
            team_ntrans[actors[i].return_name()]=actors[i].return_ntrans()
        else:
            team_scores_totals[actors[i].return_name()]+=actors[i].return_score()/actors[i].return_ntrans()
            team_size[actors[i].return_name()]+=1
            team_ntrans[actors[i].return_name()]+=actors[i].return_ntrans()

    #print "\nPlayers:"
    #for i in team_size:
    #    print "%19s: team size: %4d team transactions: %8d"%(i,team_size[i],team_ntrans[i])

    for i in team_scores_totals:
        #team_scores[i]=team_scores_totals[i]/team_size[i]
        team_scores[i]=team_scores_totals[i]

    scores=[]
    for i in range(len(strats)):
        if Nactor_list[i] != 0:
            scores.append(team_scores[strats[i].__name__])
        else:
            scores.append(0.0)

    #Update populations
    new_team_sizes=roulette_select(scores,Nactors)
    #print new_team_sizes
    Nactor_list=copy.copy(new_team_sizes)
    
    total_scores = float(sum(scores))
    rel_scores = [f/total_scores for f in scores]
    # Generate probability intervals for each individual
    probs = [sum(rel_scores[:i+1]) for i in range(len(rel_scores))]

    print("\n  Generation=%3d       Pop      Score     New Pop"%(gen))
    icount=0
    for i in range(len(strats)):
        if Nactor_list[i] != 0:
            #print "%19s: %3d %8.4f %3d %5.3f %5.3f"%(strats[i].__name__,team_size[strats[i].__name__],scores[i],new_team_sizes[i],rel_scores[i],probs[i])
            print("%19s: %3d      %8.4f     %3d"%(strats[i].__name__,team_size[strats[i].__name__],scores[i],new_team_sizes[i]))
            memo="%3d  "%(team_size[strats[i].__name__])
            outfile.write(memo)
        else:   
            #print "%19s: %3d %8.4f %3d %5.3f %5.3f"%(strats[i].__name__,0,0.0,0,0.0,0.0)
            print("%19s: %3d      %8.4f     %3d"%(strats[i].__name__,0,0.0,0))
            memo="%3d  "%(0)
            outfile.write(memo)
    outfile.write("\n")
    outfile.flush()
outfile.close()

