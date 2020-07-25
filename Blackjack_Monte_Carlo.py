# Blackjack Monte Carlo simulator program
# By Prof. M. Colvin, UC Merced
# This software is in the Public Domain for all uses

# This program implements a simple Blackjack Card game
# inside a Monte Carlo framework that allows one to compute
# statistics on the probability of winning with different strategies

# Bring in useful modules

from Blackjack_classes import *

##############################################################################################
# Simulation Settings
ntrials=10000        # Number of Monte Carlo Trials
Print_Hands=False    # Flag setting to print each hand

#Initialize win/loss counters
ndealerwins=0
nplayerwins=0

for trial in range(ntrials):
    #Create and shuffle card deck
    adeck=deck()
    adeck.shuffle()

    #Create new dealer and player
    dealer=dealerclass()
    player=playerclass()    

    #Optional--modify player's hold limit for different dealer showing card
    newholdlimit={'A':18,'2':15,'3':15,'4':15,'5':14,'6':14,'7':17,\
                      '8':17,'9':18,'10':18,'J':17,'Q':17,'K':17}
    player.set_holdlimit(newholdlimit)

    #Create dealer hand and add the first two cards to the hand
    dealerhand=hand()
    dealerhand.append(adeck.dealcard())
    dealerhand.append(adeck.dealcard())
    
    #Create dealer hand and add the first two cards to the hand
    playerhand=hand()
    playerhand.append(adeck.dealcard())
    playerhand.append(adeck.dealcard())

    #Play out the player's hand
    if playerhand.value()<player.holdvalue(dealerhand.showing_card()):
        while(player.return_status()=="Hit"):
            playerhand.append(adeck.dealcard())
            player.play(playerhand, dealerhand.showing_card())
    
    #Play out the dealer's hand
    if dealerhand.value()<dealer.holdvalue() and player.return_status()!="Busted":
        while(dealer.return_status()=="Hit"):
            dealerhand.append(adeck.dealcard())
            dealer.play(dealerhand)

    #Figure out who wins
    playeroutcome="Lost"
    dealeroutcome="Lost"
    if player.return_status()=="Busted":
        ndealerwins+=1
        dealeroutcome="Won"
    elif dealer.return_status()=="Busted":
        nplayerwins+=1
        playeroutcome="Won"
    elif playerhand.value()>dealerhand.value():
        nplayerwins+=1
        playeroutcome="Won"
    else:
        ndealerwins+=1
        dealeroutcome="Won"

    if Print_Hands:
        print("=======================================================")
        print("Player-->",playeroutcome)
        playerhand.printhand()
        print("Dealer-->", dealeroutcome)
        dealerhand.printhand()

        
print("=======================================================")
print("Ntrials=", ntrials)
print("Player wins:",nplayerwins)
print("Dealer wins:",ndealerwins)
print('Player wins: %4.2lf percent' %(nplayerwins/(nplayerwins+ndealerwins)*100))
