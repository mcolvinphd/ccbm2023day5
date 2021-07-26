# Python classes for simulating card games
# By Prof. M. Colvin, UC Merced
# This software is in the Public Domain for all uses

# Bring in useful modules
import random
import math

################ Define Python classes for simulation ########################
# The class cell holds info about a particular card
class card:
    def __init__(self, suit, kind, value):
        self.cardsuit=suit
        self.cardtype=kind
        self.cardvalue=value
    def printcard(self): print("%2s %s %2d"%(self.cardtype, self.cardsuit, self.cardvalue))
    def __repr__(self): return "%2s %s %2d"%(self.cardtype, self.cardsuit, self.cardvalue)
    def __str__(self): return "%2s %s %2d"%(self.cardtype, self.cardsuit, self.cardvalue)
    def __eq__(self, acard):
        return (self.cardvalue==acard.cardvalue and self.cardsuit==acard.cardsuit)
    def type(self):
        return self.cardtype
    def value(self):
        return self.cardvalue
    def suit(self):
        return self.cardsuit
            
# This is a utility class for making random numbers
class ran:
    def ranint(self, n):
        return int(n*random.random())
    def ranfloat(self):
        return random.random()
    
# The class stack holds a stack of cards
class hand:
    def __init__(self):
        self.stack=[]
        self.ncards=0
    def printhand(self):
        print("Cards:")
        for card in self.stack:
            card.printcard() 
    def __repr__(self): return "Hand with %d cards"%(self.ncards)
    def append(self, card): self.stack.append(card)
    def pop(self): self.stack.pop()
    def card(self, index): return self.stack[index]
    def topcard(self): return self.stack[len(self.stack)-1]
    def showing_card(self): return self.stack[0]
    def have_11_ace(self):
        #print "Looking for an Ace"
        for card in range(len(self.stack)):
            if self.stack[card].cardtype=="A" and self.stack[card].cardvalue==11:
                #print "Found Ace in search"
                self.stack[card].cardvalue=1
                return True
        return False
    def value(self):
        total=0
        for card in self.stack:
            total+=card.value()
        return total
    def matchcard(self):
        if (len(self.stack) > 2):
            return self.stack[len(self.stack)-3]
        else: return card('N','0')
    def size(self): return len(self.stack)

# The class deck holds and manipulates a complete deck of cards
class deck:
    def __init__(self):
        self.deck=[]
        self.suits=['S','C','H','D']
        self.values={'A':1,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'10':10,'J':10,'Q':10,'K':10}
        self.kinds=['A','2','3','4','5','6','7','8','9','10','J','Q','K']
        for suit in self.suits:
            for kind in self.kinds:
                self.deck.append(card(suit,kind,self.values[kind]))
    def __repr__(self): return "Deck with %d cards"%(len(self.deck))
    def printdeck(self):
        i=0
        for card in self.deck:
            print("%3d "%(i), end=' ')
            card.printcard()
            i+=1
    def shuffle(self):
        newdeck=[]
        ncards=len(self.deck)
        ranclass=ran()
        while (ncards>0):
            nextcard=ranclass.ranint(ncards)
            #print "appending:", self.deck[nextcard].printcard()
            newdeck.append(self.deck[nextcard])
            self.deck.remove(self.deck[nextcard])
            ncards-=1
        self.deck=newdeck
    def dealcard(self): return self.deck.pop()
    def cardsleft(self): return len(self.deck)

