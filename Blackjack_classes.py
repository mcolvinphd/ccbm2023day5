# Blackjack Monte Carlo simulator program
# By Prof. M. Colvin, UC Merced
# This software is in the Public Domain for all uses

# This program implements a simple Blackjack Card game
# inside a Monte Carlo framework that allows one to compute
# statistics on the probability of winning with different strategies

# Bring in useful modules

import random
import math
import copy

################ Define Python classes for simulation ########################
class playerclass:
    def __init__(self):
        self.holdlimit={'A':18,'2':15,'3':15,'4':15,'5':14,'6':14,'7':17,\
                        '8':17,'9':18,'10':18,'J':17,'Q':17,'K':17}
        self.status="Hit"
    def return_status(self):
        return self.status
    def holdvalue(self,dealer_showing_card):
        return self.holdlimit[dealer_showing_card.type()]
    def set_holdlimit(self, holdlimit):
        self.holdlimit=copy.copy(holdlimit)
    def play(self, hand, dealer_showing_card):
        self.status="Hit"
        if hand.value()>21:
            if not hand.have_11_ace():
                self.status="Busted"
                return
            elif hand.value()>=self.holdlimit[dealer_showing_card.type()]:
                self.status="Hold"
        if hand.value()>=self.holdlimit[dealer_showing_card.type()]:
        #if hand.value()>=17:
            self.status="Hold"
        return
    
# The class dealer encompasses the dealer's strategy
class dealerclass:
    def __init__(self):
        self.holdlimit=17
        self.status="Hit"
    def return_status(self):
        return self.status
    def holdvalue(self):return self.holdlimit
    def play(self, hand):
        self.status="Hit"
        if hand.value()>21:
            if not hand.have_11_ace():
                self.status="Busted"
                return
            elif hand.value()>=self.holdlimit:
                self.status="Hold"
            #   print "I have an ace"
            #   hand.printhand()
        if hand.value()>=self.holdlimit:
            self.status="Hold"
        return

# The class cell holds info about a particular card
class card:
    def __init__(self, suit, type, value):
        self.cardsuit=suit
        self.cardtype=type
        self.cardvalue=value
    def printcard(self): print(self.cardtype, self.cardsuit, self.cardvalue)
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
        suits=['S','C','H','D']
        values={'A':11,'2':2,'3':3,'4':4,'5':5,'6':6,'7':7,'8':8,'9':9,'10':10,'J':10,'Q':10,'K':10}
        types=['A','2','3','4','5','6','7','8','9','10','J','Q','K']
        for suit in suits:
            for type in types:
                self.deck.append(card(suit,type,values[type]))
    def printdeck(self):
        i=0
        for card in self.deck:
            print(i, end=' ')
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

# Utility class for building histograms
class histogram:         
    def __init__(self, nbins, binstart, binend):
        self.maxbars=30     #Maximum number bars in printed histogram
        self.nbins=nbins
        self.binstart=binstart
        self.binend=binend
        binwidth=(float(binend)-float(binstart))/float(nbins)
        self.bins=[]
        self.sum_x=0.0
        self.sum_x2=0.0
        self.n=0
        self.bins.append(binstart)
        for i in range(nbins):
            self.bins.append((i+1)*binwidth)
        self.data=[0 for i in range(nbins)]
    def add(self,data):
        if (data < self.binstart or data > self.binend):
            print("data outside of histogram range")
        i=1
        while (data >= self.bins[i]): i+=1
        self.data[i-1]+=1
        self.sum_x+=data
        self.sum_x2+=data*data
        self.n+=1
    def print_data(self):
        for i in range(self.nbins):
            print('%2.2lf-%2.2lf: %-4d' % (self.bins[i],self.bins[i+1],
                                           self.data[i]))
    def print_hist(self):
        maxdata=0;
        mean=float(self.sum_x)/float(self.n)
        stdev=math.sqrt(float(self.n)*self.sum_x2-self.sum_x*self.sum_x)/float(self.n)
        for i in range(self.nbins):
           if self.data[i]>maxdata:
               maxdata=self.data[i]
        print('Mean: %4.2lf' % (mean))
        print('Standard deviation: %4.2lf' % (stdev))
        for i in range(self.nbins):
            bars=int(self.maxbars*float(self.data[i])/float(maxdata))
            print('%2.2lf-%2.2lf: %-4d %s' % (self.bins[i],self.bins[i+1],
                                           self.data[i],bars*"="))
  
