from cards import *
ntrials=10000
awins=0
for i in range(ntrials):
    adeck=deck()
    adeck.shuffle()
    ascore=0
    bscore=0
    while adeck.cardsleft()>2:
        acard1=adeck.dealcard()
        acard2=adeck.dealcard()
        bcard=adeck.dealcard()
        if acard1.value()>bcard.value() or acard2.value()>bcard.value():
            ascore+=1
        else:
            bscore+=1
    if ascore>bscore:
        awins+=1
print("Player A win percentage=",awins/ntrials)
