import random
################ Define strategy classes for each actor ############

# Insert Waffler Class definition here from lecture notes

class defect:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="defect"
    def response(self, other):
        return "Defect"
    def inform(self, other, other_response):
        return

class cooperate:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="cooperate"
    def response(self, other):
        return "Cooperate"
    def inform(self, other, other_response):
        return

class tit_for_tat:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="tit_for_tat"
        self.responses={}
        for i in range(Nactors):
            self.responses[i]="Cooperate"
    def response(self, other):
        return self.responses[other]
    def inform(self, other, other_response):
        self.responses[other]=other_response
        return

class susp_tit_for_tat:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="susp_tit_for_tat"
        self.responses={}
        for i in range(Nactors):
            self.responses[i]="Defect"
    def response(self, other):
        return self.responses[other]
    def inform(self, other, other_response):
        self.responses[other]=other_response
        return

class grudger:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="grudger"
        self.responses={}
        for i in range(Nactors):
            self.responses[i]="Cooperate"
    def response(self, other):
        return self.responses[other]
    def inform(self, other, other_response):
        if other_response=="Defect":
            self.responses[other]=other_response
        return

class coin_flip:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="coin_flip"
        self.ran=ranclass()
    def response(self, other):
        if self.ran.ranfloat()<.5:
            return "Defect"
        else:
            return "Cooperate"
    def inform(self, other, other_response):
        return    

class biased_random:
    def __init__(self,Nactors,myid):
        self.Nactors=Nactors
        self.myid=myid
        self.name="biased_random"
        self.ran=ranclass()
        self.responses=[]
        for i in range(Nactors):
            self.responses.append((0,0))
    def response(self, other):
        ndefect,ncooperate=self.responses[other]
        if ndefect==0 and ncooperate==0:
            return "Cooperate"
        if ndefect==0:
            return "Cooperate"
        if ncooperate==0:
            return "Defect"
        ratio=float(ndefect)/float(ndefect+ncooperate)
        if self.ran.ranfloat()<ratio:
            return "Defect"
        else:
            return "Cooperate"
    def inform(self, other, other_response):
        ndefect,ncooperate=self.responses[other]
        if other_response=="Defect":
            ndefect+=1
        else:
            ncooperate+=1
        self.responses[other]=ndefect,ncooperate

################ Utility class for managing strategy functions ##############
class actor:
    def __init__(self,myid,name,strategy):
        self.strategy=strategy
        self.myid=myid
        self.name=name
        self.ntrans=0
        self.score=0
    def pay(self,income):
        self.score+=income
        self.ntrans+=1
    def return_score(self):
        return self.score
    def return_ntrans(self):
        return self.ntrans
    def return_name(self):
        return self.name
    def return_myid(self):
        return myid
    
############### Define Python classes for simulation ########################
# This is a utility class for making random numbers
class ranclass:
    def ranint(self, n):
        return int(n*random.random())
    def ranfloat(self):
        return random.random()

