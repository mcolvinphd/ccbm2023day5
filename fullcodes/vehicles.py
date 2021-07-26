class vehicle:
    def __init__(self,doors,wheels,name):
        self.doors=doors
        self.wheels=wheels
        self.name=name
    def __repr__(self):
        return "%s: %d doors and %d wheels"%(self.name,self.doors,self.wheels)
    def __eq__(self,aveh):
        return self.doors==aveh.doors and self.wheels==aveh.wheels
    def __iadd__(self,aveh):
        self.doors+=aveh.doors
        self.wheels+=aveh.wheels
        self.name+="+"
        self.name+=aveh.name
        return self
    
    
        
