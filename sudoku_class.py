### Class library for Soduko program        ###
### Solution algorithm inspired by the      ###
### book Programming Sudoku by Wei-Meng Lee ###
# This material is released under a Creative Commons License
# "Attribution-Noncommercial-Share Alike License"
# http://creativecommons.org/licenses/by-nc-sa/3.0/
# with the following attribution:
# "Original program written by Prof. Mike Colvin, UC Merced"
import sys
import random
class cell:
    def __init__(self):
        self.list=[1,2,3,4,5,6,7,8,9]
    def n(self):
        return len(self.list)
    def set(self,val):
        self.list=[val]
    def reset(self):
        self.list=[1,2,3,4,5,6,7,8,9]
    def setlist(self,vallist):
        self.list=[]
        self.list.extend(vallist)
    def val(self):
        if len(self.list)==1:
            return self.list[0]
        else:
            print("Trouble: Trying to return ambiguous value")
            sys.exit("Trouble: Trying to return ambiguous value")
    def vals(self):
        vals=[]
        vals.extend(self.list)
        return vals
    def pairs(self):
        if len(self.list)==2:
            return (self.list[0],self.list[1])
        else:
            return 0
    def triplets(self):
        if len(self.list)==3:
            return (self.list[0],self.list[1],self.list[2])
        else:
            return 0
    def remove(self, val):
        #print self.list.count(val)
        if self.list.count(val)==1:
            if len(self.list)==1:
                return
            else:
                self.list.remove(val)
            return
        else:
            return

    def add(self,val):
        if val < 1 or val > self.size:
            print("Trouble: Trying to add value outside range")
            sys.exit("Trouble: Trying to add value outside range")
        if self.list.count(val)==0:
            self.list.append(val)
            self.list.sort()
            return
        else:
            print("Trouble: Trying to add existing value")
            sys.exit("Trouble: Trying to add existing value")
    def contain(self, val):
        if self.list.count(val)==1:
            return True
        else:
            return False
    def dispchar(self):
        if (len(self.list)==1):
            return str(self.list[0])
        if (len(self.list)==0):
            return 'x'
        return '-'       
    def dispval(self):
        dispval=0
        for i in self.list:
            dispval*=10
            dispval+=i
        return dispval
            
class block:
    def __init__(self):
        self.cells=[]
    def add(self, cell):
        self.cells.append(cell)
    def find_singles(self):
        singles=[]
        for i in self.cells:
            if i.n()==1:
                singles.append(i.val())
        return singles
    def lone_val(self):
        vals=[]
        for i in self.cells:
            vals.extend(i.vals())
        #print vals
        for i in range(1,10):
            if vals.count(i)==1:
                #print "eliminating",i
                for j in self.cells:
                    if j.contain(i):
                        #print "set lv"
                        j.set(i)
                #self.display()
    def singles(self):
        singles=[]
        for i in self.cells:
            if i.n()==1:
                singles.append(i.val())
        for j in singles:
            for i in self.cells:
                i.remove(j)
    def pairs(self):
        pairs=[]
        icount=0
        for i in self.cells:
            if i.pairs()!=0:
                #print i.pairs()
                pairs.append((icount,i.pairs()))
            icount+=1
        #print "New list"
        if icount<2: return
        pairsort=sorted(pairs,key=lambda pair:pair[1])
        #print "Pair sort"
        #print pairsort
        while len(pairsort)>1:
            if pairsort[0][1]==pairsort[1][1]:
                self.rm_pairs((pairsort[0][0],pairsort[1][0]),pairsort[0][1])
                pairsort.remove(pairsort[0])
                pairsort.remove(pairsort[0])
            else:
                pairsort.remove(pairsort[0])
        return
    def triplets(self):
        triplets=[]
        icount=0
        for i in self.cells:
            if i.triplets()!=0:
                #print i.pairs()
                triplets.append((icount,i.triplets()))
            icount+=1
        #print "New list"
        if icount<3: return
        tripletsort=sorted(triplets,key=lambda pair:pair[1])
        #print "Pair sort"
        #print pairsort
        while len(tripletsort)>2:
            if tripletsort[0][1]==tripletsort[1][1]:
                if tripletsort[0][1]==tripletsort[2][1]:
                    self.rm_pairs((tripletsort[0][0],tripletsort[1][0],tripletsort[2][0]),tripletsort[0][1])
                    tripletsort.remove(tripletsort[0])
                    tripletsort.remove(tripletsort[0])
                    tripletsort.remove(tripletsort[0])
                else:
                    tripletsort.remove(tripletsort[0])
                    tripletsort.remove(tripletsort[0])
            else:
                    tripletsort.remove(tripletsort[0])
        return
    def rm_pairs(self,skip,vals):
       for i in range(len(self.cells)):
           if skip.count(i)==0:
               for j in vals:
                   self.cells[i].remove(j)
    def display(self):
        for i in self.cells:
            print("%1s "%(i.dispchar()), end=' ')
        print("")
    def testcomplete(self):
        values=[]
        for i in self.cells:
            if i.n()!=1:
                return False
            values.append(i.val())
        if len(set(values))==len(self.cells):
            return True
        return False
    def checkblock(self):
        vals=[]
        for i in self.cells:
            if i.n()==1:
                vals.append(i.val())
        for i in range(1,10):
            if vals.count(i)>1:
                print("Error count greater than 1 in a block")
                self.display()
    def digitcount(self):
        count=0
        for i in self.cells:
            count+=i.n()
        return count
    
class sudoku:
    def __init__(self):
        self.size=9
        self.cells=[]
        self.rows=[]
        self.cols=[]
        self.grids=[]
        for i in range(self.size):
            self.cells.append([])
            self.rows.append(block())
            self.cols.append(block())
            self.grids.append(block())
            for j in range(self.size):
                self.cells[i].append(cell())
        for i in range(self.size):
            for j in range(self.size):
                self.rows[i].add(self.cells[i][j])
                self.cols[i].add(self.cells[j][i])
                igrid=i//3
                jgrid=j//3
                ijblock=3*igrid+jgrid
                self.grids[ijblock].add(self.cells[i][j])
                #print "i=%d, j=%d, igrid=%d, jgrid=%d, ijblock=%d"%(i,j,igrid,jgrid,ijblock)
    def reset(self):
        self.cells=[]
        self.rows=[]
        self.cols=[]
        self.grids=[]
        for i in range(self.size):
            self.cells.append([])
            self.rows.append(block())
            self.cols.append(block())
            self.grids.append(block())
            for j in range(self.size):
                self.cells[i].append(cell())
        for i in range(self.size):
            for j in range(self.size):
                self.rows[i].add(self.cells[i][j])
                self.cols[i].add(self.cells[j][i])
                igrid=i//3
                jgrid=j//3
                ijblock=3*igrid+jgrid
                self.grids[ijblock].add(self.cells[i][j])     
    def initialize(self,initlist):
        self.reset()
        for i in initlist:
            self.cells[i[0]][i[1]].set(i[2])
    def display(self):
        print("\nSudoku cell values")
        for i in range(self.size):
            for j in range(self.size):
                print("%1s"%(self.cells[i][j].dispchar()), end=' ')
            print("")
    def dispcells(self):
        print("\nSudoku cell lists")
        for i in range(self.size):
            for j in range(self.size):
                print("%10d"%(self.cells[i][j].dispval()), end=' ')
            print("") 
    def rm_singles(self):
        for i in self.cols:
            i.singles()
        for i in self.rows:
            i.singles()
        for i in self.grids:
            i.singles()
    def rm_pairs(self):
        for i in self.cols:
            i.pairs()
        for i in self.rows:
            i.pairs()
        for i in self.grids:
            i.pairs()
    def rm_triplets(self):
        for i in self.cols:
            i.triplets()
        for i in self.rows:
            i.triplets()
        for i in self.grids:
            i.triplets()
    def rm_lone_vals(self):
        for i in self.cols:
            i.lone_val()
        #self.dispcells()
        for i in self.rows:
            i.lone_val()
        #self.dispcells()
        for i in self.grids:
            i.lone_val()
        #self.dispcells()
    def count_singles(self):
        singlecount=0
        for i in range(self.size):
            for j in range(self.size):
                if self.cells[i][j].n()==1:
                    singlecount+=1
        return singlecount
    def countdigits(self):
        count=0
        for i in self.cols:
            count+=i.digitcount()
        return count
    def lone_vals(self):
        olddigitcount=self.countdigits()
        while True:
            self.rm_lone_vals()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
    def singles(self):
        olddigitcount=self.countdigits()
        while True:
            self.rm_singles()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
    def pairs(self):
        olddigitcount=self.countdigits()
        while True:
            self.rm_pairs()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
        #if oldsinglecount-enteringsinglecount:
        #    print "pair single counts",enteringsinglecount,oldsinglecount
    def triplets(self):
        olddigitcount=self.countdigits()
        while True:
            self.rm_triplets()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
    def solve_2(self):
        olddigitcount=self.countdigits()
        while True:
            self.singles()
            self.lone_vals()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
    def solve_3(self):
        olddigitcount=self.countdigits()
        while True:
            self.singles()
            self.lone_vals()
            self.pairs()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()
    def solve_4(self):
        olddigitcount=self.countdigits()
        while True:
            self.singles()
            self.lone_vals()
            self.pairs()
            self.triplets()
            if self.countdigits()==olddigitcount:
                break
            olddigitcount=self.countdigits()  
    def solve(self):
        self.solve_4()
    def solved(self):
        return self.testcomplete()
    def testcomplete(self):
        for i in self.cols:
            if not i.testcomplete():
                #print "Col not complete",
                #i.display()
                return False
        for i in self.rows:
            if not i.testcomplete():
                #print "Row not complete",
                #i.display()
                return False
        for i in self.grids:
            if not i.testcomplete():
                #print "Grid not complete",
                #i.display()
                return False
        return True
    def try_generate(self):
        self.reset()
        for i in range(self.size):
            for j in range(self.size):
                singles=list(set(self.rows[i].find_singles()))
                singles=list(set(singles) | set(self.cols[j].find_singles()))
                igrid=i//3
                jgrid=j//3
                ijblock=3*igrid+jgrid
                singles=list(set(singles) | set(self.grids[ijblock].find_singles()))
                available=[1,2,3,4,5,6,7,8,9]
                for k in singles:
                    available.remove(k)
                #print available
                if len(available)==0:
                    #print "No choices available for cell[%d][%d]"%(i,j)
                    #self.dispcells()
                    #sys.exit("No options")
                    #10/0
                    if i>7:
                        print("i equals",i)
                        self.dispcells()
                    return False
                else:
                    self.cells[i][j].set(random.choice(available))
        return True
    def generate(self):
        while True:
            if self.try_generate():
                return
    def makepuzzle(self,n):
        self.generate()
        if not self.testcomplete():
            print("Didnt get good sudoku from generate()")
        #self.display()
        for k in range(self.size*self.size-n):
            while True:
                i=random.randint(0,8)
                j=random.randint(0,8)
                if self.cells[i][j].n()==1:
                    break
            self.cells[i][j].reset()
        #self.dispcells()
    def check(self):
        for i in self.cols:
            i.checkblock()
        for i in self.rows:
            i.checkblock()
        for i in self.grids:
            i.checkblock()
    def save(self,filename):
        f=file(filename,"w")
        for i in range(self.size):
            for j in range(self.size):
                string="%1s "%(self.cells[i][j].dispchar())
                f.write(string)
            f.write("\n")
        f.close()
    def read(self,filename):
        f=open(filename,"r")
        if not f.closed:
            for i in range(self.size):
                line=f.readline()
                vals=line.split()
                for j in range(self.size):
                    if vals[j]=="-":
                        self.cells[i][j].reset()
                    else:
                        self.cells[i][j].set(int(vals[j]))
    def savestate(self,filename):
        f=file(filename,"w")
        for i in range(self.size):
            for j in range(self.size):
                string="%1s "%(self.cells[i][j].dispchar())
                f.write(string)
            f.write("\n")
        f.close()    


    
    
