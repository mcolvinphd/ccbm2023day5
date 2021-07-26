from sudoku_class import *
a=sudoku()
a.generate()
ntrials=100
for n in range(25,71,3):
    complete=0
    for i in range(ntrials):
        a.makepuzzle(n)
        a.solve()
        if a.solved():
            complete+=1
    print("For %d initial clues, completion percentage=%5.3f"%(n,complete/ntrials))
