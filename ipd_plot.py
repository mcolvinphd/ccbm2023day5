import matplotlib.pyplot as plt
file=open("ipd_pop.txt","r")
line=file.readline()
names=line.split()
dict={}
for i in names:
    dict[i]=[]
maxpop=0    
npops=0
while 1:
    line=file.readline()
    if line=="":
        break
    values_str=line.split()
    values=[eval(i) for i in values_str]
    maxpop=max(maxpop,max(values))
    for (key,value) in zip(names,values):
        dict[key].append(value)
for i in names:
    plt.plot(dict[i],label=i)
plt.legend()
plt.show()

# name.list <- names(pops)
# pops[1,]
# ymax=max(pops)
# ymin=min(pops)
# plot(pops[,1],type="b",col=1,ylim=c(ymin,ymax),pch=1,lwd=2,ylab="Population",xlab="Generation")
# for (i in seq(2,length(name.list))) {
#     lines(pops[,i],col=i,lwd=2)
#     points(pops[,i],col=i,pch=i)
# }
# legend("topright",name.list,pch=c(1:length(name.list)),col=c(1:length(name.list)))
