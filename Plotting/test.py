from generalPlottingFunctions import *




#new = readingInComplexCSV("kalnajsTest.csv")


n = 0
old = readingInComplexCSV("evolution200.csv")
plt.plot(np.linspace(0,200, 200), np.absolute(old[:,n]), label = 'Linear calculation')
new = readingInComplexCSV("other.csv")
plt.plot(np.absolute(new[:,n]), label = 'New n body')

new = readingInComplexCSV("kalnajsTest.csv")
#plt.plot(np.absolute(9.7*new[:,n]), label = 'Old n body')
plt.legend()
plt.show()
'''
n = 0
new = readingInRealCSV("particlePositionFront.csv")
plt.plot(new[:,n+0], new[:,n+1])



new = readingInRealCSV("particlePositionBack.csv")
plt.plot(new[:,n+0], new[:,n+1])


new = readingInRealCSV("particlePositionOld.csv")
plt.plot(new[:,n+4], new[:,n+5])
plt.plot(new[:,n+0], new[:,n+1],'--')

plt.show()

n = 4
data = readingInRealOUT("../nBody/particleSamples.csv")
radius = np.power(np.power(data[:,1],2)+np.power(data[:,2],2), 0.5)
angMom = (np.multiply(data[:,1], data[:,4]) + np.multiply(data[:,2], data[:,3]))/radius
vR = (radius*(data[:,3]+data[:,4]))/(data[:,1]+data[:,2])

plt.hist(data[:,4], bins =1000)
print(np.shape(data))
data = readingInRealOUT("../nBody/Bodies/particleSamples.out")
radius = np.power(np.power(data[:,1],2)+np.power(data[:,2],2), 0.5)
angMom = (np.multiply(data[:,1], data[:,4]) + np.multiply(data[:,2], data[:,3]))/radius
vR = (radius*(data[:,3]+data[:,4]))/(data[:,1]+data[:,2])


plt.hist(data[:,4], bins =1000, alpha = .5)
print(np.shape(data))
plt.show()'''