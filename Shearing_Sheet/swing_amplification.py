from ShearingSheet import * 



def amplification_different_k():
	lambdas = np.linspace(0.1, 5, 40)
	amp = []

	for l in lambdas:
		print(l)
		sheet = ShearingSheet(1/l, 1.2)
		amp.append(sheet.delta_amplification()) 

	plt.plot(lambdas, amp)
	plt.yscale('log', base = 10)
	plt.show()



#amplification_different_k()
