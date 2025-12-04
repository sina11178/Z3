We run L=9 simulations through here:

120 Cores --> 60 for (0.1, 0) and 60 for (0.1, pi/6)

Will run for: arange(0, 30, 0.5) = [0, .5, 1, 1.5, ... 29, 29.5]

Gamma values: [0, .2, .4, .6, .8, 1]

Every 10 cores is a new Gamma value: 
	b0 - b9 --> FOR gamma = 0
	b10 - b19 --> FOR gamma = 0.2
	b20 - b29 --> FOR gamma = 0.4
	
	AND SO ON

Each core will calculate 6 data points:
	b0 --> log(t) = 0, .5, 1, 1.5, 2, 2.5
	b1 --> log(t) = 3, 3.5, 4, 4.5, 5, 5.5

	AND SO ON 


