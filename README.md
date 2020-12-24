# Echo-state-network
This is one of reservoir computing methods.

###### Julia version = 1.5.1 


### fit(x,n,τ,pwij,pwin,pwself,pwb,S)
x is time-series. n is the count of nodes. τ is the regularization parameter. pwij, pwin, pwself, pwb is the ratio of non-zero in each sparce matrix which is "reservoir layers" and "bias terms". S is the count of optimization, in other words, it is "epoch".
What we should pay attention is that the ratio of non-zero is not accurate, so the real ratio might be less or more than the ratio we decide. I plan to revise and update this someday.
