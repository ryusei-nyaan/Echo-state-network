module ESN
using LinearAlgebra
using SparseArrays

function fit(x,n,τ,pwij,pwin,pwself,pwb,S)
    
    b = sprandn(n,pwb)
    Wij = Array(sprandn(n,n,pwij))
    wij = eigvals(Wij)
    
    while maximum(abs,wij)>=1
        Wij = Array(sprandn(n,n,pwij))
        wij = eigvals(Wij)
    end
    
    X = zeros(2*n)
    Win = sprandn(n,pwin)
    Wself = sprandn(n,pwself)
    Wo = zeros(2*n,1)
    
    function esn_calc(u,win,wij,xx,wself,y,bb)
        U = win*u + wij*xx + wself*y + bb
        return tanh.(U)
    end
    
    x1 = esn_calc(x[1],Win,Wij,X[1:n],Wself,0,b)
    X[1:n] = x1
    X[n+1:2*n] = map(i->i^2,x1)
    E = Matrix(I,2*n,2*n)
    N = length(x)
    P = E*(1/τ)
    
    for i = 1:S
        for j = 1:N-1
            r = 1+transpose(X)*P*X
            k = P*X/r
            e = x[j+1] - (transpose(X)*Wo)[1]
            Wo = Wo + k*e
            P = P - (P*X*transpose(X)*P)/r
            x1 = esn_calc(x[j+1],Win,Wij,X[1:n],Wself,x[j],b)
            X[1:n] = x1
            X[n+1:2*n] = map(k->k^2,x1)
        end
    end
    
    return Win,Wij,Wself,b,Wo
end

function predict(x,xx,Win,Wij,Wself,b,Wo)
    
    n = length(Win)
    X = zeros(n*2)
    U = Win*x + Wij*X[1:n] + Wself*xx + b
    x1 = tanh.(U)
    X[1:n] = x1
    X[n+1:2*n] = map(s->s^2,x1)
    y = (transpose(X)*Wo)[1]
    return y
end

end

