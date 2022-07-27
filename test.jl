using LinearAlgebra
using Plots

function f(x, L)
    if x < 0
        return 0.0
    elseif x > L
        return 1.0
    else
        return exp(-(L/x)^2*exp(-L/(L-x)))
    end
end

N = 200
dx = 1/N
xx = 0:dx:1
L = 50*dx

plot(xx,f.(xx,L))