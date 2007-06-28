default(parisize, 100M)
default(realprecision, 100)

a = 501001011
b = 123012092

boven(x) =
{
return( max(abs(4*x*(x-a)*(x+b)), (x^2+a*b)^2) )
}

onder(x) =
{
return( max(1, abs(x)) )
}

f(x) = 
{
return( log(boven(x) / onder(x)^4) )
}

intnum(X = [-1],[1],f(X)/(1+X^2)) / Pi
