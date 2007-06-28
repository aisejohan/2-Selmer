f = x*(x-a)*(x+b)

phi = ( (x^2+a*b)^2 ) / (4*f) 

ee = ellinit([0,-a+b,0,-b*a,0],1)

twoxy = elladd(ee,[x,y],[x,y])

twox = twoxy[1]

print("The following two values should be zero:")

( 1/y^2 ) * ( \
   polcoeff(twox * y^2, 2, y)*y^2 + \
   polcoeff(twox * y^2, 0, y) \
            ) - twox


( 1/f ) * ( \
   polcoeff(twox * y^2, 2, y)*f + \
   polcoeff(twox * y^2, 0, y) \
          ) - phi

print("This proves that Lucien is right.")
