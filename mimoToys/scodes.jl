function scodes(M,Ctype)
#C = scodes(M,Ctype)
#  Return a constellation of M signals, where <Ctype> is 'PSK', 'QAM',
#  or 'PAM'.
#[C,d] = scodes(M,Ctype)
#  Also return the minimum distance between constellation points.
#[C,d,Cr] = scodes(M,Ctype)
#  Returns the underlying real PAM constellation for square (M=m^2)
#  QAM constellations

# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Fri 10 Apr 15, 10:58pm by peel
#

#if nargin==0
#  error("Bad arguments. Try again.")
#end
if Ctype=="PSK"
    C=exp(im*2*pi*[1:M]'/M);
elseif Ctype=="QAM"
    if abs(sqrt(M)-round(sqrt(M)))< 1e-3
        # Make QAM constellations if M is a square
        m = round(Int,sqrt(M))
        C = complex(zeros(m,m))
        Cr = complex([-(m-1):2:(m-1);])
        for ix = 1:m
            for jx = 1:m
                C[jx,ix] = Cr[ix] +im*Cr[jx];
            end
        end
    else
        error("Can only handle square QAM constellations.")
        #     [x y] = qaskenco(M);
        #     C = x+im*y;
    end
elseif Ctype=="PAM"
    C = [-(M-1):2:(M-1);];
    Cr = C;
else
  error("Bad args.");
end

C = C[:]; # Turn into column vector
# for QAM, gam = sqrt(2/3*(M-1))
gam = sqrt(C'*C/length(C));
C = C/gam;
#if nargout ==3
  Cr = Cr/gam;
#end

d2 = real(C[2:end]-C[1]);
id = find(x->abs(x)>0,d2);
d = minimum(abs(d2[id]));
return C, d, Cr,gam
end

