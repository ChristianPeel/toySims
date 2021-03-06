function mimoUPtoy(Ns::Int64,Mc::Int64,M::Int64,K::Int64,
                Ki::Int64,Tt::Int64,Td::Int64,
                rho::Array{Float64,1},gamma=Inf,delta=0)
#mimoUPtoy(Ns,Mc,M,K,Tt,Td,rho)
#  Simulate coding of QAM signals over a flat-fading Gaussian uplink
#  channel, for a base station with multiple receive antennas. Constelations
#  of size Mc are used with M receive antennas and K txmitters, each with
#  one transmit antenna. A power of "rho" per user is utilized. Tt is the
#  number of training samples to use, Td is the number of data symbols, with
#  Tt+Td being the training period. The resulting symbol error rate is
#  averaged over all users and Ns training periods.
#mimoUPtoy(Ns,Mc,M,K,Tt,Td,rho,gamma=Inf,delta=0)
#  If gamma is present, it is a CIR. If delta is present, it is
#  the pilot boost in db.
#
# Examples show SER
#  # Performance as a function of SNR
#  mimoUPtoy(400,4,4,2,0,12,12,[-5.0:5:20])
#  # Performance as a function of CIR
#  mimoUPtoy(400,4,4,2,1,12,12,15,[-15.0:5:25],0)
#  # Performance as a function of Tt
#  mimoUPtoy(400,4,4,2,0,[2:2:16],12,10)
#  # Performance as a function of the number of users
#  mimoUPtoy(400,4,4,[1:4],0,12,12,10)
#  # Performance as a function of the number of antennas
#  mimoUPtoy(1000,4,[2:8],2,0,4,8,10)
#  # Performance as a function of SNR
#  mimoUPtoy(1000,4,4,2,0,4,8,[-5:5:20])
#  # Performance as a function of CIR
#  mimoUPtoy(1000,4,4,2,1,4,8,10,[-15:5:25],0)
# Add an additional 'gamma' argument to the end of any of these functions to
# simulate with a single interferer of CIR gamma. 
#
# By Christian Peel  (chris.peel@ieee.org)
# Last Modified: Fri 14 Nov 14, 11:57am by cpeel

# using PyPlot

println("   Mc     M     K    Ki    Tt    Td   rho  gamma")

# Not using the following 'getAlgs()' method because it doesn't allow
# me to easily plot genie-aided receivers.

algs = getAlgs()
# It would be great to print out the alg names here
# for ax = 1:min(length(algs),4)
#     @printf("%10s ",algs[ax][2])
# end
# println

#algs = []


out = cell(0)
if length(rho)>1
  for s = 1:length(rho)
    push!(out,mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho[s],gamma,delta,algs));
  end
  xval = rho;
  xlab = L"\rho"* " (SNR in dB)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,CIR=%4.1f",
                 Mc,Ns,M,K,Ki,Tt,Td,gamma);
elseif length(gamma)>1
  for s = 1:length(gamma)
    push!(out,mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho,gamma[s],delta,algs));
  end
  xval = gamma;
  xlab = "\gamma (CIR in dB)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f",
                 Mc,Ns,M,K,Ki,Tt,Td,rho);
elseif length(K)>1
  for s = 1:length(K)
    push!(out,mimosimu(Ns,Mc,M,K[s],Ki,Tt,Td,rho,gamma,delta,algs));
  end
  xval = K;
  xlab = "K (# users)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f",
                 Mc,Ns,M,Ki,Tt,Td,rho);
elseif length(Tt)>1
  for s = 1:length(Tt)
    push!(out,mimosimu(Ns,Mc,M,K,Ki,Tt[s],Td,rho,gamma,delta,algs));
  end
  xval = Tt;
  xlab = "Tt (training length)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,M,K,Ki,Td,rho,gamma);
elseif length(Mc)>1
  for s = 1:length(Mc)
    push!(out,mimosimu(Ns,Mc[s],M,K,Ki,Tt,Td,rho,gamma,delta,algs));
  end
  xval = log2(Mc);
  xlab = "# constellation bits";
  tstr = "";
elseif length(M)>1
  for s = 1:length(M)
    push!(out,mimosimu(Ns,Mc,M[s],K,Ki,Tt,Td,rho,gamma,delta,algs));
  end
  xval = M;
  xlab = "M (# Rx antennas)";
  tstr = @sprintf("%dQAM,Ns=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,K,Ki,Tt,Td,rho,gamma);
elseif length(K)==length(M) && length(K)>1
  for s = 1:length(K)
    push!(out,mimosimu(NsMc,M[s],K[s],Ki,Tt,Td,rho,gamma,delta,algs));
  end
  xval = K;
  xlab = "M==K (# users, # antennas)";
  tstr = @sprintf("%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f",
                 Mc,Ns,M,Ki,Tt,Td,rho,gamma);
else
  out = mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho,gamma,delta,algs);
end

pColor = {"r>-","bo--","kx-.","gd-","c^--","m*-.",
          "rs--","gp-.","bv-","kh--","c+-.","m.-",};
pIdx   = 1;

Ns = size(out,1);

#clf()
ser = zeros(Ns);
for a=1:length(out[1][2])
    for k=1:Ns
        ser[k] = out[k][1][a];
#        mi(k) = out(k).mi(a);
    end
#    figure(1)

    semilogy(xval,ser,pColor[pIdx],label=out[1][2][a]);
    hold(true);
##     figure(2)
##     plot(xval,mi,pColor[pIdx]);  hold(true);
    if pIdx == length(pColor)
         pIdx = 1;
    else
         pIdx = pIdx + 1;
    end
end

hold(false);
##figure(1)
xlabel(xlab);
ylabel("SER");
legend();
grid(which="both",axis="y")
grid(which="major",axis="x")
title(tstr)

return
end


#######################################################################
function getAlgs()
#
algsOut = cell(0)

  function zf(Yt,St,Yd,C,onesKT)
    name = "ZF";
    Hhat = Yt*St'*inv(St*St');
    Wls = Hhat/(Hhat'*Hhat);
    Zd = mimo_slice(Wls'*Yd,onesKT,C);
    return Zd,name
  end
  push!(algsOut,zf)

  function ic(Yt,St,Yd,C,onesKT) 
    name = "IC";
    (M,Tt) = size(Yt)
    (K,Tt) = size(St)
    Hhat = (Yt*St')/(St*St');
    Wp = Yt-Hhat*St;   # To be used to estimate the noise variance
    Rhat = Wp*Wp'/(max(Tt-K,1)*M);
    Ri = inv(Rhat);
    W = Ri*Hhat/(Hhat'*Ri*Hhat);
    Zd = mimo_slice(W'*Yd,onesKT,C);
    return Zd,name
  end
  push!(algsOut,ic)

return algsOut
end
#######################################################################


#######################################################################
function mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rdb,gdb,ddb,algs)

delta   = 10.^(ddb/10); # pilot boost
rho     = 10.^(rdb/10);
gamma   = 10.^(gdb/10);
if K>M
    error("Can""t have K>M!!")
end
T = Tt+Td;
Ka = K+Ki;
#
# Modulation class, modulated bits
#
Mp = round(sqrt(Mc));
# Information bits, signals, Threshold.
Zreal = floor(rand(K+Ki,Td*Ns)*Mp);
Zimag = floor(rand(K+Ki,Td*Ns)*Mp);
Zi = Zreal*Mp+Zimag;
ZZ = Zi[1:K,:];
(C,d,Cr) = scodes(Mc,"QAM");
gam = sqrt(2/3*(Mc-1));
C = C[:].';

#
# Channel coefs, noise.
#
sigmaSq = 1;
rho     = rho;
# sigmaSq = 1/rho;
# rho     = 1;
HH = (randn(M,K+Ki,Ns) + im*randn(M,K+Ki,Ns))*sqrt(1/2);# Channel
HH[:,1:K,:] = HH[:,1:K,:]*sqrt(rho);                    # Signal power
HH[:,K+1:Ka,:] = HH[:,K+1:Ka,:]*sqrt(rho/gamma);        # Interference power
NN = (randn(M,T*Ns) + im*randn(M,T*Ns))*sqrt(sigmaSq/2);

#
# Training signal
#
St = repmat(eye(K),1,integer(ceil(Tt/K)));
St = St[:,1:Tt]*sqrt(delta)
SSi = sign(rand(Ki,Tt,Ns)-.5)/sqrt(Tt);

onesKT = ones(K,Td);
#eye2KTd = repmat(eye(2*K),[1 1 Td]);
#eye2KaTd = repmat(eye(2*Ka),[1 1 Td]);
ri = [1:K;]';
ii = [(K+1):(2*K);]';
ria = [1:K;]';
iia = [(Ka+1):(Ka+K);]';


Nalgs = length(algs);
#Nalgs = 2;
ZD = zeros(K,Td*Ns,Nalgs);
algNames = cell(Nalgs)

for ix = 1:Ns
    zix = (ix-1)*Td+1:ix*Td;
    N = NN[:,(ix-1)*T+1:ix*T];
    Z = Zi[:,zix];
    Si = SSi[:,:,ix];
    Sti = [St; Si];
    Sd = C[Z+1];
#    Sd = fft(Sd,[],2)/sqrt(Td);
    S = [Sti Sd];
    H = HH[:,:,ix];
    
    Y = H*S+N; # Channel
    Yt = Y[:,1:Tt];
    Yd = Y[:,Tt+1:Tt+Td];
    Ydc = [real(Yd); imag(Yd)];
    #
    # Genie channel
    #
    Hdes = H[:,1:K];
    Hc = [real(Hdes) -imag(Hdes); imag(Hdes) real(Hdes)];
    HcAll = [real(H) -imag(H); imag(H) real(H)];
    #
    # Training
    #
    Hhat = Yt*St'*inv(St*St');
    Hhatc = [real(Hhat) -imag(Hhat); imag(Hhat) real(Hhat)];
    
    #
    # Receiver algorithms
    #
    for ax = 1:Nalgs
        (Zd,name) = algs[ax](Yt,St,Yd,C,onesKT);
        ZD[:,zix,ax] = Zd;
        algNames[ax] = name;
    end


    # ax=0;

    # # Showing some toy example receivers

    # ax = ax+1;
    # algNames[ax] = "ZF";
    # Hhat = Yt*St'*inv(St*St');
    # Wls = Hhat/(Hhat'*Hhat);
    # Zd  = mimo_slice(Wls'*Yd,onesKT,C);
    # ZD[:,zix,ax] = Zd;

    # ax = ax+1;
    # algNames[ax] = "IC";
    # Hhat = (Yt*St')/(St*St');
    # Wp = Yt-Hhat*St;   # To be used to estimate the noise variance
    # Rhat = Wp*Wp'/(max(Tt-K,1)*M);
    # Ri = inv(Rhat);
    # W = Ri*Hhat/(Hhat'*Ri*Hhat);
    # Zd = mimo_slice(W'*Yd,W'*Hhat*onesKT,C);
    # ZD[:,zix,ax] = Zd;


end
@printf("%5d %5d %5d %5d %5d %5d %5.1f %5.1f",
          Mc,  M,  K, Ki, Tt, Td,  rdb,  gdb)

ser = zeros(Nalgs);
#Calculate symbol error rate
for ax = 1:Nalgs
    ser[ax] = mean(map(!=,ZZ,ZD[:,:,ax]));
    if ax<6
        @printf("    %7.4f",ser[ax])
    end 
#
end
@printf("\n")

return ser, algNames
end

#######################################################################
function mimo_slice(Y,H,C)
M = size(Y,1);
Td = size(Y,2);
Zd = zeros(M,Td);
Mc = length(C);
dy = zeros(M,Mc);
for t=1:Td
    for i = 1:Mc
        dd = Y[:,t]-H[:,t]*C[i];
        dy[:,i] = abs(dd).^2
    end
    for m = 1:M
        (mn,zd) = findmin(dy[m,:],2);
        Zd[m,t] = zd[1]-1.0;
    end
end
return Zd
end
#######################################################################


