function mimoUPtoy(Ns,Mc,M,K,Ki,Tt,Td,rho,gamma,delta)
%mimoUPtoy(Ns,Mc,M,K,Tt,Td,rho)
%  Simulate coding of QAM signals over a flat-fading Gaussian uplink
%  channel, for a base station with multiple receive antennas. Constelations
%  of size Mc are used with M receive antennas and K txmitters, each with
%  one transmit antenna. A power of "rho" per user is utilized. Tt is the
%  number of training samples to use, Td is the number of data symbols, with
%  Tt+Td being the training period. The resulting symbol error rate is
%  averaged over all users and Ns training periods.
%mimoUPtoy(Ns,Mc,M,K,Tt,Td,rho,gamma,delta)
%  If gamma is present, it is a CIR. If delta is present, it is
%  the pilot boost in db.
%
% Examples show SER
%  % Performance as a function of SNR
%  mimoUPtoy(400,4,4,2,0,12,144,[-5:5:20])
%  % Performance as a function of CIR
%  mimoUPtoy(400,4,4,2,1,12,144,15,[-15:5:25],0)
%  % Performance as a function of Tt
%  mimoUPtoy(400,4,4,2,0,[2:2:16],144,10)
%  % Performance as a function of the number of users
%  mimoUPtoy(400,4,4,[1:4],0,12,144,10)
%  % Performance as a function of the number of antennas
%  mimoUPtoy(1000,4,[2:8],2,0,4,8,10)
%  % Performance as a function of SNR
%  mimoUPtoy(1000,4,4,2,0,4,8,[-5:5:20])
%  % Performance as a function of CIR
%  mimoUPtoy(1000,4,4,2,1,4,8,10,[-15:5:25],0)
% Add an additional 'gamma' argument to the end of any of these functions to
% simulate with a single interferer of CIR gamma. 
%
% By Christian Peel  (chris.peel@ieee.org)
% Last Modified: Fri 14 Nov 14, 11:57am by cpeel

if nargin==0
  error('Bad arguments. Try again.'); 
end
if nargin<=8
    gamma =inf;    % Turn interferer off
end
if nargin<=9
    delta =0;
end

disp('    Mc     M     K    Ki    Tt    Td   rho  gamma')

if length(rho)>1
  for s = 1:length(rho)
    algs(s) = mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho(s),gamma,delta);
  end
  xval = rho;
  xlab = '\rho (SNR in dB)';
  tstr = sprintf('%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,CIR=%4.1f',...
                 Mc,Ns,M,K,Ki,Tt,Td,gamma);
elseif length(gamma)>1
  for s = 1:length(gamma)
    algs(s) = mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho,gamma(s),delta);
  end
  xval = gamma;
  xlab = '\gamma (CIR in dB)';
  tstr = sprintf('%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f',...
                 Mc,Ns,M,K,Ki,Tt,Td,rho);
elseif length(K)>1
  for s = 1:length(K)
    algs(s) = mimosimu(Ns,Mc,M,K(s),Ki,Tt,Td,rho,gamma,delta);
  end
  xval = K;
  xlab = 'K (# users)';
  tstr = sprintf('%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f',...
                 Mc,Ns,M,Ki,Tt,Td,rho);
elseif length(Tt)>1
  for s = 1:length(Tt)
    algs(s) = mimosimu(Ns,Mc,M,K,Ki,Tt(s),Td,rho,gamma,delta);
  end
  xval = Tt;
  xlab = 'Tt (training length)';
  tstr = sprintf('%dQAM,Ns=%d,M=%d,K=%d,Ki=%d,Td=%d,SNR=%4.1f,CIR=%4.1f',...
                 Mc,Ns,M,K,Ki,Td,rho,gamma);
elseif length(Mc)>1
  for s = 1:length(Mc)
    algs(s) = mimosimu(Ns,Mc(s),M,K,Ki,Tt,Td,rho,gamma,delta);
  end
  xval = log2(Mc);
  xlab = '# constellation bits';
  tstr = '';
elseif length(M)>1
  for s = 1:length(M)
    algs(s) = mimosimu(Ns,Mc,M(s),K,Ki,Tt,Td,rho,gamma,delta);
  end
  xval = M;
  xlab = 'M (# Rx antennas)';
  tstr = sprintf('%dQAM,Ns=%d,K=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f',...
                 Mc,Ns,K,Ki,Tt,Td,rho,gamma);
elseif length(K)==length(M) && length(K)>1
  for s = 1:length(K)
    algs(s) = mimosimu(NsMc,M(s),K(s),Ki,Tt,Td,rho,gamma,delta);
  end
  xval = K;
  xlab = 'M==K (# users, # antennas)';
  tstr = sprintf('%dQAM,Ns=%d,M=%d,Ki=%d,Tt=%d,Td=%d,SNR=%4.1f,CIR=%4.1f',...
                 Mc,Ns,M,Ki,Tt,Td,rho,gamma);
else
  algs = mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rho,gamma,delta);
end

pColor = {'r>-','bo--','kx-.','gd-','c^--','m*-.',...
          'rs--','gp-.','bv-','kh--','c+-.','m.-',};
pIdx   = 1;
% figure(1); clf
% figure(2); clf

for a=1:length(algs(1).name)
    for k=1:length(algs)
        ser(k) = algs(k).ser(a);
    end
    semilogy(xval,ser,pColor{pIdx});  hold on;
    if   pIdx == length(pColor), pIdx = 1;
    else pIdx = pIdx + 1;
    end
end
hold off;
%figure(1)
xlabel(xlab);
ylabel('SER');
legend(algs(1).name,4);
grid on;
title(tstr)

disp(sprintf('\a'))
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function algs = mimosimu(Ns,Mc,M,K,Ki,Tt,Td,rdb,gdb,ddb)

delta   = 10.^(ddb/10); % pilot boost
rho     = 10.^(rdb/10);
gamma   = 10.^(gdb/10);
if K>M
    error('Can''t have K>M!!')
end
T = Tt+Td;
Ka = K+Ki;
%
% Modulation class, modulated bits
%
Mp = round(sqrt(Mc));
% Information bits, signals, Threshold.
Zreal = floor(rand(K+Ki,Td*Ns)*Mp);
Zimag = floor(rand(K+Ki,Td*Ns)*Mp);
Zi = Zreal*Mp+Zimag;
ZZ = Zi(1:K,:);
[C d Cr] = scodes(Mc,'QAM');
gam = sqrt(2/3*(Mc-1));
C = C(:).';

%
% Channel coefs, noise.
%
sigmaSq = 1;
rho     = rho;
% sigmaSq = 1/rho;
% rho     = 1;
HH = (randn(M,K+Ki,Ns) + j*randn(M,K+Ki,Ns))*sqrt(1/2); % Signal Channel
HH(:,1:K,:) = HH(:,1:K,:)*sqrt(rho);                    % Signal power
HH(:,K+1:Ka,:) = HH(:,K+1:Ka,:)*sqrt(rho/gamma);        % Interference power
NN = (randn(M,T*Ns) + j*randn(M,T*Ns))*sqrt(sigmaSq/2);

%
% Training signal
%
St = repmat(eye(K),[1 ceil(Tt/K)]);
St = St(:,1:Tt)*sqrt(delta);
SSi = sign(rand(Ki,Tt,Ns)-.5)/sqrt(Tt);

onesKT = ones(K,Td);
eye2KTd = repmat(eye(2*K),[1 1 Td]);
eye2KaTd = repmat(eye(2*Ka),[1 1 Td]);
ri = [1:K]';
ii = [(K+1):(2*K)]';
ria = [1:K]';
iia = [(Ka+1):(Ka+K)]';

algs = struct('name',[]);
for ix = 1:Ns
    zix = (ix-1)*Td+1:ix*Td;
    N = NN(:,(ix-1)*T+1:ix*T);
    Z = Zi(:,zix);
    Si = SSi(:,:,ix);
    Sti = [St; Si];
    Sd = C(Z+1);
%    Sd = fft(Sd,[],2)/sqrt(Td);
    S = [Sti Sd];
    H = HH(:,:,ix);
    
    Y = H*S+N; % Channel
    Yt = Y(:,1:Tt);
    Yd = Y(:,Tt+1:Tt+Td);
    Ydc = [real(Yd); imag(Yd)];
    %
    % Genie channel
    %
    Hdes = H(:,1:K);
    Hc = [real(Hdes) -imag(Hdes); imag(Hdes) real(Hdes)];
    HcAll = [real(H) -imag(H); imag(H) real(H)];
    %
    % Training
    %
    Hhat = Yt*St'*inv(St*St');
    Hhatc = [real(Hhat) -imag(Hhat); imag(Hhat) real(Hhat)];
    
    %
    % Receiver algorithms
    %
    a=0;

    % Showing some toy example receivers

    a = a+1;
    algs.name{a} = 'ZF';
    Hhat = Yt*St'*inv(St*St');
    Wls = Hhat/(Hhat'*Hhat);
    Zd = mimo_slice(Wls'*Yd,onesKT,C);
    algs.Zd{a}(:,zix) = Zd;

    a = a+1;
    algs.name{a} = 'IC';
    Hhat = (Yt*St')/(St*St');
    Wp = Yt-Hhat*St;   % To be used to estimate the noise variance
    Rhat = Wp*Wp'/(max(Tt-K,1)*M);
    Ri = inv(Rhat);
    W = Ri*Hhat/(Hhat'*Ri*Hhat);
    Zd = mimo_slice(W'*Yd,W'*Hhat*onesKT,C);
    algs.Zd{a}(:,zix) = Zd;

end
disp([Mc M K Ki Tt Td round(rdb) round(gdb)])

%Calculate symbol error rate
for a = 1:length(algs.name)
    algs.ser(a) = mean(mean(ZZ~=algs.Zd{a}));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Zd = mimo_slice(Y,H,C)
[M Td] = size(Y);
Zd = zeros(M,Td);
Mc = length(C);
dy = zeros(M,Mc);
for t=1:Td
    for i = 1:Mc
        dy(:,i) = abs(Y(:,t)-H(:,t)*C(i));
    end
    [m zd] = min(dy,[],2);
    Zd(:,t) = zd-1;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,d,Cr] = scodes(M,type)
%C = scodes(M,type)
%  Return a constellation of M signals, where <type> is 'PSK', 'QAM',
%  or 'PAM'.
%[C,d] = scodes(M,type)
%  Also return the minimum distance between constellation points.
%[C,d,Cr] = scodes(M,type)
%  Returns the underlying real PAM constellation for square (M=m^2)
%  QAM constellations
%
% By Christian Peel  (chris.peel@ieee.org)
% Last Modified: Tue 23 Jan 07, 11:02am by cpeel
%
if nargin==0
  error('Bad arguments. Try again.'); 
end

if strcmp(type,'PSK')
  C=exp(j*2*pi*[1:M]'/M);
elseif strcmp(type,'QAM')
  if abs(sqrt(M)-round(sqrt(M)))< 1e-3
    % Make QAM constellations if M is a square
    m = round(sqrt(M));
    Cr = [-(m-1):2:(m-1)];
    for ix = 1:m
      for jx = 1:m
        C(jx,ix) = Cr(ix) +j*Cr(jx);
      end
    end
  else
    error('Can only handle square QAM constellations.')
%     [x y] = qaskenco(M);
%     C = x+j*y;
 end
elseif strcmp(type,'PAM')
  C = [-(M-1):2:(M-1)];
else
  error('Bad args.');
end

C = C(:);
% for QAM, gam = sqrt(2/3*(M-1))
gam = sqrt(C'*C/length(C));
C = C/gam;
if nargout ==3
  Cr = Cr/gam;
end

d2 = real(C(2:end)-C(1));
id = find(abs(d2)>0);
d = min(abs(d2(id)));
return

