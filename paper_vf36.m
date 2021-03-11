% ex4c.m             
%
% Fitting all elements of the admittance matrix Y of a six-terminal system 
% (power system distribution network)
%
% -Reading frequency admittance matrix Y(s) from disk.
% -Stacking the (21) elements of the lower triangle of Y into a single
%  column f(s)
% -Fitting g(s)=sum(f(s)) --> new initial poles
%   -Initial poles: 25 linearly spaced complex pairs (N=50)
%   -5 iterations
% -Fitting elements of vector f(s) using a common pole set 
%  (3 iterations)
% -Converting resulting state-space model into a model for the full Y using tri2full.m.
% -Converting state space model (A,B,C,D,E) into a pole-residue model (a,R,D,E) using ss2pr.m
%
% This example script is part of the vector fitting package (VFIT3.zip) 
% Last revised: 08.08.2008. 
% Created by:   Bjorn Gustavsen.
%
clear all


N=55; %order of approximation
Niter1=5; %Fitting column sum: n.o. iterations
Niter2=3; %Fitting column: n.o. iterations

disp('Reading data from file ...') %--> s(1,Ns), bigY(Nc,Nc,Ns)
fid1=fopen('fdne.txt','r');
Nc=fscanf(fid1,'%f',1);
Ns=fscanf(fid1,'%f',1);
bigY=zeros(Nc,Nc,Ns); s=zeros(1,Ns);
for k=1:Ns
  s(k)=fscanf(fid1,'%e',1);
  for row=1:Nc
    for col=1:Nc
      dum1=fscanf(fid1,'%e',1);
      dum2=fscanf(fid1,'%e',1);   
      bigY(row,col,k)=dum1+j*dum2;
    end
  end
end
s=i*s;
fclose(fid1);

tic
disp('-----------------S T A R T--------------------------')

disp('****Stacking matrix elements (lower triangle) into single column....')
tell=0;
for col=1:Nc
  for row=col:Nc
    tell=tell+1;
    f(tell,:)=squeeze(bigY(row,col,:)).'; %stacking elements into a single vector
  end
end


AA=sparse([]); BB=sparse([]); CC=[]; DD=[]; EE=[];
bigf=[]; bigfit=[];


%Complex starting poles :
w=s/i;
bet=linspace(w(1),w(Ns),N/2);
poles=[];
for n=1:length(bet)
  alf=-bet(n)*1e-2;
  poles=[poles (alf-i*bet(n)) (alf+i*bet(n)) ]; 
end


%weight=ones(1,Ns);
%weight=1./abs(f);
weight=1./sqrt(abs(f));

%Fitting options
opts.relax=1;      %Use vector fitting with relaxed non-triviality constraint
opts.stable=1;     %Enforce stable poles
opts.asymp=3;      %Fitting includes D and E (%William Changed)
opts.spy1=0; 
opts.spy2=1; 
opts.logx=0; 
opts.logy=1; 
opts.errplot=1;
opts.phaseplot=1;

opts.skip_pole=0; 
opts.skip_res=1;
opts.cmplx_ss=1;  %=1 --> Will generate state space model with diagonal A
opts.legend=1;

  

%Forming (weighted) column sum:
g=0;
for n=1:Nc
  %g=g+f(n,:); %unweighted sum     
  g=g+f(n,:)/norm(f(n,:));
  %g=g+f(n,:)/sqrt(norm(f(n,:)));     
end
weight_g=1./abs(g);

disp('****Calculating improved initial poles by fitting column sum ...')
for iter=1:Niter1
   disp(['   Iter ' num2str(iter)])
   if iter==Niter1,opts.skip_res=0; end
   [SER,poles,rmserr,fit]=vectfit3(g,s,poles,weight_g,opts);  
end

 

disp(['****Fitting column ...'])
opts.skip_res=1;
for iter=1:Niter2
  disp(['   Iter ' num2str(iter)])
  if iter==Niter2, opts.skip_res=0; end
  [SER,poles,rmserr,fit]=vectfit3(f,s,poles,weight,opts);  
end

disp('vf runtime')
%toc
A=SER.A; B=SER.B; C=SER.C; D=SER.D; E=SER.E; 

%disp('****Transforming model of lower matrix triangle into state-space model of full matrix....')
%[SER]=tri2full(SER);


%disp('****Generating pole-residue model....')
%[R,a]=ss2pr(SER.A,SER.B,SER.C);


disp('-------------------E N D----------------------------')
toc

freq=s/(2*1i*pi);

%fit the data with aaa
%reoganize the data
k=1; F=[];
for ii=1:6
    for j=1:6
        F(k,:)=bigY(ii,j,:);
        k=k+1;
    end
end
%remove dupliate values(Why do these exists?)
[s,ia,~]=unique(s);
F=F(:,ia);


%Approximation of a single function from the power systems example
[aaaf,poles_aaa,~,~] = aaa(f(2,:),s,'tol',1e-6);
disp('AAA Runtime')
mm = @() aaa(f(2,:),s,'tol',1e-6); % handle to function
timeit(mm)


%Stable approximation of the single function
[saaaf,wj,~,zj,wj1,fj,err] = smiaaa(f(2,:),s,1e-6,false,0,1);
wden = wj(1:length(zj)); wnum = wj(length(zj)+1:end);
[spoles_aaa,sres_aaa,pfaaaf,~,bestpoly]=properrational(zj.',wnum,wden,fj.',saaaf,s);
disp('smiAAA Runtime')
mm = @() smiaaa(f(2,:),s,1e-6,false,0,1); % handle to function
timeit(mm)

figure();hold on;
set(gca,'Box','on');
scatter(real(poles_aaa),imag(poles_aaa),'rx');       %poles from aaa
scatter(real(spoles_aaa),imag(spoles_aaa),'bo');     %poles from stableaaa
legend('25 AAA poles','28 smiAAA poles','Location','Northwest')
title('Poles for AAA and smiAAA')

figure();
loglog(freq,abs(f(2,:)-aaaf(s)+1e-13),freq,abs(f(2,:)-saaaf+1e-13))
legend('AAA','smiAAA')
title('AAA vs smiAAA tol=1e-6')

%Stable approximation of multiple functions with common poles can also be
%constructed similar to vector fitting (this does not include the symetry
%of vector fitting)
[saaaf,wj,~,zj,wj1,fj,err] = smiaaa(f,s,1e-6,false,0,1);
wden = wj(1:length(zj)); wnum = wj(length(zj)+1:end);
[spoles_aaa,sres_aaa,pfaaaf,~,bestpoly]=properrational(zj.',wnum,wden,fj.',saaaf,s);
disp('smiAAA 21 functions Runtime')
mm = @() smiaaa(f,s,1e-6,false,0,1); % handle to function
timeit(mm)

figure();
loglog(freq,abs(f-saaaf)+1e-13)
title('smiAAA tol=1e-6 21 functions')

%Fully Symetric Problem for direct comparison with VF
[SER,rmserr,bigYfit,opts2]=VFdriver(bigY,s,poles,opts);
linf_vf=max(abs(bigY-bigYfit),[],'all')

[lpaaaf,pwj,paaaf,pzj,~,pfj] = smiaaasymetric(F,s,1e-3,false,10,1);
nn=length(pwj)/2;
[ppoles_aaa,pres_aaa,lppfaaaf,~,pbestpoly]=properrational(pzj.',pwj(nn+1:end),pwj(1:nn),pfj.',F,s);
%[fnew,aaa_res,D,E]= residueparty(F,s,ppoles_aaa);
rmserr_aaal=comp_error(F,lpaaaf);
rmserr_aaa=comp_error(F,paaaf);
linf_aaa=max(abs(paaaf-F),[],'all');
linf_aaal=max(abs(lpaaaf-F),[],'all');


