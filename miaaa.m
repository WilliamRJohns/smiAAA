function [bestbcr,bestw,bcr,z,w,fz,err,r] = miaaa(f,Z,tol,normalize,iter)
%Accepts a matrix of function values over Z with each row being a single
%function's values
%Computes Multifunction AAA rational approximations over a common set of poles
%in barycentric and proper rational form such that each functions approximation 
%is with tol in inf-norm sense, where iter = number of Lawson optimization iterations
%and normalize=true, normalizes the functions before computing the approximation and renormalizes after

%Returns the Lawson optimized barycentric approximation bestbcr and its
%associated wieghts bestw. The pre optimization approximation bcr and its
%weights w and the support points z and function values at the support
%poitns fz.


mmax=100;         %Max number of Support Points
k=size(f,1);     %The number of functions

if(normalize)    %Nomalize Each function
    norms=vecnorm(f,2,2);
    
    f=f(1:k,:)./norms(1:k);
end

%Initilize Variables
z=[];           %Support Points
fz=[];          %functions Values at support points
M=length(Z);    
bcr=mean(f,2);  %Initial "barycentric" "rational" 'approximation' for selection of first support point
J=1:M;          %Index vector of Samples not chosen as support points
Jz=[];          %Index vector of the Support Points
C=[];           %Cauchy Matrix
sm={};          %Spares Matrices for Loewner Matrix
err=[];         %Max error at each iteration
for i=1:k
    sm{i}=spdiags(f(i,:).',0,M,M);
end

%Primary Iteration to compute the barycentric rational approximations
for n=1:mmax+1
    %Compute the deviations at each points in Z
     dev=(abs(f-bcr));
     err=[err max(dev,[],'all')];
    %William looking for monotonicity in error problem
   % fprintf('Error with %d support points = %d\n',n-1,max(sum(dev,1)))
    
%     figure()
%     subplot(2,2,1)
%     hold on
%     plot(Z,f(1,:),'b-',Z,bcr(1,:),'r-o')
%     title(sprintf('%d Support points',length(z)))
%     if(~isempty(z))
%         scatter(z,fz(1,:),'cx')
%     end
%     subplot(2,2,2)
%     hold on
%     if(~isempty(z))
%         [pol_pr,res_pr,pra,prhandle]=properrational(z',w,w,fz',bcr,Z);
%         plot(Z,pra(1,:),'c-x',Z,bcr(1,:));
%         legend('pra1','bcr1')
%         disp('Poles');
%         disp(pol_pr);
%     end
%     
%     
%     subplot(2,2,3);
%     hold on
%     plot(Z,f(2,:),'b-',Z,bcr(2,:),'r-o')
%     if(~isempty(z))
%         scatter(z,fz(2,:),'cx')
%     end
%     subplot(2,2,4)
%     hold on
%      if(~isempty(z))
%        plot(Z,pra(2,:),'c-x',Z,bcr(2,:));
%         legend('pra2','bcr2')
%     end
%     hold off
    
    if(max(dev)<tol)
        %disp('Under Tolerance')
        break
    end
    %Find the next support points
    totaldev=sum(dev,1);    
    [~,I]=max(totaldev);    %Find the Point with the largest sum of deviations
    z=[z Z(I)];             %Update the Support Points
    fz=[fz f(:,I)];         %Update the function values at Support Points
    J(J==I)=[];             %Update the Sample point indices 
    Jz=[Jz I];              %Update the Support Point idices(used later in Lawson)
    
    %Build Loewner Matrices for each function
    C = [C 1./(Z.' - Z(I))]; %Next row of Cauchy Matrix
    L=[];
    for i=1:k
        Li=sm{i}*C - C*diag(fz(i,:)); 
        L=[L;Li(J,:)];
    end
    
    [~, S, W] = svd(L, 0);      %Solve for the common weights
    w=W(:,n);
    
    bcr = f;                    %Sets approximation values = functions values at support points
    D = C*w;                    %Calculate the Denominator
    for i=1:k
        Nk = C*(w.*fz(i,:).');  %Calculate the Numerator
        bcr(i,J)=Nk(J)./D(J);   %Sets approximation values at non support points
    end  
%% For Williams powerpoint, will be deleted    
%     figure()
%     
%     loglog(Z./(2*pi*1i),abs(f(1,:)),Z./(2*pi*1i),abs(bcr(1,:)),'LineWidth',2.0)
%     
%     ylim([10^-4,10^1])
%     title(sprintf('f vs aaaf iternation %d',n))
%     hold on
%     freqz=z./(2*pi*1i);
%     scatter(z./(2*pi*1i),abs(fz),'cx','LineWidth',3.0)
%     hold off
%     legend('f','aaaf','support points')
    
    
end
wj=w; % Store unoptimized weights
%Lawson Optimization (iterative wighted least squares)
if(iter>=0)
%Initialize Variables of Lawson Iteration
m=length(z);
gama=1;         %Initial Lawson Exponent
lw=ones(1,M);   %Initial Lawson Weigths
lw=lw/norm(lw);
bestbcr=bcr;    %Best Approximation so far
bestw=repmat(w,2,1);        %w associated with best approximation;
lbcr=[];        %Approximation at current Lawson iteration
maxerror=max(abs(f-bestbcr),[],'all');    %Max abs error of best approx

L=[];           %Build Non-interpolatory Loewner Matrix
%replacement for the interpolation condition r(ti)~f(ti)
Lsupp=[diag(ones(m,1)) -1*diag(ones(m,1))];
for i=1:k
    Li=[sm{i}*C -C*diag(fz(i,:))];
    Li(Jz,:)=Lsupp;
    L=[L;Li];   
end

%figure(); hold on; title('Abs Error of Lawson Iterations'); legend();h=[];
lerr=[];
for l=1:iter %Lawson Iterations
lws=repmat(lw,1,k);
d=spdiags(sqrt(lws).',0,M*k,M*k);  %This should be sparse
%Get the Lawson weighted approximation weight vector
[~, ~, W] = svd(d*L,0);
w=W(:,end);
D = C(J,:)*w(1:m);                                  %Calculate the Denominator
    for i=1:k
        Nk = C(J,:)*(w(m+1:end).*fz(i,:).');        %Calculate the Numerator for each funtion
        lbcr(i,J)=Nk./D;                            %Set approximation values at non support points
        lbcr(i,Jz)=w(m+1:end)'.*f(i,Jz)./w(1:m)';   %Set approximation values at support points
    end
    
%figures for comparing errors
%figure()
%hi=plot(Z./(2*pi*1i),sum(abs(lbcr-f),1))
%h=[h hi];
%legend('Privious','Lawsoned Bcr')

%Compare with the old approximation
lmaxerror=max(abs(f-lbcr),[],'all');
% figure()
% loglog(Z,abs(f),'b',Z,abs(lbcr),'r')
% legend('f','lbcr')
lerr=[lerr ;max(abs(f-lbcr),[],'all')];
if(maxerror<=lmaxerror)
    %gama=gama;
    %disp('changed gama')
else
    fprintf('Optimized from %d to %d\n',maxerror,lmaxerror);
    maxerror=lmaxerror;
    bestbcr=lbcr;
    bestw=w;
end

%Update the Lawson wieghts(I extended to multiple functions with a summation)
testlw=lw.*((max(abs(f-lbcr),[],1)).^gama);
testlw(find(testlw==Inf))=max(testlw(testlw~=Inf));
testlw(isnan(testlw))=mean(testlw(~isnan(testlw)));
if(~isempty(find(testlw==NaN)) | ~isempty(find(testlw==inf)))    %If all values were NaN or Inf they can't be fixed
    disp('Lawson terminated,Weights could not be fixed')
    break;
end
testlw(find(testlw==0))=mean(testlw(testlw~=0));   %avoid any 0's from perfect interpolation
testlw=testlw/norm(testlw);
if(norm(testlw-lw,inf)<1e-8) %This Tolerance should be tested
    break;
end
lw=testlw;
    
end
%figure()
%plot(Z./(2*pi*1i),lerr(1,:),Z./(2*pi*1i),lerr(3,:),'--',Z./(2*pi*1i),lerr(5,:),'--',Z./(2*pi*1i),lerr(end,:),'Linewidth',1.7)
%yline(max(lerr(1,:)),'b','Linewidth',1.5)
%yline(max(lerr(end,:)),'Color','#7E2F8E','Linewidth',1.5)
%legend('MIAAA','Iteration 3','Iteration 5','Final Iteration')
%title('Abs Laswon Iteration Errors')
end

%Remove any support points with zero weight
io=find(w(1:m)==0);
io2=find(w(m+1:end)==0);
io=intersect(io,io2);
z(io)=[];
fz(:,io)=[];
w(io)=[];w(io+m-1)=[];
if(~isempty(io))
    disp('zero wieght removed')
end


if(normalize)    %Re-nomalize Each function
    bcr=bcr(1:k,:).*norms(1:k);
    bestbcr=bestbcr(1:k,:).*norms(1:k);
    fz(1:k,:)=fz(1:k,:).*norms(1:k);
end
%uses non lawson weights
r=@(zz)reval(zz,z.',fz,wj);
end