%function properrational takes a barycentric rational approximation as
%support points(z),weights(wnum and wden),function values at the support points(fz) over doamin Z 
%Returns the proper rational approximation (poles,residues,proper rational
%approx) and a function handle for each approximation
function [poles,res,bestpra,prhandle,bestpoly]=properrational(z,wnum,wden,fz,bcf,Z)
maxpolydegree=1; %A better way to decide this degree probably exists

k=size(bcf,1);

%Compute the poles with new and improved prz routine
%[ta, poles] = zpf(z, wden, 1e-10, 1e-10);
%fprintf('Final Pole accuracy = %d\n',ta);


%Use poles from deflated prz e-value problem
%poles2=przl(z,wden);
poles=przd(z,wden);

%Compute the residues via Cauchy Matrices
Cnum=bsxfun(@minus,poles,z.').^(-1);
Cden=-1*Cnum.^2;


res=((Cnum*(wnum.*fz))./(Cden*wden));


%Recalculate the partial fraction part of the approximation in proper rational form
CC=bsxfun(@minus,Z,poles).^(-1);
pra=res.'*CC;
%fprintf('Cauchy Matrix Conditions are %d %d %d\n',cond(Cnum),cond(Cden),cond(CC));

%Calculate the polynomial parts by fitting the remaider
remainder=bcf-pra;
bestpra=pra;
bestpoly ={};
for i=1:k
besterr=max(abs(remainder(i,:)));
%fprintf('Pre-Polynomial fitting error = %d\n',besterr);
errvec=[];
if(besterr>=1e-12)
for polydegree = maxpolydegree:maxpolydegree
    
    polypart=polyfit(Z,remainder(i,:),polydegree);
    pra2=pra(i,:)+polyval(polypart,Z);
        
    err=max(abs(bcf(i,:)-pra2));
    errvec=[errvec err];
   % fprintf('Polydegree = %d with err = %d\n',polydegree,err);
    if(err<besterr)
        bestpra(i,:)=pra2;
        besterr=err;
        bestpoly{i} = polypart;
    end
end

%disp('Best error acheived was ');
[G,I]=min(errvec);
%disp(G);
%disp('With a polynomial of degree');
%disp(I-1);
end
end
%Build function handles for each approximation
prhandle={};
for i=1:k
    prhandle{i} = @(ww) pfeval(ww,poles,res(:,i),bestpoly{i});
end
end