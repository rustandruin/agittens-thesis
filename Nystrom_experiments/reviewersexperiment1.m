function reviewerexperiment1()

n=201; 
m=21;
numTrials=60;
eps=1e-14; % smallest singular value
rho0=1e-15; % smallest rho
rho1=1; % largest rho

l=100; % pick 100 columns
numrhos=30;
randomSign=false;
doNystrom=true; % if true, pick indR=indC

fprintf('Varying regularization parameter\n');
fprintf('n=%d l=%d\n',n,l);
fprintf('numTrials=%d\n',numTrials);

n2=(n-1)/2; 
m2=(m-1)/2;

x=(0:n-1)/n; 
xi=-n2:n2;

[tmp1,tmp2]=ndgrid(x,xi);

X=exp(2*pi*1i*tmp1.*tmp2)/sqrt(n); % Fourier matrix

rhoList=exp(linspace(log(rho0),log(rho1),numrhos));
rhoList=rhoList(:);

sigma=ones(n,1)*eps;
sigma(1:m)=exp(linspace(log(1),log(eps),m));
sigma=sigma(randperm(n));
if randomSign
    sigma=sigma.*sign(rand(size(sigma))*2-1);
end

data=zeros(3,numTrials,length(rhoList));

for j=1:numTrials
    A=X*diag(sigma)*X';
    indC=randperm(n,l);
    if doNystrom
        indR=indC;
    else
        indR=randperm(n,l);
    end
    
    AC=A(:,indC);
    AR=A(indR,:);
    ARC=A(indR,indC);
    
    for i=1:length(rhoList)
        rho=rhoList(i); % U is coupling matrix
        U=pinv(ARC,rho); % svd truncate
        data(1,j,i)=norm(A-AC*U*AR);
        U=pinv(ARC+rho*eye(l),0); % tikhonov
        data(2,j,i)=norm(A-AC*U*AR);
        data(3,j,i)=norm(A-AC*pinv(ARC,0)*AR);
    end
    
    if mod(j,10)==0
        fprintf('numTrials done=%d\n',j);
        data2=reshape(mean(data(:,1:j,:),2),[3,numrhos]);
        
        clf; 
        hold on;
        plot(rhoList,data2(1,:),'ro');
        plot(rhoList,data2(2,:),'bx');
        plot(rhoList,data2(3,:),'m+');
        set(gca,'yscale','log');
        set(gca,'xscale','log');
        title(sprintf('trials=%d',j));
        xlabel('nrho');
        ylabel('Spectral error');
        
        legend('SVDtruncate','Tikhonov','Simple');
        
        drawnow;        
    end
end
end
