clc
clear all
initialprob=.65;
r=[.02 .07];
sigma=[.52 1.86];
beta=.4;
s0=50;
T=1;
n=300;
M=500;
lambda=5;
dt=T/n;
t=[0 cumsum(ones(1,n-1)*dt)];
Q=[-.20 .20;.6 -.6];
strike=30;
K=10;
flag=1;
mu=zeros(K,1);
X=ones(K,K);
tic
while flag==1;
    X=-abs(rand(K,K)./K);
    for i=1:K
        X(i,i)=1;
        X(i+1:K,i)=X(i,i+1:K);
    end
    if min(eig(X))>0
        flag=0;
    end
end
R=mvnrnd(mu,X,n*M);
T1=toc;
s_N=zeros(n,M);
s_N(1,:)=s0;s_KAVE=s_N;
tic
for jj=1:M
    pd = makedist('Binomial','N',1,'p',1-initialprob);
    ber = random(pd,1,1);
    sim=7;
    while sim>6
        sim=poissrnd(lambda);
    end
    flag=0;
    while flag==0
        for j=1:sim
            if mod(j,2)==0
                k=1;
            else
                k=2;
            end
            u(j)=exprnd(-Q(k,k));
        end
        L=[0 cumsum(u) T];
        for i=1:sim+1
            L1(i)=L(i+1)-L(i);
        end
        if min(L1)<dt
            flag=0;
        else
            flag=1;
        end
    end
    for i=1:sim+1
        if ber==0
            k=abs(((-1)^(i+1)-3)/2);
        else
            k=abs(((-1)^(i)-3)/2);
        end
        sigma1(i)=sigma(k);
        r1(i)=r(k);
    end
    for j=2:n
        t1=t(j);
        for i=1:sim+1
            if L(i)<t1&&t1<L(i+1)
                sigma_zt(j-1)=sigma1(i);
                r_zt(j-1)=r1(i);
                break;
            end
        end
        s_N(j,jj)=s_N(j-1,jj)+r_zt(j-1)*s_N(j-1,jj)*dt+sigma_zt(j-1)*s_N(j-1,jj)^beta*sqrt(dt)*randn;
    end
    for i=2:n-1
        r_N(i-1)=((r_zt(i)+r_zt(i-1))*dt)/2;
    end
    rr_N(jj)=exp(-sum(r_N));
end
T2=toc;
s_kave=[];
tic
for kk=1:K
    p=1;R2=R(:,kk);
    for i=1:n
        for j=1:M
            R1(i,j)=R2(p);
            p=p+1;
        end
    end
    for jj=1:M
        pd = makedist('Binomial','N',1,'p',1-initialprob);
        ber = random(pd,1,1);
        sim=7;
        while sim>6
            sim=poissrnd(lambda);
        end
        flag=0;
        while flag==0
            for j=1:sim
                if mod(j,2)==0
                    k=1;
                else
                    k=2;
                end
                u(j)=exprnd(-Q(k,k));
            end
            L=[0 cumsum(u) T];
            for i=1:sim+1
                L1(i)=L(i+1)-L(i);
            end
            if min(L1)<dt
                flag=0;
            else
                flag=1;
            end
        end
        for i=1:sim+1
            if ber==0
                k=abs(((-1)^(i+1)-3)/2);
            else
                k=abs(((-1)^(i)-3)/2);
            end
            sigma1(i)=sigma(k);
            r1(i)=r(k);
        end
        for j=2:n
            t1=t(j);
            for i=1:sim+1
                if L(i)<t1&&t1<L(i+1)
                    sigma_zt(j-1)=sigma1(i);
                    r_zt(j-1)=r1(i);
                    break;
                end
            end
            s_KAVE(j,jj)=s_KAVE(j-1,jj)+r_zt(j-1)*s_KAVE(j-1,jj)*dt+sigma_zt(j-1)*s_KAVE(j-1,jj)^beta*sqrt(dt)*R1(j-1,jj);
        end
        
        for i=2:n-1
            r_KAVE(i-1)=((r_zt(i)+r_zt(i-1))*dt)/2;
        end
        rr_KAVE(jj,kk)=exp(-sum(r_KAVE));
    end
    s_kave=[s_KAVE s_kave];
end
T3=toc;
for i=1:M
    pay_N(i)=rr_N(i)*max(mean(s_N(:,i))-strike,0);
end
g=1;
for k=1:K
    s_kave1=s_kave(1:n,g:k*M);
    for i=1:M
        pay_KAVE1(i,k)=rr_KAVE(i,k)*(max(mean(s_kave1(:,i))-strike,0));
    end
    g=k*M+1;
end
for i=1:M
    pay_KAVE(i)=1/K*sum(pay_KAVE1(i,:));
end
Mean_Naive=mean(pay_N)
STD_Naive=std(pay_N)
Mean_KAVE=mean(pay_KAVE)
STD_KAVE=std(pay_KAVE)

