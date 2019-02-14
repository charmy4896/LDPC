clc;
clear all;
%message bit of K bit and N total bit and parity bits
K=16;
N=25;
P=N-K;
n=sqrt(N);
%randomly create K bits message bit(2^K possibility)
M = mod(randi(10,K,1),2);
%creating generator matrix by seeing 3*3 matix representation of (9,4) code
%G=[1 0 0 0; 0 1 0 0;1 1 0 0; 0 0 1 0;0 0 0 1;0 0 1 1;1 0 1 0;0 1 0 1;1 1 1 1];
starts =1 ;
messagebits=1;
ends =1;
G=zeros(N,K);
for i=1:N 
        if (mod(i,n)==0)
            ends=ends-1;
            for k=starts:ends
                G(i,k)=1;
            end
            starts=starts+n-1;
             if(messagebits>K)
                 break;
             end
        else
            for j=1:K
                if (j==messagebits)
                    G(i,j)=1;
                else
                    G(i,j)=0;
                end
            end
            messagebits=messagebits+1;
        end
        ends=ends+1;
    
end
x=1;
for l=i+1:N-1
    a=x;
    for m=1:N-n
        
        if(m==x)
            G(l,x)=1;
            x=x+n-1;
        if(x>K)
            break;
        end
        %else
        %  G(l,m)=0;
        end
    end
    x=a+1;
end
for i=1:K
    G(N,i)=1;
end
    
% creating codeword using G*M mod 2
C = transpose(mod((G*M),2));
%received data
p=0.4;
R = bsc(C,p); 
Q=R;
Q=reshape(Q,n,n);
%xor R with C to get erasure bits
errarr = xor(R,C);
%replace erasure with x=2 in R
for i=1:N
    if (errarr(i)==1)
        R(i)=2;
    end
end
%ans=array2table(R);
%matrix(n,n,);
%R=[2,2,1;1,2,0;0,2,2];
Rnot=reshape(R,n,n);
R=reshape(R,n,n);
for i=1:25
    for j=1:n
        %row wise checking
        count=0;
        for k=1:n
            if(R(j,k)==2)
                count=count+1;
                index=k;
            end
        end
        if (count==1)
            xors=0;
            for k=1:n
                if (R(j,k)~=2)
                    xors=xor(R(j,k),xors);
                end
            end
            R(j,index)=xors;
        end
        
        %column wise checking
        count=0;
        for k=1:n
            if(R(k,j)==2)
                count=count+1;
                index=k;
            end
        end
        if (count==1)
            xors=0;
            for k=1:n
                if (R(k,j)~=2)
                    xors=xor(R(k,j),xors);
                end
            end
            R(index,j)=xors;
        end
        
        
        
    end
end
% create parity check matrix for (9,4) code
%H=[1 1 1 0 0 0 0 0 0;0 0 0 1 1 1 0 0 0;1 0 0 1 0 0 1 0 0;0 1 0 0 1 0 0 1 0;0 0 1 0 0 1 0 0 1];
i=1
for j=1:N
        H(i,j)=1      
        if(mod(j,n)==0)
            i=i+1;
        end
end
for i= n:P
    for j=1:N
        if(mod(j,n)==i-n+1)
            H(i,j)=1;
        else
            H(i,j)=0;
        end
        if(mod(j,n)==0 && i==P)
            H(i,j)=1;
        end
    end
end
%creating bipartite graph for H of our (N,K) code 
cn=P; % cn= check nodes
vn=N; % vn= variable nodes
tanner_graph = [zeros(cn,cn), H;H', zeros(vn,vn)];     
h = plot(graph(tanner_graph));
h.YData(1:cn) = 2;
h.YData((cn+1):end) = 1;
h.XData(1:cn) = linspace(0,1,cn);
h.XData((cn+1):end) = linspace(0,1,vn);
