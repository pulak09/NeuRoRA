function [R,RR,I,a] = RandomSO3Graph(N,Completeness,Sigma,Outlier)
% [RR] = function RandomSO3Graph(N)
% INPUT  : N            = Number of nodes or View Points or Observation Points
%          Completeness = Fraction [0,1] of edges that will have Rij edges
%          Sigma        = Noise added to Relative Rotations in degree
%          Outlier      = Fraction of Outliers
% OUTPUT : R            = Absolute Rotations. R(:,:,1)=eye(3);
%          RR,I         = Relative Rotation Matrices and the Indices. RR(:,:,i) is
%                         the relative rotation between the nodes I(1,i) and I(2,i)
%                         RR(:,:,i) = R(:,:,I(2,i))*R(:,:,I(1,i))' 
Sigma=Sigma*pi/180/sqrt(3);
R=zeros(3,3,N);
R(:,:,1)=eye(3); tt = zeros(1, N); 
for i=2:N
    w=randn(3,1);w=w/norm(w)*pi*rand(1)/18; y = 2*(rand(1)-0.5); 
    w2 =  [y, sqrt(1 - y^2), 0.0]; w2 = w2/norm(w2)*pi*randn/4; tt(i) = w2(2); 
    R(:,:,i)= w2R(w2)*w2R(w)'; 
end
[~, ii] = sort(tt); 
R = R(:, :, ii); 
I=zeros(2,N*(N-1)/2);RR=zeros(3,3,N*(N-1)/2); %angle = zeros(1,N*(N-1)/2); 
k=1;
for i=1:N-1
    for j=i+1:min(N, 2*(i+2))
        I(:,k)=[i;j];
        RR(:,:,k)=R(:,:,j)*R(:,:,i)';
%         angle(k) = acos(0.5*(trace(RR(:,:,k)) - 1))*180/pi; 
        k=k+1;
    end
end

% angle = angle(1:k-1); 
% [, id] = sort(angle); 
p = randperm(k-1);
I = I(:, p); 
RR = RR(:, :, p); 

% Outlier=round(size(I,2)*Completeness*Outlier);

if(nargin>1&&Completeness~=1)
    G=zeros(N,N);d=zeros(1,size(I,2));
    k=1;
    for i=1:N-1
        for j=i+1:min(N, 2*(i+2))
            G(j,i)=acos((trace(RR(:,:,k))-1)/2)*180/pi;
            d(1,k)=G(j, i);
            k=k+1;
        end
    end
    [MST, ~]=graphminspantree(sparse(G));
    k= false(1,size(I,2));
    for i=1:size(I,2)
        if(MST(I(2,i),I(1,i))~=0);    k(i)=1;     d(i)=inf;    end;
    end
    e=min(ceil(Completeness*size(I,2)), sum(d < 45+rand*15));
    if(e>N-1);    [~,dd]=sort(d);    k(dd(1:e-N+1))=1;    end;
    I=I(:,k);    RR=RR(:,:,k);
end

r2 = Sigma; out = Outlier; 
if(nargin>2)
    for i=1:size(RR,3)
%         q1 = R2q(R(:, :, I(1, i))); q1 = quat2axang(q1);         
%         q2 = R2q(R(:, :, I(2, i))); q2 = quat2axang(q2); 
%         r1 = rand; ss = rand < out; ss_out = rand*pi/2; y1 = 2*(rand - 0.5); y2 = 2*(rand - 0.5); 
%         w = q1(1:3); w = (w.*(1-ss) + 0.4*randn(1, 3)); w = w/norm(w); w1 =  w*(r1*0.05+r2*0.1*rand+ss*ss_out);
%         w = q2(1:3); w = (w.*(1-ss) + 0.4*randn(1, 3)); w = w/norm(w); w2 =  w*(r1*0.05+r2*0.1*rand+ss*ss_out); 
%         RR(:,:,i)=w2R(w2)*RR(:,:,i)*w2R(w1)';         
        
        w=randn(3,1);w=w/norm(w)*Sigma*randn(1)/10; y = 2*(rand(1)-0.5); 
        w2 =  [y, sqrt(1 - y^2), 0]*Sigma*randn(1); 
        RR(:,:,i)=RR(:,:,i)*w2R(w2)*w2R(w); 
    end
end
Outlier = round(Outlier*(size(RR,3))); 
if(nargin>3)
    [~,i]=sort(rand(size(RR,3),1));i=i(1:Outlier);
    for j=1:size(i,1)
        RR(:,:,i(j))=RR(:,:,i(j))*w2R((randn*90*pi/180/sqrt(3))*randn(1,3));
    end
end

aa = ones(size(RR, 3), 1); 
aa(i) = 0; 

%% Turns out, its supper impotant for message to flow in the graph 
I2 = [I(2, :); I(1, :)]; 
R2 = RR; 
for i = 1:size(RR, 3)
    R2(:, :, i) = RR(:, :, i)'; 
end

I = [I, I2]; 
RR = cat(3, RR, R2); 
a = [aa; aa]; 

end


