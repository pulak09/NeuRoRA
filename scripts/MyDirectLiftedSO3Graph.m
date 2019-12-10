
function  [R,Iteration,Time]=MyDirectLiftedSO3Graph(RR,I,Rinit,maxIters)

tic;

if(nargin<4 || isempty(maxIters));maxIters=100;end

N=max(max(I));%Number of cameras or images or nodes in view graph

QuaternionIP=(size(RR,1)==4);
if(~QuaternionIP)
    %Convert Rij to Quaternion form without function call
    QQ=[RR(1,1,:)+RR(2,2,:)+RR(3,3,:)-1, RR(3,2,:)-RR(2,3,:),RR(1,3,:)-RR(3,1,:),RR(2,1,:)-RR(1,2,:)]/2;
    QQ=reshape(QQ,4,size(QQ,3),1)';
    QQ(:,1)=sqrt((QQ(:,1)+1)/2);
    QQ(:,2:4)=(QQ(:,2:4)./repmat(QQ(:,1),[1,3]))/2;
else
    QQ=RR';
end

if(nargin>2 && (~isempty(Rinit)))
    if(size(Rinit,1)==3)
        Q=[Rinit(1,1,:)+Rinit(2,2,:)+Rinit(3,3,:)-1, Rinit(3,2,:)-Rinit(2,3,:),Rinit(1,3,:)-Rinit(3,1,:),Rinit(2,1,:)-Rinit(1,2,:)]/2;
        Q=reshape(Q,4,size(Q,3),1)';
        Q(:,1)=sqrt((Q(:,1)+1)/2);
        Q(:,2:4)=(Q(:,2:4)./repmat(Q(:,1),[1,3]))/2;
    else
        Q=Rinit';
    end
else
    Q=repmat([1,0,0,0],[N,1]);
    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);
    %[~,a]=max(hist(sort(I(:)),[1:5530]));
    a=1;
    i(a)=1;
    while(sum(i)<N)
       SpanFlag=0;
        for j=1:size(I,2)
            if(i(I(1,j))==1&&i(I(2,j))==0)
                %Rinit(:,:,I(2,j))=RR(:,:,j)*Rinit(:,:,I(1,j));
                Q(I(2,j),:)=[ (QQ(j,1).*Q(I(1,j),1)-sum(QQ(j,2:4).*Q(I(1,j),2:4),2)),...  %scalar terms
                    repmat(QQ(j,1),[1,3]).*Q(I(1,j),2:4) + repmat(Q(I(1,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(1,j),4)-QQ(j,4).*Q(I(1,j),3),QQ(j,4).*Q(I(1,j),2)-QQ(j,2).*Q(I(1,j),4),QQ(j,2).*Q(I(1,j),3)-QQ(j,3).*Q(I(1,j),2)] ];   %cross product terms
                i(I(2,j))=1;
                SpanFlag=1;
            end
            if(i(I(1,j))==0&&i(I(2,j))==1)
                %Rinit(:,:,I(1,j))=RR(:,:,j)'*Rinit(:,:,I(2,j));
                Q(I(1,j),:)=[ (-QQ(j,1).*Q(I(2,j),1)-sum(QQ(j,2:4).*Q(I(2,j),2:4),2)),...  %scalar terms
                    repmat(-QQ(j,1),[1,3]).*Q(I(2,j),2:4) + repmat(Q(I(2,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(2,j),4)-QQ(j,4).*Q(I(2,j),3),QQ(j,4).*Q(I(2,j),2)-QQ(j,2).*Q(I(2,j),4),QQ(j,2).*Q(I(2,j),3)-QQ(j,3).*Q(I(2,j),2)] ];   %cross product terms
                i(I(1,j))=1;
                SpanFlag=1;
            end
        end
        if(SpanFlag==0&&sum(i)<N)
            fprintf('Relative rotations DO NOT SPAN all the nodes in the VIEW GRAPH');
            fprintf('Number of nodes in Spanning Tree = %d\n',sum(i));
            fprintf('Connected Nodes are given as output\n');
            fprintf('Remove extra nodes and retry\n');
            R=i;
            return;
        end
    end    
end

m=size(I,2);
i=[(1:m);(1:m)];i=i(:);
j=I(:);
s=repmat([-1;1],[m,1]);
k=(j~=1);

n=N-1;    i=I(1,:)-1;    j=I(2,:)-1;
r=[n*(i-1)+i,n*(j-1)+j,n*(i-1)+j,n*(j-1)+i]';
c=[1:m,1:m,1:m,1:m]';
s=[ones(2*m,1);-ones(2*m,1)];
k=[i&i,j&j,i&j,j&i];
r=r(k,1);  c=c(k,1);  s=s(k,1);
AtA=sparse(r,c,s,n*n,m);

%% 
i=I(1,:);j=I(2,:);

w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
    repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
    [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

w(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
    repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
    [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms

w = w(:,1:4).*w(:,1:4); 
s2=(sum(w,2));
w=w./repmat(s2,[1,4]);
    
sig1 = log(sum((0.1 ./(0.01+(w(:, 2:4)))), 2)) - 1;
% sig1 = sig1/max(sig1); 
% w=zeros(size(QQ,1),4);W=zeros(N,4);

score=inf;    Iteration=0;    L1Step=2;

fprintf('Iteration: %4d; Time: %7.1f; MaxChange: %8.4f',0,toc,0);

optmincon = optimoptions(@fmincon, 'MaxIter', 200, 'TolFun', 1.0000e-5); 
optmincon = optimoptions(optmincon,'GradObj','off','GradConstr','on', 'Algorithm','interior-point');
optmincon = optimoptions(optmincon, 'TolX', 1.0000e-5, 'Display',  'iter', 'Diagnostics', 'on'); 
optmincon = optimoptions(optmincon, 'MaxFunEvals', 1000000);%, 'Hessian', 'user-supplied', 'HessFcn', @hessinterior);%, 'TolX', 0, 'TolFun', 0);

options.ub = [ones(numel(Q)-4, 1); Inf*ones(m, 1)];   % Lower bound on the variables.
options.lb = [-ones(numel(Q)-4, 1); -Inf*ones(m, 1)];  % Upper bound on the variables.
options.MaxFunctionEvaluations = 100000;

%     hPixelO = invK*[linesPts, ones(size(linesPts, 1), 1)]'; 
%     [xopt] = fmincon(@objective_fmincon, x0, A, Ab, [], [], options.lb, options.ub, [], optmincon);%@constraints_fmincon, optmincon);

Q2 = Q(2:end, :); 
x0 = [Q2(:); sig1];  
% x0 = [Q2(:); 10000.0*ones(m, 1)];  
xoptm = fmincon(@objective_fmincon, x0, [], [], [], [], options.lb, options.ub, [], optmincon);
 
Q = [[1, 0, 0, 0]; reshape(xoptm(1:end-m), [N-1, 4])]; 

Q1 = Q(:,1:4).*Q(:,1:4); 
s1=sqrt(sum(Q1,2));
Q=Q./repmat(s1,[1,4]);
    
if(~QuaternionIP)
    R=zeros(3,3,N);
    for i=1:size(Q,1)
        R(:,:,i)=real(q2R(Q(i,:)));
    end
else
    R=Q';
end

if(Iteration>=maxIters);fprintf(' (Max Iteration)');end;fprintf('\n');
if(nargout==3);Time=toc;else toc;end;


function [f, g] = objective_fmincon(x) 
    Qq = [[1, 0, 0, 0]; reshape(x(1:end-m), [N-1, 4])]; 
    sig = 1.0 ./ (1.0 + exp(-x(end-m+1:end))); threshold = 0.005; 
    
    w(:,:)=[ (QQ(:,1).*Qq(i,1)-sum(QQ(:,2:4).*Qq(i,2:4),2)),...  %scalar terms
        repmat(QQ(:,1),[1,3]).*Qq(i,2:4) + repmat(Qq(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
        [QQ(:,3).*Qq(i,4)-QQ(:,4).*Qq(i,3),QQ(:,4).*Qq(i,2)-QQ(:,2).*Qq(i,4),QQ(:,2).*Qq(i,3)-QQ(:,3).*Qq(i,2)] ];   %cross product terms

    w(:,:)=[ (-Qq(j,1).*w(:,1)-sum(Qq(j,2:4).*w(:,2:4),2)),...  %scalar terms
        repmat(-Qq(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Qq(j,2:4) + ...   %vector terms
        [Qq(j,3).*w(:,4)-Qq(j,4).*w(:,3),Qq(j,4).*w(:,2)-Qq(j,2).*w(:,4),Qq(j,2).*w(:,3)-Qq(j,3).*w(:,2)] ];   %cross product terms

    w = w(:,1:4).*w(:,1:4); 
    s2=(sum(w,2));
    WW=w(:, 2:4)./repmat(s2,[1,3]);
    
%     C = repmat([1, 0, 0, 0], [m, 1]); 
    diff = sum(WW, 2); 
    f =  sig.^2 .* diff + threshold*(1 - sig.^2); 
    f = sum(f(:)); 
    g = []; 
end

end
