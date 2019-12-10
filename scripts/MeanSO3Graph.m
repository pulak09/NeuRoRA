%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This implementation is based on V. M. Govindu's Lie Algebraic Motion
% Averaging on SO3
% This code minimizes the geodesic angular distance and NOT the chordal
% distance

% IMPORTANT NOTES:
% We assume the underlying model or equation to be: 
% X'=R*X; Rij=Rj*inv(Ri) i.e. we use camera centered coordinate system
% and NOT the geocentered coordinate for which the underlying equations are
% X'=inv(R)*X; Rij=inv(Ri)*Rj. 
% To use geocentered coordinate please transpose the rotations or change
% the sign of the scalar term of the quaternions before feeding into the
% code and also after getting the result.
%
% Feeding of not connected graph is not recomended.
%
% This code is able to handle inputs in both Rotation matrix as well as
% quaternion format. The Format of output is same as that of the input.

% Programmer: AVISHEK CHATTERJEE
%             PhD Student (S. R. No. 04-03-05-10-12-11-1-08692)
%             Learning System and Multimedia Lab
%             Dept. of Electrical Engineering
%             INDIAN INSTITUTE OF SCIENCE

% Dated:      October 2011

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [R,Iteration,Time]=MeanSO3Graph(RR,I,Rinit,maxIters)
% function [R] = MeanSO3Graph(RR,I,Rinit)
% INPUT  : RR = 'm' number of 3 X 3 Relative Rotation Matrices (R_ij) 
%                       stacked as a 3 X 3 X m Matrix
%                       OR
%                       'm' number of 1 X 4 Relative Quaternions (R_ij) 
%                       stacked as a m X 4  Matrix
%               I  =  Index matrix (ij) of size (2 X m) such that RR(:,:,p)
%                       (OR RR(p,:) for quaternion representation)  is
%                       the relative rotation from R(:,:,I(1,p)) to R(:,:,I(2,p))
%                       (OR R(I(1,p),:) and  R(I(2,p),:) for quaternion representation)
%
% OUTPUT : R  = 'n' number of 3 X 3 Absolute Rotation matrices stacked as
%                        a  3 X 3 X n Matrix 
%                        OR
%                       'n' number of 1 X 4 Relative Quaternions (R_ij) 
%                       stacked as a n X 4  Matrix
tic;

if(nargin<4);maxIters=250;end
changeThreshold=1e-3;

N=max(max(I));%Number of cameras or images or nodes in view graph

QuaternionIP=(size(RR,2)==4);
if(~QuaternionIP)
    %Convert Rij to Quaternion form without function call
    QQ=[RR(1,1,:)+RR(2,2,:)+RR(3,3,:)-1, RR(3,2,:)-RR(2,3,:),RR(1,3,:)-RR(3,1,:),RR(2,1,:)-RR(1,2,:)]/2;
    QQ=reshape(QQ,4,size(QQ,3),1)';
    QQ(:,1)=sqrt((QQ(:,1)+1)/2);
    QQ(:,2:4)=(QQ(:,2:4)./repmat(QQ(:,1),[1,3]))/2;
else
    QQ=RR;
end

if(nargin>2&&(~isempty(Rinit)))
    Q=[Rinit(1,1,:)+Rinit(2,2,:)+Rinit(3,3,:)-1, Rinit(3,2,:)-Rinit(2,3,:),Rinit(1,3,:)-Rinit(3,1,:),Rinit(2,1,:)-Rinit(1,2,:)]/2;
	Q=reshape(Q,4,size(Q,3),1)';
	Q(:,1)=sqrt((Q(:,1)+1)/2);
	Q(:,2:4)=(Q(:,2:4)./repmat(Q(:,1),[1,3]))/2;
else
    Q=repmat([1,0,0,0],[N,1]);
    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);    i(1)=1;
    while(sum(i)<N)
       SpanFlag=0;
        for j=1:size(I,2)
            if(i(I(1,j))==1&&i(I(2,j))==0)
                %Rinit(:,:,I(2,j))=RR(:,:,j)*Rinit(:,:,I(1,j)); %Dont Uncomment
                Q(I(2,j),:)=[ (QQ(j,1).*Q(I(1,j),1)-sum(QQ(j,2:4).*Q(I(1,j),2:4),2)),...  %scalar terms
                    repmat(QQ(j,1),[1,3]).*Q(I(1,j),2:4) + repmat(Q(I(1,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(1,j),4)-QQ(j,4).*Q(I(1,j),3),QQ(j,4).*Q(I(1,j),2)-QQ(j,2).*Q(I(1,j),4),QQ(j,2).*Q(I(1,j),3)-QQ(j,3).*Q(I(1,j),2)] ];   %cross product terms
                i(I(2,j))=1;
                SpanFlag=1;
            end
            if(i(I(1,j))==0&&i(I(2,j))==1)
                %Rinit(:,:,I(1,j))=RR(:,:,j)'*Rinit(:,:,I(2,j)); %Dont Uncomment
                Q(I(1,j),:)=[ (-QQ(j,1).*Q(I(2,j),1)-sum(QQ(j,2:4).*Q(I(2,j),2:4),2)),...  %scalar terms
                    repmat(-QQ(j,1),[1,3]).*Q(I(2,j),2:4) + repmat(Q(I(2,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(2,j),4)-QQ(j,4).*Q(I(2,j),3),QQ(j,4).*Q(I(2,j),2)-QQ(j,2).*Q(I(2,j),4),QQ(j,2).*Q(I(2,j),3)-QQ(j,3).*Q(I(2,j),2)] ];   %cross product terms
                i(I(1,j))=1;
                SpanFlag=1;
            end
        end
        if(SpanFlag==0&&sum(i)<N)
            warning('Relative rotations DO NOT SPAN all the nodes in the VIEW GRAPH');
            pause();break;
        end
    end    
end

% Formation of A matrix.
m=size(I,2);
i=[(1:m);(1:m)];i=i(:);
j=I(:);
s=repmat([-1;1],[m,1]);
k=(j~=1);
Amatrix=sparse(i(k),j(k)-1,s(k),m,N-1);

w=zeros(size(QQ,1),4);W=zeros(N,4);

score=inf;    Iteration=0;

fprintf('Iteration: %4d; Time: %7.1f; MaxChange: %8.4f',0,toc,0);

while((score>changeThreshold)&&(Iteration<maxIters))

    i=I(1,:);j=I(2,:);

    % w=Qij*Qi
    w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
        repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
        [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

    % w=inv(Qj)*w=inv(Qj)*Qij*Qi
    w(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
        repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
        [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms


    s2=sqrt(sum(w(:,2:4).*w(:,2:4),2));
    w(:,1)=2*atan2(s2,w(:,1));
    i=w(:,1)<-pi;  w(i,1)=w(i,1)+2*pi;  i=w(:,1)>=pi;  w(i,1)=w(i,1)-2*pi;
    B=w(:,2:4).*repmat(w(:,1)./s2,[1,3]);
% Here is an alternative solution for the above 4 lines. This may be
% marginally faster. But use of this is not recomended as the domain of
% acos is bounded which may result in truncation error when the solution
% comes near optima. Usage of atan2 justifies omition of explicit
% quaternion normalization at every stage.
%     i=w(:,1)<0;w(i,:)=-w(i,:);
%     theta2=acos(w(:,1));
%     B=((w(:,2:4).*repmat((2*theta2./sin(theta2)),[1,3])));
    
    
    B(isnan(B))=0;% This tackles the devide by zero problem.
  
    W(1,:)=[1 0 0 0];
    
    W(2:end,2:4)=Amatrix\B;
    
    score=max(sqrt(sum(W(2:end,2:4).*W(2:end,2:4),2)));
    
    theta=sqrt(sum(W(:,2:4).*W(:,2:4),2));
    W(:,1)=cos(theta/2);
    W(:,2:4)=W(:,2:4).*repmat(sin(theta/2)./theta,[1,3]);
    
    W(isnan(W))=0;
    
    Q=[ (Q(:,1).*W(:,1)-sum(Q(:,2:4).*W(:,2:4),2)),...  %scalar terms
        repmat(Q(:,1),[1,3]).*W(:,2:4) + repmat(W(:,1),[1,3]).*Q(:,2:4) + ...   %vector terms
        [Q(:,3).*W(:,4)-Q(:,4).*W(:,3),Q(:,4).*W(:,2)-Q(:,2).*W(:,4),Q(:,2).*W(:,3)-Q(:,3).*W(:,2)] ];   %cross product terms

    Iteration=Iteration+1;
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\bIteration: %4d; Time: %7.1f; MaxChange: %8.4f',Iteration,toc,score);
end;
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\n');

if(~QuaternionIP)
    R=zeros(3,3,N);
    for i=1:size(Q,1)
        R(:,:,i)=q2R(Q(i,:));
    end
else
    R=Q;
end

if(Iteration>=maxIters);disp('Max iterations reached');end;
if(nargout==3);Time=toc;else toc;end;
end
