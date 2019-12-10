
%% ** PLEASE DO NOT DISTRIBUTE **  

%%% Multiple rotation averaging using Weiszfeld algorithm http://users.cecs.anu.edu.au/~trumpf/pubs/Hartley_Aftab_Trumpf_CVPR2011.pdf

function R=AverageSO3Graph_weisz(RR,I,Rinit)

N=max(I(:));
[s,c]=graphconncomp(sparse([I(1,:),I(2,:)],[I(2,:),I(1,:)],ones(1,2*size(I,2)),N,N));
m=0;j=0;for i=1:s;k=sum(c==i);if(k>m);m=k;j=i;end;end;
i=find(c==j);  [~,I]=ismember(I,i);  j=all(I);
fprintf('#Cameras = %d; #Edges = %d\n',length(i),sum(j));
Time=zeros(2,1);

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

% The main algorithm 

iter=1;
S(:,:, iter) = Q; dist = 100; distn = 1; 
while isreal(S(:,:,iter)) & iter < 50  
     dist = distn; 
    i=I(1,:);j=I(2,:);
    Q = S(:,:,iter); 
    New_Q = Q; 
    for ii = 1:size(Q, 1)
        id = i == ii; id2 = j == ii; 
        i_nbd = j(id); %j_nbd = i(id2); 
        R_i = Q(ii, :); 
        R_nbdi = compute_nbdi(QQ, Q, id, i_nbd, R_i); %R_nbdj = compute_nbdj(QQ, Q, id2, j_nbd, R_i);
        R_nbd = R_nbdi; %[ R_nbdj]; 
        New_Q(ii, :) = weiszfeld(R_nbd, R_i); 
        
    end
    iter = iter + 1; 
    S(:,:,iter) = New_Q; 
    distn = S(1, :, iter) - S(1, :, iter-1); 
    distn = sum(norm(distn(:))); 
end

Q=S(:,:,iter-1);

if(~QuaternionIP)
    R=zeros(3,3,N);
    for i=1:size(Q,1)
        R(:,:,i)=real(q2R(Q(i,:)));
    end
else
    R=Q';
end

fprintf('Total Computation Time %d seconds\n',round(sum(Time)));

end

function R_nbd = compute_nbdi(QQ, Q, id, i_nbd, R_i)

        R_nbd = [ (-Q(i_nbd,1).*QQ(id,1)-sum(Q(i_nbd,2:4).*QQ(id,2:4),2)),...  %scalar terms
        repmat(-Q(i_nbd,1),[1,3]).*QQ(id,2:4) + repmat(QQ(id,1),[1,3]).*Q(i_nbd,2:4) + ...   %vector terms
        [Q(i_nbd,3).*QQ(id,4)-Q(i_nbd,4).*QQ(id,3),Q(i_nbd,4).*QQ(id,2)-Q(i_nbd,2).*QQ(id,4),Q(i_nbd,2).*QQ(id,3)-Q(i_nbd,3).*QQ(id,2)] ];   %cross product terms
        R_nbd(:, 1) = - R_nbd(:, 1); 
        if R_i(1, 1) > 0
                R_nbd((R_nbd(:, 1) < 0), :) = -R_nbd((R_nbd(:, 1) < 0), :); 
        else
            R_nbd((R_nbd(:, 1) > 0), :) = -R_nbd((R_nbd(:, 1) > 0), :); 
        end
end

% function R_nbd = compute_nbdj(QQ, Q, id, i_nbd, R_i)
% 
%         R_nbd = [ (QQ(id,1).*Q(i_nbd,1)-sum(QQ(id,2:4).*Q(i_nbd,2:4),2)),...  %scalar terms
%                     repmat(QQ(id,1),[1,3]).*Q(i_nbd,2:4) + repmat(Q(i_nbd,1),[1,3]).*QQ(id,2:4) + ...   %vector terms
%                     [QQ(id,3).*Q(i_nbd,4)-QQ(id,4).*Q(i_nbd,3),QQ(id,4).*Q(i_nbd,2)-QQ(id,2).*Q(i_nbd,4),QQ(id,2).*Q(i_nbd,3)-QQ(id,3).*Q(i_nbd,2)] ];   %cross product terms
% 
%         if R_i(1, 1) > 0
%                 R_nbd((R_nbd(:, 1) < 0), :) = -R_nbd((R_nbd(:, 1) < 0), :); 
%         else
%             R_nbd((R_nbd(:, 1) > 0), :) = -R_nbd((R_nbd(:, 1) > 0), :); 
%         end
% end

function Q = weiszfeld(R_i, R)

Q(1, :, 1) = R;
iter = 1; dist = 1; 
while(dist > 0.0001) & iter < 50 
    eps1 = 0.01; 
    lambda = 1./(eps1+sum((R_i - repmat(Q(1, :, iter), [size(R_i, 1), 1, 1])).^2, 2)); 
    iter = iter + 1;
    Q(1, :, iter) = sum(repmat(lambda, [1, 4]).*R_i, 1)/(eps1+ sum(lambda)); 
    Q(1, :, iter) = Q(1, :, iter)/(eps+norm(Q(1, :, iter))); 
    dist = norm(Q(1, :, iter) - Q(1, :, iter-1)); 
end
Q = Q(1, :, iter); 
end

