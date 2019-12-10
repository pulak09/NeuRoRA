
%% ** NOT A RELEASE, PLEASE DO NOT SHARE/DISTRIBUTE **  
% Included with supplementary * CVPR * Paper - ID 2243

clear; 
no_sample = 1250; 
err1 = []; 

for ind = 1:no_sample
    d = round(rand * 600); % 200 - 500 images 
    f = rand*0.25;       % pairwise measurements 10 - 50 % 
    k = rand*15;         % Error in the mesurements 5 - 10 covariences  
    s = rand*10/100;    % Outlier percentage 10 - 50 % 
    InitQ = []; Q = []; I = []; QQ = []; Qest = []; 
    [R,RR,I,out]=RandomSO3Graph(400+d, .25+f, 15+k, 0.1+s); 
    
    counts = hist(I(:),1:max(I(:))); 
    [~, idd] = sort(counts, 'descend'); 
    inv_idd = 1:numel(idd); 
    inv_idd(idd) = inv_idd; 
    R = R(:, :, idd); R1 = R(:, :, 1)'; 
    for i = 1:size(R, 3)
        R(:, :, i) = R(:, :, i)*R1; 
    end
    I = inv_idd(I); 
    
    counts = hist(I(:),1:max(I(:))); 
%     plot((counts)/max(counts)); pause(0.5);     
    fprintf('Validating Input data w.r.t Ground truth\n');
%     ValidateSO3Graph(R,RR,I);close all;
    tic; 
    Rest=AverageSO3Graph(RR,I, 'switch', 0);
    t = toc; 
    [Ebest,eall,R2] = CompareRotationGraph(Rest, R); 
    err1 = [err1; [Ebest, t]];
    
    for i=1:size(RR,3);QQ(i, :)=R2q(RR(:,:,i)); end
    for i=1:size(R,3); Q(i, :)=R2q(R(:,:,i)); end
    for i=1:size(R,3); Qest(i, :)=R2q(Rest(:,:,i)); end
    
    out = compute_outliers(QQ, I, Q); 
    disp(['Percentage of outliers = ' num2str(round(sum(1-out)*100/numel(out))) '%']); 
    disp(['Percentage of Edges = ' num2str(round(numel(out)*100/size(R, 3)^2)) '%']);
%     w = compute_relative(Q, QQ, I);
%     ww = sum(sqrt(abs(w.^2 - repmat([1, 0, 0, 0], [size(w, 1), 1]))), 2); 
    InitQ = initialize_quaternions(QQ, I, size(Q, 1)); 
   for ii=1:size(InitQ,1); RO(:,:,ii) = q2R(InitQ(ii, :)); end
   [Ebest,eall,R2] = CompareRotationGraph(RO, R); 

    N = isnan(InitQ) + isnan(Q); % + isnan(QQ); 
    if sum(N(:)) + sum(isnan(QQ(:))) + sum(isnan(I(:))) > 0
        continue;
    end
    
%     no_images = max(I(:)); 
%     new_index = [ [ones(1, no_images-1); 2:no_images], [2:no_images; ones(1, no_images-1)]]; 
%     new_attr = [InitQ(2:end, :); [InitQ(2:end, 1), - InitQ(2:end, 2:4)]]; 
%     mask = 1 - (I(1, :) == 1) | (I(2, :) == 1); 
%     QQ = [QQ(mask, :); new_attr]; 
%     I = [I(:, mask), new_index]; 
    
    data(ind).x = InitQ'; 
    data(ind).xt = Qest'; 
    data(ind).o = out'; 
    data(ind).y = Q';
    data(ind).edge_index = I-1;
    data(ind).edge_feature = QQ'; 
    
end

%%
disp(mean(err1)); 

% filename = '../data/gt_graph_random_large_outliers.h5'; 
% 1.2498    1.1505    1.3920    0.0912 % for the fine algorithm 
% 5.0197    4.2014    6.5590    0.6174 % For large graphs 
%    3.3153    2.5571    5.1500    5.8286 % Large 2
% 6.1606    4.3770   11.4373    0.0463 
% 1.6062    1.1480    3.1679    3.9136

delete(filename); 
h5create(filename, ['/data/' num2str(1) '/x'],  size(data(1).x)); 
h5write(filename, ['/data/' num2str(1) '/x'], data(1).x);
h5create(filename, ['/data/' num2str(1) '/xt'],  size(data(1).xt)); 
h5write(filename, ['/data/' num2str(1) '/xt'], data(1).xt);
h5create(filename, ['/data/' num2str(1) '/y'],  size(data(1).y)); 
h5write(filename, ['/data/' num2str(1) '/y'], data(1).y);
h5create(filename, ['/data/' num2str(1) '/o'],  size(data(1).o)); 
h5write(filename, ['/data/' num2str(1) '/o'], double(data(1).o));
h5create(filename, ['/data/' num2str(1) '/edge_index'],  size(data(1).edge_index)); 
h5write(filename, ['/data/' num2str(1) '/edge_index'], data(1).edge_index);
h5create(filename, ['/data/' num2str(1) '/edge_feature'],  size(data(1).edge_feature)); 
h5write(filename, ['/data/' num2str(1) '/edge_feature'], data(1).edge_feature);
for ind = 2:no_sample
    h5create(filename, ['/data/' num2str(ind) '/x'],  size(data(ind).x)); 
    h5write(filename, ['/data/' num2str(ind) '/x'], data(ind).x);
    h5create(filename, ['/data/' num2str(ind) '/xt'],  size(data(ind).xt)); 
    h5write(filename, ['/data/' num2str(ind) '/xt'], data(ind).xt);
    h5create(filename, ['/data/' num2str(ind) '/y'],  size(data(ind).y)); 
    h5write(filename, ['/data/' num2str(ind) '/y'], data(ind).y);
    h5create(filename, ['/data/' num2str(ind) '/o'],  size(data(ind).o)); 
    h5write(filename, ['/data/' num2str(ind) '/o'], double(data(ind).o));
    h5create(filename, ['/data/' num2str(ind) '/edge_index'],  size(data(ind).edge_index)); 
    h5write(filename, ['/data/' num2str(ind) '/edge_index'], data(ind).edge_index);
    h5create(filename, ['/data/' num2str(ind) '/edge_feature'],  size(data(ind).edge_feature)); 
    h5write(filename, ['/data/' num2str(ind) '/edge_feature'], data(ind).edge_feature);
end
% datat = hdf5read(filename, 'data/0'); 
h5disp(filename) 
disp(mean(err1)); 




%%
function o = compute_outliers(QQ, I, Q)
    i=I(1,:);j=I(2,:);

    w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
        repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
        [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

    w(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
        repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
        [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms

    s2=sqrt(sum(w(:,2:4).*w(:,2:4),2));
    w(:,1)=2*atan2(s2,w(:,1));
    ii=w(:,1)<-pi;  w(ii,1)=w(ii,1)+2*pi;  ii=w(:,1)>=pi;  w(ii,1)=w(ii,1)-2*pi;
    
%     vv = abs(w(:,1))*180/pi; 
    
    o = abs(w(:,1)) < pi/9; 

end

function [w, att] = compute_relative(Q, QQ, I)
    i=I(1,:);j=I(2,:);

    w(:,:)=[ (QQ(:,1).*Q(i,1)-sum(QQ(:,2:4).*Q(i,2:4),2)),...  %scalar terms
        repmat(QQ(:,1),[1,3]).*Q(i,2:4) + repmat(Q(i,1),[1,3]).*QQ(:,2:4) + ...   %vector terms
        [QQ(:,3).*Q(i,4)-QQ(:,4).*Q(i,3),QQ(:,4).*Q(i,2)-QQ(:,2).*Q(i,4),QQ(:,2).*Q(i,3)-QQ(:,3).*Q(i,2)] ];   %cross product terms

    w(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
        repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
        [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms
    
    
    att(:,:)=[ (-Q(j,1).*w(:,1)-sum(Q(j,2:4).*w(:,2:4),2)),...  %scalar terms
        repmat(-Q(j,1),[1,3]).*w(:,2:4) + repmat(w(:,1),[1,3]).*Q(j,2:4) + ...   %vector terms
        [Q(j,3).*w(:,4)-Q(j,4).*w(:,3),Q(j,4).*w(:,2)-Q(j,2).*w(:,4),Q(j,2).*w(:,3)-Q(j,3).*w(:,2)] ];   %cross product terms
end

%%
function Q = initialize_quaternions(QQ, I, N)
    Q=repmat([1,0,0,0],[N,1]);
   
%     No_edges = size(QQ, 1)/2; 
%     I = I(:, 1:No_edges); 
%     QQ = QQ(1:No_edges, :); 
%     
    counts = hist(I(:),1:N); 
%     plot(sort(counts)/max(counts)); 
    [vv, ind] = sort(counts(I(1, :)).*counts(I(2, :)), 'descend'); 
    QQ = QQ(ind, :); 
    I = I(:, ind); 
    
    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);
    %[~,a]=max(hist(sort(I(:)),[1:5530]));
    a=I(1, 1);
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

