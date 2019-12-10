
%% ** NOT A RELEASE, PLEASE DO NOT SHARE/DISTRIBUTE **  
% Included with supplementary * CVPR * Paper - ID 2243

clear; 

%file = 'gt_graph_random_fine'; 
 file = 'gt_graph_random_large_outliers';
file = 'gt_graph_random_large_outliers_test';
% file = 'gt_graph_random_large_outliers_real';
% file = 'gt_graph_random_large_outliers2'; 
% file = 'gt_graph_random_large_landmark10k2'; 
% file = 'gt_graph_random_large_landmark10k'; 
filename = ['../data/', file,'_pred_rot.h5']; 
%h5disp(filename); 

filename2 = ['../data/', file, '.h5']; 

no_images = 100; 
val_pred = []; 
val_pred1 = []; 
val_pred2 = []; 

err1 = []; 
err2 = []; 

ell1 = []; 
ell2 = []; 
count = 1; 
for ind = 1:no_images 
    R = [];
    RO = []; 
    RR = [];
    out_predicted_org = double(h5read(filename, ['/data/', num2str(ind), '/ot'])); 
    out_predicted = double(out_predicted_org > 0.5); 
    out_original = double(h5read(filename, ['/data/', num2str(ind), '/o'])); 
    prct = 100*sum(abs(out_predicted - out_original)) / numel(out_original); 
%     figure(1), hist(out_predicted)
%     figure(2), hist(out_original)
    ind1 = (out_original == 0); 
    ind2 = (out_original == 1); 
    prct1 = 100*sum(abs(out_predicted(ind1) - out_original(ind1))) / sum(ind1); 
    prct2 = 100*sum(abs(out_predicted(ind2) - out_original(ind2))) / sum(ind2); 
    val_pred = [val_pred, prct]; 
    val_pred1 = [val_pred1, prct1]; 
    val_pred2 = [val_pred2, prct2]; 
    
    QQ_mod = double(h5read(filename, ['/data/', num2str(ind), '/refined_qq']))'; 
    Q = double(h5read(filename, ['/data/', num2str(ind), '/y']))'; 
    I = double(h5read(filename, ['/data/', num2str(ind), '/edge_index']))'+1; 
    QQ = double(h5read(filename, ['/data/', num2str(ind), '/edge_feature']))'; 
    
%     InitQ = double(h5read(filename2, ['/data/', num2str(ind), '/x']))'; 
%     a = 1; 
    InitQ = []; th = 0.25; 
%     [InitQ, a] = initialize_quaternions(QQ_mod, I, size(Q, 1), out_predicted_org, th);
    [InitQ, a] = initialize_quaternions_SPT(QQ_mod, I, size(Q, 1), out_predicted_org, th);
    Qest = double(h5read(filename, ['/data/', num2str(ind), '/xt']))'; 

%     Qest = AverageSO3Graph_weisz(QQ', I)'; 
%     while(numel(InitQ) < 1 && th > 0.5) 
%         count = 0; 
%         InitQ = initialize_quaternions(QQ_mod, I, size(Q, 1), out_predicted_org, th); 
%         th = th - 0.05; 
%         if numel(InitQ) < 1
%             disp(count); count = count+1; 
%         end
%     end
    RO = []; R = []; RT = []; 
    for ii=1:size(InitQ,1); RO(:,:,ii) = q2R(InitQ(ii, :)); end
    for ii=1:size(InitQ,1); R(:,:,ii) = q2R(Q(ii, :)); end
    for ii=1:size(Qest,1); RT(:,:,ii) = q2R(Qest(ii, :)); end
    
    for ii=1:size(QQ,1); RR(:,:,ii) = q2R(QQ_mod(ii, :)); end
    
    R1 = RO(:, :, a)'; R2 = RT(:, :, a)'; R3 = R(:, :, a)'; 
    for i = 1:size(RO, 3)
        RO(:, :, i) = RO(:, :, i)*R1; 
        RT(:, :, i) = RT(:, :, i)*R2; 
        R(:, :, i) = R(:, :, i)*R3; 
    end
    
    for i=1:size(RO,3); InitQ(i, :)=R2q(RO(:,:,i)); end
    for i=1:size(RT,3); Qest(i, :)=R2q(RT(:,:,i)); end
    for i=1:size(R,3); Q(i, :)=R2q(R(:,:,i)); end
    
%     tic; 
%     Rest=AverageSO3Graph(RR,I, 'switch', 0);
%     t = toc; 
    disp([ind]); 
    [Ebest1,eall1,R2] = CompareRotationGraph(RO, R); 
    [Ebest2,eall2,R2] = CompareRotationGraph(RT, R); 
    
    err1 = [err1; Ebest1];
    err2 = [err2; Ebest2];
    ell1 = [ell1; eall1];
    ell2 = [ell2; eall2];

    if Ebest1(1) > 90
        disp(ind)
        continue;
    end
        
    No_edges = numel(out_predicted)/2; 
    out_predicted = out_predicted_org(1:No_edges).*out_predicted_org(No_edges+1:end) > 0.1; 
    out_predicted = [out_predicted, out_predicted]; 
        
    data(count).x = InitQ'; 
    data(count).xt = Qest'; 
    data(count).o = out_original(out_predicted); 
    data(count).y = Q';
    data(count).edge_index = I(:, out_predicted)-1;
    data(count).edge_feature = QQ(out_predicted, :)'; 
    count = count + 1; 
end

% disp([mean(val_pred), mean(val_pred1), mean(val_pred2)]); 
disp(mean(err1)); 
disp(mean(err2)); 

%%
figure, 
n = 180; 
[n, xout] = hist(ell1, [1:5:n]); 
bar(xout, n+1, 'barwidth', 1, 'basevalue', 1, 'FaceColor',[0.9 0 .1],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
set(gca,'YScale','log');
legend('CleanNet'); 

%%

figure, 
n = 180; 
[n, xout] = hist(ell2, [1:5:n]); 
bar(xout, n+1, 'barwidth', 1, 'basevalue', 1, 'FaceColor',[.5 0 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5);
set(gca,'YScale','log');
legend('Weiszfeld'); 


%%
no_images = count - 1; 

filename = ['../data/', file, '_updated_nocleannet.h5']; 
% 1.2498    1.1505    1.3920    0.0912 % for the fine algorithm 
% 1.0780    0.9828    1.1820

delete(filename); 
h5create(filename, ['/data/' num2str(1) '/x'],  size(data(1).x)); 
h5write(filename, ['/data/' num2str(1) '/x'], data(1).x);
h5create(filename, ['/data/' num2str(1) '/xt'],  size(data(1).xt)); 
h5write(filename, ['/data/' num2str(1) '/xt'], data(1).xt);
h5create(filename, ['/data/' num2str(1) '/y'],  size(data(1).y)); 
h5write(filename, ['/data/' num2str(1) '/y'], data(1).y);
h5create(filename, ['/data/' num2str(1) '/o'],  size(data(1).o)); 
h5write(filename, ['/data/' num2str(1) '/o'], data(1).o);

h5create(filename, ['/data/' num2str(1) '/onode'],  [1, size(data(1).x, 2)]); 
h5write(filename, ['/data/' num2str(1) '/onode'], ones([1, size(data(1).x, 2)]));

h5create(filename, ['/data/' num2str(1) '/edge_index'],  size(data(1).edge_index)); 
h5write(filename, ['/data/' num2str(1) '/edge_index'], data(1).edge_index);
h5create(filename, ['/data/' num2str(1) '/edge_feature'],  size(data(1).edge_feature)); 
h5write(filename, ['/data/' num2str(1) '/edge_feature'], data(1).edge_feature);
for ind = 2:no_images
    h5create(filename, ['/data/' num2str(ind) '/x'],  size(data(ind).x)); 
    h5write(filename, ['/data/' num2str(ind) '/x'], data(ind).x);
    h5create(filename, ['/data/' num2str(ind) '/xt'],  size(data(ind).xt)); 
    h5write(filename, ['/data/' num2str(ind) '/xt'], data(ind).xt);
    h5create(filename, ['/data/' num2str(ind) '/y'],  size(data(ind).y)); 
    h5write(filename, ['/data/' num2str(ind) '/y'], data(ind).y);
    h5create(filename, ['/data/' num2str(ind) '/o'],  size(data(ind).o)); 
    h5write(filename, ['/data/' num2str(ind) '/o'], data(ind).o);
    
    h5create(filename, ['/data/' num2str(ind) '/onode'],  [1, size(data(ind).x, 2)]); 
    h5write(filename, ['/data/' num2str(ind) '/onode'], ones([1, size(data(ind).x, 2)]));

    h5create(filename, ['/data/' num2str(ind) '/edge_index'],  size(data(ind).edge_index)); 
    h5write(filename, ['/data/' num2str(ind) '/edge_index'], data(ind).edge_index);
    h5create(filename, ['/data/' num2str(ind) '/edge_feature'],  size(data(ind).edge_feature)); 
    h5write(filename, ['/data/' num2str(ind) '/edge_feature'], data(ind).edge_feature);
end
% datat = hdf5read(filename, 'data/0'); 
% h5disp(filename) 
% disp(mean(err1)); 

    %%
    function Q = initialize_quaternions_MST(QQ, I, N, out_predicted, th)
    Q=repmat([1,0,0,0],[N,1]);

    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);
    %[~,a]=max(hist(sort(I(:)),[1:5530]));
    G=zeros(N,N); 
    for j=1:size(I,2)
        G(I(1, j), I(2, j)) =  1 - out_predicted(j); 
    end
    [MST, ~]=graphminspantree(sparse(G));
    a=N;
    i(a)=1;
    while(sum(i)<N)
       SpanFlag=0;
        for j=1:size(I,2)
            if(i(I(1,j))==1&&i(I(2,j))==0&&MST(I(1,j), I(2,j)) > 0)
                %Rinit(:,:,I(2,j))=RR(:,:,j)*Rinit(:,:,I(1,j));
                Q(I(2,j),:)=[ (QQ(j,1).*Q(I(1,j),1)-sum(QQ(j,2:4).*Q(I(1,j),2:4),2)),...  %scalar terms
                    repmat(QQ(j,1),[1,3]).*Q(I(1,j),2:4) + repmat(Q(I(1,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(1,j),4)-QQ(j,4).*Q(I(1,j),3),QQ(j,4).*Q(I(1,j),2)-QQ(j,2).*Q(I(1,j),4),QQ(j,2).*Q(I(1,j),3)-QQ(j,3).*Q(I(1,j),2)] ];   %cross product terms
                i(I(2,j))=1;
                SpanFlag=1;
            end
            if(i(I(1,j))==0&&i(I(2,j))==1&&MST(I(1,j), I(2,j)) > 0)
                %Rinit(:,:,I(1,j))=RR(:,:,j)'*Rinit(:,:,I(2,j));
                Q(I(1,j),:)=[ (-QQ(j,1).*Q(I(2,j),1)-sum(QQ(j,2:4).*Q(I(2,j),2:4),2)),...  %scalar terms
                    repmat(-QQ(j,1),[1,3]).*Q(I(2,j),2:4) + repmat(Q(I(2,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(2,j),4)-QQ(j,4).*Q(I(2,j),3),QQ(j,4).*Q(I(2,j),2)-QQ(j,2).*Q(I(2,j),4),QQ(j,2).*Q(I(2,j),3)-QQ(j,3).*Q(I(2,j),2)] ];   %cross product terms
                i(I(1,j))=1;
                SpanFlag=1;
            end
        end
%         if(SpanFlag==0&&sum(i)<N)
%             fprintf('Relative rotations DO NOT SPAN all the nodes in the VIEW GRAPH');
%             fprintf('Number of nodes in Spanning Tree = %d\n',sum(i));
%             fprintf('Connected Nodes are given as output\n');
%             fprintf('Remove extra nodes and retry\n');
%             R=i;
% %             Q = []; 
%             return;
%         end
    end    
    end

function [Q, a] = initialize_quaternions(QQ, I, N, out_predicted, th)
    Q=repmat([1,0,0,0],[N,1]);
    No_edges = numel(out_predicted)/2; 
    out_predicted = out_predicted(1:No_edges).*out_predicted(No_edges+1:end); 
    I = I(:, 1:No_edges); 
    out_predicted(isnan(out_predicted)) = 1; 

%     counts = hist(I(:),1:max(I(:))); 
%     [~, idd] = sort(counts, 'descend'); 
%     inv_idd = 1:numel(idd); 
%     inv_idd(idd) = inv_idd; 
%     I = inv_idd(I); 
%     QQ = inv_idd(QQ); 
%     out_predicted = out_predicted(ind); 
 
    counts = hist(I(:),1:N); 
%     plot(sort(counts)/max(counts)); 
    [vv, ind] = sort(counts(I(1, :)).*counts(I(2, :)), 'descend'); 
    QQ = QQ(ind, :); 
    I = I(:, ind); 
    out_predicted = out_predicted(ind); 
%     
    
    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);
    %[~,a]=max(hist(sort(I(:)),[1:5530]));
    a=I(1, 1);
    i(a)=1;
    while(sum(i)<N)
       SpanFlag=0;
        for j=1:size(I,2)
            if(i(I(1,j))==1&&i(I(2,j))==0&&out_predicted(j) > th)
                %Rinit(:,:,I(2,j))=RR(:,:,j)*Rinit(:,:,I(1,j));
                Q(I(2,j),:)=[ (QQ(j,1).*Q(I(1,j),1)-sum(QQ(j,2:4).*Q(I(1,j),2:4),2)),...  %scalar terms
                    repmat(QQ(j,1),[1,3]).*Q(I(1,j),2:4) + repmat(Q(I(1,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(1,j),4)-QQ(j,4).*Q(I(1,j),3),QQ(j,4).*Q(I(1,j),2)-QQ(j,2).*Q(I(1,j),4),QQ(j,2).*Q(I(1,j),3)-QQ(j,3).*Q(I(1,j),2)] ];   %cross product terms
                i(I(2,j))=1;
                SpanFlag=1;
            end
            if(i(I(1,j))==0&&i(I(2,j))==1&&out_predicted(j) > th)
                %Rinit(:,:,I(1,j))=RR(:,:,j)'*Rinit(:,:,I(2,j));
                Q(I(1,j),:)=[ (-QQ(j,1).*Q(I(2,j),1)-sum(QQ(j,2:4).*Q(I(2,j),2:4),2)),...  %scalar terms
                    repmat(-QQ(j,1),[1,3]).*Q(I(2,j),2:4) + repmat(Q(I(2,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(2,j),4)-QQ(j,4).*Q(I(2,j),3),QQ(j,4).*Q(I(2,j),2)-QQ(j,2).*Q(I(2,j),4),QQ(j,2).*Q(I(2,j),3)-QQ(j,3).*Q(I(2,j),2)] ];   %cross product terms
                i(I(1,j))=1;
                SpanFlag=1;
            end
        end
        if(SpanFlag==0&&sum(i)<N)
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
%             Q = Q(idd, :); 
%             Q = []; 
            return;
        end
        end
    end    
%     Q = Q(idd, :); 
end


function [Q, a] = initialize_quaternions_SPT(QQ, I, N, out_predicted, th)
    
    Q=repmat([1,0,0,0],[N,1]);
    No_edges = numel(out_predicted)/2; 
    out_predicted = out_predicted(1:No_edges).*out_predicted(No_edges+1:end); 
    I = I(:, 1:No_edges); 
    QQ = QQ(1:No_edges, :); 
    out_predicted(isnan(out_predicted)) = 0; 

%     counts = hist(I(:),1:max(I(:))); 
%     [~, idd] = sort(counts, 'descend'); 
%     inv_idd = 1:numel(idd); 
%     inv_idd(idd) = inv_idd; 
%     I = inv_idd(I); 
%     QQ = inv_idd(QQ); 
%     out_predicted = out_predicted(ind); 
 
   out_predicted = double((out_predicted > 0.5)); 
    counts = hist(I(:),1:N); 
    [~, id] = max(counts); 
    
    G = graph(I(1, :), I(2, :), 0.1 + 1.0 - out_predicted);
    TR = shortestpathtree(G, 'all', id); 
    
    final_edges = table2array(TR.Edges);
    [vv, ind] = sort(counts(final_edges(:, 2)), 'descend'); 
    
    final_edges = final_edges(ind, :);
    
    %Compute initial Q from a Spanning Tree
    i=zeros(N,1);
    %[~,a]=max(hist(sort(I(:)),[1:5530]));
    a = id;

    for count = 1:size(final_edges, 1)
       SpanFlag=0;
       j = find(I(1, :)*10000+I(2, :) == [final_edges(count, 2)*10000 + final_edges(count, 1)]);  
       if numel(j) > 0
                Q(I(2,j),:)=[ (QQ(j,1).*Q(I(1,j),1)-sum(QQ(j,2:4).*Q(I(1,j),2:4),2)),...  %scalar terms
                    repmat(QQ(j,1),[1,3]).*Q(I(1,j),2:4) + repmat(Q(I(1,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(1,j),4)-QQ(j,4).*Q(I(1,j),3),QQ(j,4).*Q(I(1,j),2)-QQ(j,2).*Q(I(1,j),4),QQ(j,2).*Q(I(1,j),3)-QQ(j,3).*Q(I(1,j),2)] ];   %cross product terms
                i(I(2,j))=1;
                SpanFlag=1;
       end
       j = find(I(1, :)*10000+I(2, :) == [final_edges(count, 1) * 10000 + final_edges(count, 2)]);  
       if numel(j) > 0
                %Rinit(:,:,I(1,j))=RR(:,:,j)'*Rinit(:,:,I(2,j));
                Q(I(1,j),:)=[ (-QQ(j,1).*Q(I(2,j),1)-sum(QQ(j,2:4).*Q(I(2,j),2:4),2)),...  %scalar terms
                    repmat(-QQ(j,1),[1,3]).*Q(I(2,j),2:4) + repmat(Q(I(2,j),1),[1,3]).*QQ(j,2:4) + ...   %vector terms
                    [QQ(j,3).*Q(I(2,j),4)-QQ(j,4).*Q(I(2,j),3),QQ(j,4).*Q(I(2,j),2)-QQ(j,2).*Q(I(2,j),4),QQ(j,2).*Q(I(2,j),3)-QQ(j,3).*Q(I(2,j),2)] ];   %cross product terms
                i(I(1,j))=1;
                SpanFlag=1;
            
       end

    end    
%     Q = Q(idd, :); 
end


