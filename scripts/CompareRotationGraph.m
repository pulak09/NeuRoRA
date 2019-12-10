function [Ebest,eall,R2]=CompareRotationGraph(R1,R2,AlignMethod)
i=min([size(R1,3),size(R2,3)]);R1(:,:,i+1:end)=[];R2(:,:,i+1:end)=[];
i=~(any(any(isnan(R1),1),2)|sum(sum(R1.*R1,1),2)==0|any(any(isnan(R2),1),2)|sum(sum(R2.*R2,1),2)==0);
R1=R1(:,:,i);R2=R2(:,:,i);I=i;

if(nargin<3);AlignMethod='median';end;
SIGMA2=(5*pi/180)^2;
N=size(R1,3);

% e=zeros(N,1);
% for i=1:N-1
%     if(~mod(i,10000));disp(num2str(i));end;
%     for j=i+1:N
%         a=acos(max(min((R1(1,:,i)*R1(1,:,j)'+R1(2,:,i)*R1(2,:,j)'+R1(3,:,i)*R1(3,:,j)'-1)/2,1),-1));
%         e(i,1)=e(i,1)+a;
%         e(j,1)=e(j,1)+a;
%     end
% end
% [~, j]=min(e);

Emeanbest=inf;E=[0 0 0];Ebest=E;e=zeros(N,1);ebest=e;
for t=1:4
    j=randi(N,1);
    R=R1(:,:,j)';  for i=1:N;R1(:,:,i)=R1(:,:,i)*R;end;
    R=R2(:,:,j)';  for i=1:N;R2(:,:,i)=R2(:,:,i)*R;end;

    W=zeros(N,3);d=inf; count = 1; 
    while(d>1e-5 && count < 20)
        for i=1:N
            W(i,:)=R2w(R2(:,:,i)'*R1(:,:,i));
        end
        if(strcmp(AlignMethod,'mean'))
            w=mean(W);
            d=norm(w);    R=w2R(w);
        elseif(strcmp(AlignMethod,'median'))
            w=median(W);
            d=norm(w);    R=w2R(w);
        elseif(strcmp(AlignMethod,'robustmean'))
            w=1./sqrt((sum((W.*W),2)+SIGMA2));w=w/sum(w);
            w=mean(repmat(w,[1,3]).*W);
            d=norm(w);    R=w2R(w);
        else
            error('Undefined AlignMethod');
        end
        for i=1:N
            R2(:,:,i)=R2(:,:,i)*R;
        end
        count = count + 1; 
    end

    for i=1:N
        %acos(trace(R2*R1'))
        e(i,1)=acos(max(min((R1(1,:,i)*R2(1,:,i)'+R1(2,:,i)*R2(2,:,i)'+R1(3,:,i)*R2(3,:,i)'-1)/2,1),-1));
    end
    e=e*180/pi;
    E=[mean(e) median(e) sqrt(e'*e/size(e,1))];
    if(E(2)<Emeanbest);
        ebest=e;Ebest=E; Emeanbest = E(2); 
    end
end
fprintf('#Common=%d; Angular Error (Degrees): Mean=%.2f; Median=%.2f; RMS=%.2f\n',length(ebest),round(Ebest(1,1)*100)/100,round(Ebest(1,2)*100)/100,round(Ebest(1,3)*100)/100);
% hist(ebest,180);
eall=nan(length(I),1);eall(I)=ebest;
end
