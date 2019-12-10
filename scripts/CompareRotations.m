function [E,e]=CompareRotations(R1,R2)

e=zeros(size(R1,3),1);
for i=1:size(R1,3)
    if(sum(sum((R1(:,:,i)==0)))==9||any(any(isnan(R1(:,:,i))))||...
        sum(sum((R2(:,:,i)==0)))==9||any(any(isnan(R2(:,:,i)))));
        e(i,1)=NaN;
    else
        e(i,1)=acos(max(min((R1(1,:,i)*R2(1,:,i)'+R1(2,:,i)*R2(2,:,i)'+R1(3,:,i)*R2(3,:,i)'-1)/2,1),-1));
    end
end
e=e*180/pi;
i=~isnan(e);
E=[mean(e(i)) median(e(i)) sqrt(e(i)'*e(i)/sum(i))];
fprintf('Angular Error (Degrees): Mean=%.2f; Median=%.2f; RMS=%.2f\n',round(E(1,1)*100)/100,round(E(1,2)*100)/100,round(E(1,3)*100)/100);
hist(e,180);
end
