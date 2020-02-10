function [w]=R2w(R)

%R2W	Rotation matrix to vector
%
%	Usage : [w]=R2w(R)
%
%	R is the 3x3 rotation matrix
%       w : 3 dimensional matrix such that R=expm(crossmat(w))

w=[R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]/2;
s=norm(w);
if(s)
    w=w/s*atan2(s,(trace(R)-1)/2);
end
end

%R2W	Rotation matrix to vector
%
%	Usage : [w]=R2w(R)
%
%	R is the 3x3 rotation matrix
%       w : 3 dimensional matrix such that R=expm(crossmat(w))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%            ORIGINAL CODE BY Prof. VENU MADHAV GOVINDU             %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [w]=R2w(R)
% 
% %R2W	Rotation matrix to vector
% %
% %	Usage : [w]=R2w(R)
% %
% %	R is the 3x3 rotation matrix
% %       w : 3 dimensional matrix such that R=expm(crossmat(w))
% 
% 
% 
% Omega=acos((trace(R)-1)/2);
% 
% n=-[R(2,3)-R(3,2),R(3,1)-R(1,3),R(1,2)-R(2,1)];
% n=n/2;
% n=div(n,sin(Omega));
% % changed from n=n/sin(Omega) to remove division by 0
% 
% w=n*Omega;
% 
% %verify=R-expm(crossmat(w))
% 
%%%%-------------------------------------------------------------------%%%%
% 
% function y=div(a,b)
% 
% %DIV	Division (division by zero is zero)
% %
% %	Usage: y=div(a,b)
% 
% mask=b==0;
% b=b+mask;
% 
% y=a./b;
% 
% y=y.*(1-mask);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%