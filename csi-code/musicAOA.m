%estimate theta by csi information, the labels are one a line
%h: csi information, M*1 vector, ith item is the csi of ith antenna
%l: the distance of two near label; lambda: lambda of the wave
%M: number of labels; D: number of eigenvector to throw away
%P: a 1*180 vector -- ith item represent the possibility of the angle in i
%degree (unnormalized)
function [P] = musicAOA(h)

global c channel_frequency d 
%global d theta_the l_the  delta_f  ;
M = 3;
lambda = c/channel_frequency;

Rxx = h*h'; 
[V, Val]=eig(Rxx);
if ~issorted(diag(Val))
        [~,I] = sort(diag(Val));
        V = V(:, I);
end

% Find max eigenvalue for normalization
max_eigenvalue = Val(3,3);

% normalization
for ii = 1:size(Val, 1)
    Val(ii, ii) = Val(ii, ii) / max_eigenvalue;
end
 
if Val(3,3)/Val(2,2) > Val(2,2)/Val(1,1)
    D = 2;
else
    D = 1;
end

E = V(:, 1:M-D);
P = zeros(1, 180);
E_product = E*E';

for theta=1:180
    a_theta{theta} = exp(-1i*((M:-1:1)-1)*2*pi*d/lambda*cos(deg2rad(theta)));
    P(theta)=abs(1./(a_theta{theta}*E_product*(a_theta{theta})'));
end
end

