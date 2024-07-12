function [n,v,mu,mu1,nn]=integrador(m)
n=[1:(m+1)/2];
v=[1:(m+1)/2];
mu=[0:(m)/2];
mu1=[1:2:(m+1)];
nn=(m+1)-n;
end