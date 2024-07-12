close all
% function [Gma,fi_n]=slope(m)
%hoja de ajuda
% clear all
% clc
% m=7;
n=[1:(m+1)/2];
v=[1:(m+1)/2];
mu=[0:(m)/2];
mu1=[1:2:(m+1)];
nn=(m+1)-n;


fi_n=(n*pi)/(m+1);
fi_v=(v*pi)/(m+1);
fi_mu=(mu*pi)/(m+1);
fi_mu1=(mu1*pi)/(m+1);
fi_nn=(nn*pi)/(m+1); 
b_v_v=(m+1)./(4*sin(fi_v));
b2_v_v=2*b_v_v;

sin_fi_n=sin(fi_n);

cos_fi_n=cos(fi_n);
cos_fi_mu=cos(fi_mu);
etat=cos_fi_mu;
cos_fi_v=cos(fi_v);
eta=cos_fi_v;

sin_mu=sin(mu.*fi_n);
cos_mu=cos(mu.*fi_mu);
sin_mu1=sin(mu1.*fi_n);
cos_mu1=cos(mu1.*fi_mu1);

S=335.5;
% AR=6;
% con=1;
b=sqrt(AR*S);
ar_v=1./((2*(1-(1-con).*eta))/(AR*(1.+con)));
corda=b./ar_v;
% enfle=0;
tan_enfle=tand(enfle);


for i=1:length(v)
    for j=1:length(n)
         B1(i,j)=(((sin(fi_n(j)))./((cos(fi_n(j))-cos(fi_v(i))).^2))*(1-((-1).^(n(j)-v(i))))./(2*(m+1)));
         B2(i,j)=(((sin(fi_nn(j)))./((cos(fi_nn(j))-cos(fi_v(i))).^2))*(1-((-1).^((m+1-n(j))-v(i))))./(2*(m+1)));  
         
         if j~=((m+1)/2)
            B_v_n(i,j)=(B1(i,j)+B2(i,j));
        else
            B_v_n(i,j)=B1(i,j);
        end

    end

end
B_v_n(isnan(B_v_n))=0;   
        

for i=1:length(n)
    for j=1:length(mu)
BB(i,j)=(cos(mu1(i)*fi_mu(j)));
AA(i,j)=mu1(j)*(sin(mu1(j)*fi_n(i)));
    end 
end
f_n_mu=(2/(m+1)).*AA*BB;

for i=1:length(n)
    for j=1:length(mu)
       
        if i~=((m+1)/2) && j~=1
        ft_n_mu(i,j)=2*f_n_mu(i,j);
        else if i==(m+1)/2 && j~=1
        ft_n_mu(i,j)=f_n_mu(i,j);
                else if i~=(m+1)/2 && j==1
                ft_n_mu(i,j)=f_n_mu(i,j);
                    else if i==(m+1)/2 && j==1
                    ft_n_mu(i,j)=f_n_mu(i,j)/2;
                        end
                    end
            end
end
    end
end

%%

for i=1:length(v)
    for j=1:length(mu)
       
        dif_eta(i,j)=eta(i)-etat(j);
        sum_eta(i,j)=eta(i)+etat(j);
        A1(i,j)=(1+ar_v(i).*dif_eta(i,j).*tan_enfle).^2;
        A2(i,j)=(ar_v(i).*dif_eta(i,j)).^2;
        A2_2(i,j)=(ar_v(i).*sum_eta(i,j)).^2;
        A3(i,j)=ar_v(i)*dif_eta(i,j);
        A4(i,j)=ar_v(i)*sum_eta(i,j);
        A5(i)=1+2.*ar_v(i).*eta(i).*tan_enfle;
        A6(i)=(1+ar_v(i).*eta(i).*tan_enfle).^2;
        A7(i)=(ar_v(i).*eta(i)).^2;
        A8(i,j)=(eta(i)^2)-(etat(j)^2);
         mat1(i,j)=(sqrt(A1(i,j)+A2(i,j)))./A3(i,j);
         mat2(i,j)=(1./(A5(i))).*((sqrt(A1(i,j)+A2_2(i,j)))./A4(i,j));
         vet1(i)=(2*tan_enfle/A5(i)).*sqrt(A6(i)+A7(i));
         mat3(i,j)=(2*etat(j))./(ar_v(i).*A8(i,j));
         
         L(i,j)=mat1(i,j)-mat2(i,j)-vet1(i)-mat3(i,j);
         L(isnan(L))=0;
       
         
    end 
end

%%Observação
gt_v_v=(-1/(2*(m+1)))*(sum((ft_n_mu.*L),2)');
gt_v_mu=(-1/(2*(m+1)))*(ft_n_mu*L')';

Bv=b2_v_v+(ar_v.*gt_v_v);

for i=1:length(v)
    for j=1:length(n)
        
    C1(i,j)=(ar_v(i)*gt_v_mu(i,j));
    end
end

Bvn=-(2.*B_v_n)+C1;

for i=1:length(v)
    for j=1:length(n)
        
if i==j
    avn(i,j)=Bv(i);
else
    avn(i,j)=Bvn(i,j);
end
    end
end
alfa=ones(((m+1)/2),1);
Gma=inv(avn) * (alfa);

% slope= ((pi*AR)/(m+1))*(Gma((m+1)/2,:)+sum(Gma(1:1:(m-1)/2).*(2*(sin(fi_n(1:1:(m-1)/2)))')))
% 
% 
figure(1)
plot(eta,Gma')
grid on
 xlabel('eta')
 ylabel('Gamma')
title('Distribuicao de circulação')

% scatter(eta,Gma)
