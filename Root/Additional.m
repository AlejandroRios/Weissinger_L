%% Additional loading
clear all
clc

%ft_n_mu amandaerak
m=7;
M=m;
v=[1 2 3  4];
n=[1 2 3 4];
mu_1=[1 3 5 7];
mu=[0 1 2 3];
nn=(m+1)-n;


fi_n=zeros(size(4,4));
fi_mu=zeros(size(4,4));
f_n_mu=zeros(size(4,4));
ft_n_mu=zeros(size(4,4));

for i=1:length(n)
    for j=1:length(mu)
        fi_n(i)=(n(i)*pi)/(m+1);
        fi_mu(j)=(mu(j)*pi)/(m+1);


f_n_mu(i,j)=(2/(m+1))*sum(mu_1.*sin(mu_1.*fi_n(i)).*cos(mu_1.*fi_mu(j)));

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

%% B_v_n

fi_v=zeros(size(4,4));
fi_n=zeros(size(4,4));
fi_nn=zeros(size(4,4));
b_v_v=zeros(size(4,4));
B1=zeros(size(4,4));
B2=zeros(size(4,4));
B_v_n=zeros(size(4,4));


for i=1:length(v)
    for j=1:length(n)
        fi_v(i)=(v(i)*pi)./(m+1);
        fi_n(j)=(n(j)*pi)./(m+1);
        fi_nn(j)=(nn(j)*pi)./(m+1);  
        
             b_v_v(i,j)=(m+1)./(4*sin(fi_v(i)));

        
        B1(i,j)=(((sin(fi_n(j)))./((cos(fi_n(j))-cos(fi_v(i))).^2))*(1-((-1).^(n(j)-v(i))))./(2*(m+1)));
        B2(i,j)=(((sin(fi_nn(j)))./((cos(fi_nn(j))-cos(fi_v(i))).^2))*(1-((-1).^((m+1-n(j))-v(i))))./(2*(m+1)));  
       
        if j~=((m+1)/2)
            B_v_n(i,j)=(B1(i,j)+B2(i,j));
        else
            B_v_n(i,j)=B1(i,j);
        end

    end

end

 B_v_n(isnan(B_v_n))=0; %bizu
 b_v_v=diag(diag(b_v_v));
 
 %% La_v_mu
AR=2.99;
con=0.376;
GAMA=-45.12;
A1=zeros(size(4,4));
A2=zeros(size(4,4));
A3=zeros(size(4,4));
A4=zeros(size(4,4));
A5=zeros(size(4,4));
A6=zeros(size(4,4));
A7=zeros(size(4,4));
A8=zeros(size(4,4));
A9=zeros(size(4,4));
A10=zeros(size(4,4));
A11=zeros(size(4,4));
A12=zeros(size(4,4));
A13=zeros(size(4,4));
A14=zeros(size(4,4));
A15=zeros(size(4,4));
A16=zeros(size(4,4));
A17=zeros(size(4,4));
A18=zeros(size(4,4));
A19=zeros(size(4,4));
A20=zeros(size(4,4));
A21=zeros(size(4,4));
A22=zeros(size(4,4));
A23=zeros(size(4,4));
A24=zeros(size(4,4));
A25=zeros(size(4,4));
A26=zeros(size(4,4));
L_v_n=zeros(size(4,4));
gt_v_n=zeros(size(4,4));

for i=1:length(v)
    for j=1:length(mu)
        
eta_v(i)=cos(fi_v(i));
etat_mu(j)=cos(fi_mu(j));

A1(i)=cos(fi_v(i));
A2(j)=cos(fi_mu(j));
% dif(i,j)=A1(i)-A2(j);
A3(i,j)=(AR*(1+con))./(2*(1-abs(A1(i))*(1-con)));
A4=tand(GAMA);
A5(i,j)=2*A4;
A6(i,j)=A1(i)-A2(j);
A7(i,j)=A1(i)+A2(j);
A8(i,j)=A6(i,j).^2;
A9(i,j)=A7(i,j).^2;
A10(i,j)=A3(i,j).^2;
A11(i,j)=A1(i).^2;
A12(i,j)=1./(A3(i,j).*A6(i,j));
A13(i,j)=1+(A3(i,j).*A6(i,j).*A4);
A14(i,j)=A13(i,j).^2;
A15(i,j)=A10(i,j).*A8(i,j);
A16(i,j)=1./(A3(i,j).*A7(i,j));
A17(i,j)=A10(i,j)*A9(i,j);
A18(i,j)=1+(2*A3(i,j).*A1(i).*A4);
A19(i,j)=(1+(A3(i,j).*A1(i).*A4)).^2;
A20(i,j)=sqrt(A14(i,j)+A15(i,j));
A21(i,j)=sqrt(A14(i,j)+A17(i,j));
A22(i,j)=A21(i,j)./A18(i,j);
A23(i,j)=sqrt(A19(i,j)+(A10(i,j).*A11(i,j)));
A24(i,j)=(A5(i,j).*A23(i,j))./A18(i,j);
A25(i,j)=A12(i,j).*(A20(i,j)-1);
A26(i,j)=A16(i,j).*(A22(i,j)-1);
L_v_n(i,j)=A25(i,j)-A26(i,j)-A24(i,j);
L_v_n(isnan(L_v_n))=0 %bizu
% gt_v_n(i,j)=((A27(i,j).*ft_n_mu(i,j)))
% gt_v_n_s(i,j)=
% (-1/(2*(m+1)))
     end
end
gt_v_n=(-1/(2*(m+1)))*((L_v_n.*ft_n_mu'))'
%% avn

a_v_n=zeros(4,4);

for i=1:length(v)
    for j=1:length(n)
        
        if j==i 
            a_v_n(i,j)=(2*b_v_v(i,j))+A3(i,j)*gt_v_n(i,j);
        else
            a_v_n(i,j)=(-2*B_v_n(i,j))+A3(i,j)*gt_v_n(i,j);
        end
    end
end
avn=a_v_n;
a_v_n_inv=inv(a_v_n);

A=[1; 1; 1; 1];
% A=1
B=a_v_n_inv*A;

C=sin(fi_n');


T=((pi*AR)/(m+1))*(B(4,:)+(2*(sum(B.*C))))









