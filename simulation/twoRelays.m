clear all
close all
clc
sigma=1;
k=10000;
no_channel=5;%
y=raylrnd(sigma,k,1);
[pr,chan]=hist(y,no_channel);
prob=pr/k;
y1=raylrnd(sigma,k,1);
[pr1,chan1]=hist(y1,no_channel);
prob1=pr1/k;
% chan=[0.001, 0.005];
% chan1=chan;
% prob=[0.5, 0.5];
% prob1=prob;
joint_prob=(prob'*prob1);
a=zeros(size(joint_prob));
b=zeros(size(joint_prob));
c=zeros(size(joint_prob));
d=zeros(size(joint_prob));
del_zero=zeros(size(joint_prob));
del_one=zeros(size(joint_prob));
ca=zeros(size(joint_prob));
cb=zeros(size(joint_prob));
pow=zeros(size(joint_prob));
%rate1=zeros(size(joint_prob));
lamda=50; % lagrange multiplier
g_sr=(1/(50)^2)*chan;
g_rd=1/(50)^2;
g_sd=(1/(70)^2)*chan1;
g_sr_ones=g_sr'*ones(1,no_channel);
g_sd_ones=ones(no_channel,1)*g_sd;
noise=10e-3;
relay_pow=10e-3;
band_width=10;
%pow=0:0.1:7;
source_pow=10e-3;
avg_pow=source_pow*ones(size(joint_prob));
total_pow=3;
omega=conj(-(-1)^(1/3));
omega2=-(-1)^(1/3);
n=1;
step_size=100;
while((total_pow-source_pow)>0.00001 || n<20)
    pow=zeros(size(joint_prob));
    for i=1:no_channel
        for j=1:no_channel
            if (((g_sr(i)*g_rd*relay_pow)/(noise^4))>(((lamda*reallog(2)/band_width)-((g_sd(j))/(noise)^2)))*(1+(g_rd*relay_pow)/(noise)^2))
%                 if(j==1)
%                   g_sr(i)
%                   g_sd(j)
%                  end   
            a(i,j)=(lamda*reallog(2)/band_width)*(((g_sr(i))^2)*(g_sd(j))/(noise^6));
            b(i,j)=((lamda*reallog(2)/band_width)*(((g_sr(i))^2)/(noise^4))*(1+(relay_pow*g_rd/(noise^2)))) + (((g_sd(j))/(noise^2))*((g_sr(i))/(noise)^2))*((lamda*2*reallog(2)/band_width)*(1+g_rd*relay_pow/(noise^2))-(g_sr(i))/(noise^2));
            c(i,j)=((lamda*reallog(2)/band_width)-((g_sd(j))/(noise)^2))*((2*(g_sr(i))/(noise)^2)+(2*(g_sr(i))*g_rd*relay_pow/(noise)^4)) + ((lamda*reallog(2)/band_width)*((g_sd(j))/(noise)^2))*((2*relay_pow*g_rd/(noise)^2)+((relay_pow^2)*(g_rd^2)/(noise)^4)+1) + ((lamda*reallog(2)/band_width)*((g_rd*relay_pow)/(noise)^2)*(1+((g_rd*relay_pow)/(noise)^2))*((g_sr(i))/(noise^2)));
            d(i,j)=((((lamda*reallog(2)/band_width)-((g_sd(j))/(noise)^2)))*(1+(g_rd*relay_pow)/(noise)^2)-((g_sr(i)*g_rd*relay_pow)/(noise^4)))*(1+(g_rd*relay_pow)/(noise)^2);
            del=18*a(i,j)*b(i,j)*c(i,j)*d(i,j)-4*(b(i,j)^3)*d(i,j)+(b(i,j)^2)*(c(i,j)^2)-4*a(i,j)*(c(i,j)^3)-27*(a(i,j)^2)*(d(i,j)^2)
            del_zero(i,j)=((b(i,j))^2)- 3*(a(i,j))*(c(i,j));
            del_one(i,j)= 2*((b(i,j))^3)-9*(a(i,j))*(b(i,j))*(c(i,j))+27*((a(i,j))^2)*(d(i,j));
            ca(i,j)=((del_one(i,j)+sqrt(((del_one(i,j))^2)-4*((del_zero(i,j))^3)))/2)^(1/3);
            cb(i,j)=((del_one(i,j)-sqrt(((del_one(i,j))^2)-4*((del_zero(i,j))^3)))/2)^(1/3);
            if(del<0)
            r=(-1/(3*a(i,j)))*(b(i,j)+omega*(ca(i,j)+cb(i,j)));
            if(abs(imag(r))<=0.00000001)
            root=real(r);
            end
            else
            r=(-1/(3*a(i,j)))*(b(i,j)+omega*(ca(i,j))+omega2*(cb(i,j)));
            if(abs(imag(r))<=0.00000001)
            root=real(r)
            end
%             if(root<0)   
%             r=(-1/(3*a(i,j)))*(b(i,j)+omega*(ca(i,j))+omega2*(cb(i,j)));
%             if(abs(imag(r))<=0.00000001)
%             root=real(r)
%             end
%             end
            end
            pow(i,j)=root;           
%            rate1(i,j)=band_width*log2(1+(g_sd(j)*pow(i,j)/((noise)^2))+((g_rd*relay_pow*g_sr(i)*pow(i,j))/(noise)^4)/(1+(g_rd*relay_pow/(noise)^2)+(g_sr(i)*pow(i,j)/(noise)^2)));
            end
        end  
    end
    pow
    total_pow=sum(sum(joint_prob.*pow))
    if (n<1000)
        lamda=lamda+(step_size)*(total_pow-source_pow)
    else
        lamda=lamda+(step_size/n)*(total_pow-source_pow)
    end
   % lamda=lamda+(step_size)*(total_pow-source_pow)
    if(lamda<0)
        lamda=0.1/n;
    end
    power(n,:)=total_pow;
    rate1(n,:)=sum(sum(joint_prob.*(band_width*log2(1+(g_sd_ones.*avg_pow/((noise)^2))+((g_rd*relay_pow*g_sr_ones.*avg_pow)/(noise)^4)./(1+(g_rd*relay_pow/(noise)^2)+(g_sr_ones.*avg_pow/(noise)^2))))));
    rate(n,:)=sum(sum(joint_prob.*(band_width*log2(1+(g_sd_ones.*pow/((noise)^2))+((g_rd*relay_pow*g_sr_ones.*pow)/(noise)^4)./(1+(g_rd*relay_pow/(noise)^2)+(g_sr_ones.*pow/(noise)^2))))));
    iteration(n,:)=n;
    n=n+1;
end
% bincenters = {g_sr,g_sd};
% hist3([0,0],bincenters);
% Hsurface = get(gca,'children');
%// apply bin count format
% pad = [0 0 0 0 0;0 1 1 0 0;0 1 1 0 0;0 0 0 0 0;0 0 0 0 0]; %// padding for each point
% powtrans=kron(reshape(pow,[no_channel,no_channel]),pad); %// apply padding to each point
% %// update plot
% set(Hsurface,'ZData',powtrans)
% %// to set colour based on bar height
% colormap('autumn');
% set(Hsurface,'CData',powtrans,'FaceColor','interp')
% xlabel('Gsr')
% ylabel('Gsd')
% zlabel('Power allocation')
% title(['Power allocation for Moderate SNR at Relay'])
% figure
% [AX,H1,H2] = plotyy(iteration,power,iteration,[rate,rate1],'plot');
% set(get(AX(1),'Ylabel'),'String','power') 
% %set(H1, 'color','b')
% set(get(AX(2),'Ylabel'),'String','rate(MHz)')
% %set(H2, 'color','r')
% legend({'power' 'power control rate' 'average power rate'}, 'Location','NorthEast')
% xlabel('iteration')  
% title('Optimal power and rate vs iteration')
plot(iteration,power)
xlabel('iteration') 
ylabel('source power') 
%legend({'bad channel optimal rate' 'good channel optimal rate' 'bad channel average power rate' 'good channel average power rate'}, 'Location','NorthWest')
title('Source power and rate')
