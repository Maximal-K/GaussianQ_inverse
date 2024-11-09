%% 参数
clc
clear
close
x=linspace(0,15,26);
t=10.^(x./20);



%%
Qlbkw2=exp((-2.*(3.^0.5)./pi).*t.^2)./6+exp((-3.^0.5./pi).*t.^2)./6;
%%
c=(sqrt(t.^4+6.*t.^2+1)+t.^2+1)./4;
Q_Wu_vtc=sqrt(exp(1)./pi.*c)./(2.*c+1).*exp(-(2.*c+1)./(4.*c).*t.^2);


%% 实际的q函数
exactQ=qfunc(t);


%% 画画
figure
semilogy(x,Qlbkw2,'k^-')
hold on
semilogy(x,Q_Wu_vtc,'b+-')

plot(x,exactQ,'rx-')


xlabel('SNR(dB)  x(dB)=20*log_{10}_ x',"FontName","Times New Roman")
ylabel('Q(x)',"FontName","Times New Roman")
legend('Q_{LB-KW-2}_ (x)','Q_{LB-CL}_ (x)',' exact Q',...
'Location','southwest','FontSize',13,"FontName","Times New Roman")
set(gca,'FontSize',13,'Fontangle','italic',"FontName",'Times New Roman');
% set(gca, 'FontSize', 13,"FontName",'Times New Roman');
%  axis([0 8 6*10.^(-3) 0.2])
 axis([0 15 5.*10^(-7) 2.*10.^(-1)])