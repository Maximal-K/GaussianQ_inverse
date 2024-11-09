clear
clc
close
%% 计算部分
x_dB = 0:0.5:10;
xdb=x_dB;
x = 10.^(x_dB/20);

Q_exact = qfunc(x);
finv_exact = qfuncinv(Q_exact);
semilogx(Q_exact,20*log10(finv_exact),'r:')
hold

%%


f_exact=@(t)qfunc(t);
f_exact_derivation=@(t)-exp(-t.^2./2)./sqrt(2.*pi);
iterations_2=[];
x_approximate=[];
%% 精度为-1
j=1;
for x_dB = 0:0.5:10
    a=qfunc(x);
    y=a(j);
    x0 = sqrt(-pi/2.*log(4.*y));
    for i=1:10000
        if i==1
            x_approximate(1)=x0;
        end
        x_approximate(i+1)= x_approximate(i)-...
            (f_exact(x_approximate(i))-y)./f_exact_derivation(x_approximate(i));
        if x_approximate(i+1)-x_approximate(i)<10^(-1)
            xf2(j)=x_approximate(i+1);
            iterations_3(j)=i;
            break
        else
            if i==10000
                xf2(j)=x_approximate(100);
                iterations_3(j)=i;
                break
            end
        end
        
    end
    j=j+1;
end
semilogx(Q_exact ,20*log10(xf2), 'k-x')
%% 精度为-12
j=1;
for x_dB = 0:0.5:10
    a=qfunc(x);
    y=a(j);
    x0 = sqrt(-pi/2.*log(4.*y));
    for i=1:10000
        if i==1
            x_approximate(1)=x0;
        end
        x_approximate(i+1)= x_approximate(i)-...
            (f_exact(x_approximate(i))-y)./f_exact_derivation(x_approximate(i));
        if x_approximate(i+1)-x_approximate(i)<10^(-6)
            xf1(j)=x_approximate(i+1);
            iterations_2(j)=i;
            break
        else
            if i==10000
                xf1(j)=x_approximate(100);
                iterations_2(j)=i;
                break
            end
        end
    end
    j=j+1;
end
semilogx(Q_exact ,20*log10(xf1), 'b+')
%



xlabel('y',"FontName","Times New Roman");
ylabel('20log_{10}(x)',"FontName","Times New Roman");
legend('Exact Q^{-1}(x)','Exact Q^{-1}(x) by NT with  \epsilon =10^{-1}',...
    'Exact Q^{-1}(x) by NT with {\epsilon} =10^{-12}',"FontName","Times New Roman",'FontSize', 13, 'Location', 'southwest');
set(gca, 'FontSize', 13)
ylim([0 10])
set(gca,'FontSize',13,'Fontangle','italic');
set(gca, 'FontSize', 13,"FontName",'Times New Roman');
xlim([0.001 0.15])