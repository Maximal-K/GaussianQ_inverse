%% 参数
clc
clear
syms m n
t=linspace(0.5,3,600);
interval_time=1;
%% 实际的q函数
exactQ=qfunc(t);
%%
c=((sqrt(m.^4+6.*m.^2+1)+m.^2+1)./4);                                              %含参函数c
f=sqrt(exp(1)./pi.*((sqrt(m.^4+6.*m.^2+1)+m.^2+1)./4))./(2.*...
    ((sqrt(m.^4+6.*m.^2+1)+m.^2+1)./4)+1).*exp(-(2.*...
    ((sqrt(m.^4+6.*m.^2+1)+m.^2+1)./4)+1)./(4.*((sqrt(m.^4+6.*m.^2+1)+m.^2+1)./4)).*m.^2);          %定义函数
y1=subs(f,m,t);                                                      %函数中的字母用
diff(f);                                                             %对函数进行求导
f_lb_derivation=@(m)(exp(-(m^2*((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2))/...   %得到的导函数
    ((m^4 + 6*m^2 + 1)^(1/2) + m^2 + 1))*((3896766506551243*m)/9007199254740992 ...
    + (3896766506551243*(4*m^3 + 12*m))/(36028797018963968*(m^4 + 6*m^2 + 1)^(1/2))))/...
    (2*((3896766506551243*(m^4 + 6*m^2 + 1)^(1/2))/18014398509481984 + ...
    (3896766506551243*m^2)/18014398509481984 + 3896766506551243/18014398509481984)...
    ^(1/2)*((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2)) - (exp(-(m^2*((m^4 + 6*m^2 + 1)^...
    (1/2)/2 + m^2/2 + 3/2))/((m^4 + 6*m^2 + 1)^(1/2) + m^2 + 1))*((m^2*(m + (4*m^3 + 12*m)/...
    (4*(m^4 + 6*m^2 + 1)^(1/2))))/((m^4 + 6*m^2 + 1)^(1/2) + m^2 + 1) + (2*m*((m^4 + 6*m^2 + 1)^...
    (1/2)/2 + m^2/2 + 3/2))/((m^4 + 6*m^2 + 1)^(1/2) + m^2 + 1) - (m^2*(2*m + (4*m^3 + 12*m)/...
    (2*(m^4 + 6*m^2 + 1)^(1/2)))*((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2))/((m^4 + 6*m^2 + 1)^...
    (1/2) + m^2 + 1)^2)*((3896766506551243*(m^4 + 6*m^2 + 1)^(1/2))/18014398509481984 + ...
    (3896766506551243*m^2)/18014398509481984 + 3896766506551243/18014398509481984)^(1/2))/...
    ((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2) - (exp(-(m^2*((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2))...
    /((m^4 + 6*m^2 + 1)^(1/2) + m^2 + 1))*(m + (4*m^3 + 12*m)/(4*(m^4 + 6*m^2 + 1)^(1/2)))*...
    ((3896766506551243*(m^4 + 6*m^2 + 1)^(1/2))/18014398509481984 + (3896766506551243*m^2)/...
    18014398509481984 + 3896766506551243/18014398509481984)^(1/2))/((m^4 + 6*m^2 + 1)^(1/2)/2 + m^2/2 + 3/2)^2;
i=1;
a=1;
xf_2 = sqrt(-pi/2.*log(4*qfunc(a)));                        % 设定初值 假设为元素为Q函数x=2的点
xn_2(1) = xf_2;                                             % 将初值放入数组第一个元素
xn_2(2) =(xf_2-(subs(f,m,(xf_2))-qfunc(a))./f_lb_derivation(xf_2));%进行牛顿迭代，得到第一次的迭代的值
n=2;                                                                %每次数组迭代的次数
xxf(i)=(xf_2);                                                      %
while abs(xn_2(n)-xn_2(n-1))>10^(-3)                                %相邻两次差值跳出判断条件
    xf_2=xn_2(n);                                                   %把上一个循环的结果赋给xf_2
    xn_2(n+1)=xf_2-((subs(f,m,xn_2(n)))-qfunc(a))./f_lb_derivation(xn_2(n));% 利用上次的迭代结果进行这次的迭代计算
    n=n+1;
end
x_f(i)=xn_2(n);
iteration(i)=n-2;
%% 先画原函数 暂停
figure
% pause(interval_time*5)
base_line=plot([0,10],[qfunc(a),qfunc(a)],'k:');
hold on
fir1=plot(t,y1,'r','linewidth',1);
axis([0.83 0.98 qfunc(a)-0.01 qfunc(a)+0.025])
pause(interval_time)
%%
for i=1:length(xn_2)
    f_k=f_lb_derivation(xn_2(i));
    y_0=subs(f,m,xn_2(i));
    y_tangent_line=f_k.*(t-xn_2(i))+y_0;
    if i==1
        pause(interval_time)
        scatter(xn_2(1),y_0,'g')
        %                 text(xn_2(1),y_0+0.0008,'\fontsize{10} 起始点 X_{a}(0)');
        annotation('textarrow',[0.233229166666667,0.213229166666666],[0.888274428274429,0.878274428274429],...
            'String','\fontsize{10} 起始点X_{a}(0)','Fontangle','italic')
        pause(interval_time)
        fir2=plot(t,y_tangent_line,'g');
    else
        if i==2
            pause(interval_time)
            %             dim1=[0.7 0.63 0.05 0.1];
            %             annotation('ellipse',dim1)
            %             annotation('arrow',[0.64 0.71],[0.71 0.71])
            %             text(11,9,'\fontsize{18}R=0');
            first = plot([xn_2(i) xn_2(i)],[0 1],'b');
            pause(interval_time)
            delete(first);
            scatter(xn_2(i),y_0,'b') %text(xn_2(i),y_0,'o','color','b')
            pause(interval_time)
            %             text(xn_2(i),y_0+0.001,'\fontsize{10}X_{a}(1)','Fontangle','italic');
            annotation('textarrow',[0.726,0.726],[0.385,0.36],...
                'String','\fontsize{10} 第一次迭代X_{a}(1)','Fontangle','italic')
            pause(interval_time)
            fir3=plot(t,y_tangent_line,'b');
        else
            if i==3
                pause(interval_time)
                second=  plot([xn_2(i) xn_2(i)],[0 1],'m');
                pause(interval_time)
                delete(second);
                scatter(xn_2(i),y_0,'m')
                %text(xn_2(i),y_0-0.0015,'\fontsize{10}X_{a}(2)','Fontangle','italic');
                annotation('textarrow',[0.72,0.736],[0.333,0.345],...
                    'String','\fontsize{10} 第二次迭代X_{a}(2)','Fontangle','italic')
                pause(interval_time)
                fir4=plot(t,y_tangent_line,'c');
            else
                if i==4
                    pause(interval_time)
                    third =  plot([xn_2(i) xn_2(i)],[0 1],'k');
                    pause(interval_time)
                    delete(third);
                    pause(interval_time);
                    scatter(xn_2(i),y_0,'k')
                    %                     text(xn_2(i)+0.0004,y_0+0.001,'\fontsize{10}X_{n}(3)','Fontangle','italic');
                    annotation('textarrow',[0.76,0.74],[0.355,0.345],...
                        'String','\fontsize{10} 第三次迭代X_{a}(3)','Fontangle','italic')
                    pause(interval_time)
                    fir5=plot(t,y_tangent_line,'g');
                else
                    %                     plot(t,y_tangent_line,'c')
                    %                     text(xn_2(i),y_0,'o','color','c')
                    %                     plot([xn_2(i),xn_2(i)],[0,1],'c:')
                    %                     axis([1.9 2 qfunc(a)-0.006 qfunc(a)+0.006])
                end
            end
        end
    end
end
%% 画画
legend([fir1 base_line],{'Q_{LB-CL}_ (x)','Y=Qfunc(2)'},...
    'Location','northeast','FontSize',13,"FontName","Times New Roman")
set(gca,'FontSize',13,'Fontangle','italic');
set(gca, 'FontSize', 13,"FontName",'Times New Roman');



