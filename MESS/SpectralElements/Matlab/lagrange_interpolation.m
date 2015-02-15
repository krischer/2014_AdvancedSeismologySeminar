
clear 
close all



fs=10; % Fontsize
lw=1. % LineWidth

% Initialize arbitrary function
n = 5
a = -[ .5 1 -3 -2 -5 4] ;
x=-1:.001:1;

b=1;
f=x*0+b;
for i=1:n
f = f + sin(pi/a(i)*x);
end 



% Interpolate this function with Lagrange Polynomials of various degree,
% calculate the misfic energy


ii=1;

for N = 2:4:6,
    
    [xi,w] = gll(N);
    
    
    
    fi=xi*0+b;
    for i=1:n
    fi = fi+ sin(pi/a(i)*xi);
    end 

    for i=1:length(x);
        for j=0:N,
            lp(j+1,i)=lagrange(N,j,x(i));
        end
    end

    s=x*0;
    for i=1:N+1
    %plot(x,lp(i,:)*fi(i),'k-')
    
    s=s+lp(i,:)*fi(i);
    end

    
    subplot(2,1,ii)
    ii=ii+1;
    
    plot(x,s,'k--','LineWidth',lw)
    set(gca,'FontSize',fs,'FontName','Times New Roman')   
    hold on
    plot(x,f,'k-','LineWidth',lw)
            
            plot(xi,fi,'s','LineWidth',lw,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','b',...
                'MarkerSize',5)

%             h = legend(' u(x) ',' $\overline{u}$(x) ',...
%                 'Location','NorthWest')
%             set(h,'Interpreter','Latex');

    
    text(-.9,3,sprintf(' N = %i ',N),'FontSize',fs,'FontName','Times New Roman')
    
    xlabel('$\xi$','Interpreter','Latex')
    ylabel('$u^e (\xi)$','Interpreter','Latex')

hold off
axis square
%set(gca,'YTick',[])
    
end

print -dpng fig_lpint