function []=PLSPlotXY(x,y,maxSI,Groups, colors,figureName)
    figure('Name',figureName);
    for m=1:maxSI
        subplot(maxSI,1,m );
        hold on;
        for i=1:size(x,1)
            plot(x(i,m),y(i,m),colors{ mod( Groups(i),size(colors,2) ) });
        end
        title(num2str(m));
        hold off
    end
end