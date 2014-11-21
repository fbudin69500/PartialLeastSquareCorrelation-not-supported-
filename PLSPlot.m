function []=PLSPlot(x,maxSI,Groups, colors,figureName)
figure('Name',figureName);
    for m=1:maxSI-1
        for n=m+1:maxSI
            subplot(maxSI-1,maxSI-1,(n-1)+(m-1)*(maxSI-1) );
            hold on;
            for i=1:size(x,1)
                plot(x(i,m),x(i,n),colors{ mod( Groups(i),size(colors,2) ) });
            end
            title([num2str(m) ' vs ' num2str(n)]);
            hold off
        end
    end
end