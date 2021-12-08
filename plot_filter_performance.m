function [f1,f2] = plot_filter_performance(...
    x_filt,P_filt,x_true,ts,...
    x_labels,x_units,name)
    %plot_filter_performance Plots of filter confidence and accuracy
    % 
    % Inputs:
    %   x_filt: n-by-N array of filtered state vectors
    %   P_filt: x-by-n-by-N array of state estimate covariance matrices
    %   x_true: n-by-N array of truth data
    %   ts: N-length vector of timestamps
    %   x_labels: n-length array of labels for state vector variables
    %   x_units: n-length array of units for state vector variables
    %   name: string label for the dataset

    plot_2sigma = @(x,y,sigma) set(fill([x,fliplr(x)],...
        [y+2*sigma,fliplr(y-2*sigma)],'b'),...
        'facealpha',0.5,'edgealpha',0);

    n = size(x_filt,1);
%     N = size(x_filt,2);
    xlims = [ts(2),ts(end)];
    f1 = figure('Units','inches','Position',[0,0,12,16]);
    f2 = figure('Units','inches','Position',[0,0,12,16]);
    f3 = figure('Units','inches','Position',[0,0,12,16]);
    for i = 1:n
        figure(f1);
        subplot(n,1,i)
        sigma = sqrt(squeeze(P_filt(i,i,:)))';
        plot_2sigma(ts,x_filt(i,:),sigma)
        hold on
        plot(ts,x_filt(i,:),'r')
        plot(ts,x_true(i,:),'k')
        hold off
        ylabel(sprintf('%s [%s]',x_labels(i),x_units(i)),...
            "Interpreter","latex")
        legend("$2\sigma$",...
            sprintf('%s Filter Estimate',x_labels(i)),...
            sprintf('%s Truth',x_labels(i)),...
            "Interpreter","latex","Location","eastoutside")
        xlim(xlims)
        figure(f2);
        subplot(n,1,i)
        plot(ts,2*sigma,'k','LineWidth',1)
        ylabel(sprintf('%s [%s]',x_labels(i),x_units(i)),...
            "Interpreter","latex")
        xlim(xlims)
        figure(f3);
        subplot(n,1,i)
%         plot_2sigma(ts,zeros(1,N),sigma)
        plot(ts,2*sigma,'b--','LineWidth',1)
        hold on
        plot(ts,abs(x_filt(i,:)-x_true(i,:)),'k')
        hold off
%         set(gca,'YScale','log')
        ylabel(sprintf('%s Error [%s]',x_labels(i),x_units(i)),...
            "Interpreter","latex")
        xlim(xlims)
    end
    figure(f1)
    xlabel("Time [s]","Interpreter","latex")
    sgtitle(sprintf('%s State Estimate',name),...
        "Interpreter","latex")
    figure(f2)
    xlabel("Time [s]","Interpreter","latex")
    sgtitle(sprintf('%s Estimate $2\\sigma$ Error',name),...
        "Interpreter","latex")
    figure(f3)
    xlabel("Time [s]","Interpreter","latex")
    sgtitle(sprintf('%s State Estimate Error',name),...
        "Interpreter","latex")

%     set(f1,"Visible","on")
%     set(f2,"Visible","on")
%     set(f3,"Visible","on")

end