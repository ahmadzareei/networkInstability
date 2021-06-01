

function plot_graph(A,points,data)

    if size(points,1)==2
        h = plot(A,'XData',points(1,:),'YData',points(2,:));
    else
        h =  plot(A,'XData',points(1,:),'YData',points(2,:),'ZData', points(3,:));
    end
    h.LineWidth = 2.0;
    
    if ~isempty(data)
        h.EdgeCData=data;%   weightP;
        h.LineWidth= 0.1; % A.Edges.Diameter./max(A.Edges.Diameter)*20;
        h.EdgeAlpha = 1.0;
        colormap(bluewhitered(256));
        h.LineWidth= 2;
        axis equal
        box off
        axis off
        colorbar
        set(gca,'color','k');
    end    

        h.NodeLabel='';
        h.Marker = 'none';               
end