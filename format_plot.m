function format_plot(xl,yl)
    xlabel(xl,'Interpreter','Latex','FontSize',15)
    ylabel(yl,'Interpreter','Latex','FontSize',15)
    
    box on;
    set(gcf,'color','w')
    set(gca,'FontName','Times','FontSize',15)
    set(gca,'TickLabelInterpreter','latex')
    xa = gca;
    xa.TickLength = [0.025,0.025];
    xa.LineWidth = 1.5;
end