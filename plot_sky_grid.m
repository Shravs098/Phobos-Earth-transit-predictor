function plot_sky_grid()
% PLOT_SKY_GRID  Draws a polar "sky dome" background on the current axes:
% concentric elevation rings (30/60/90 deg from zenith) and N/E/S/W
% compass labels. Convention: zenith at center, horizon at the outer
% ring; azimuth measured clockwise from North (up) through East (right),
% matching how compass headings are normally read.
theta = linspace(0,2*pi,200);
for rr = [30 60 90]
    plot(rr*sin(theta), rr*cos(theta), 'Color',[0.85 0.85 0.85], 'HandleVisibility','off');
    hold on;
end
text(0, 93, 'N', 'HorizontalAlignment','center', 'FontWeight','bold');
text(93, 0, 'E', 'HorizontalAlignment','center', 'FontWeight','bold');
text(0, -93, 'S', 'HorizontalAlignment','center', 'FontWeight','bold');
text(-93, 0, 'W', 'HorizontalAlignment','center', 'FontWeight','bold');
text(3, 62, '60°', 'FontSize',8, 'Color',[0.5 0.5 0.5]);
text(3, 32, '30°', 'FontSize',8, 'Color',[0.5 0.5 0.5]);
text(3, 2,  'zenith', 'FontSize',8, 'Color',[0.5 0.5 0.5]);
axis equal; axis([-100 100 -100 100]); axis off;
end