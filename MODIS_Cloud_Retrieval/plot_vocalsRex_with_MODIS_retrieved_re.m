%% Plot the VOCALS-Rex in-situ data along with a MODIS retrieved droplet size and optical depth for a single pixel


% By Andrew John Buggee

%%

function plot_vocalsRex_with_MODIS_retrieved_re(vocalsRex, modis, pixel_2Plot)


% ------------------------------------------------------------------
% ---------------- Compute Liquid Water Path -----------------------
% ------------------------------------------------------------------
LWP_vocals = vocalsRex.lwp;                 % grams of water/m^2



% ----- Compute the CDP uncertainty -----
re_uncertainty = cloud_droplet_probe_uncertainty_estimate(vocalsRex.re);



% I want a subplot with the number concentration and altitude, and the
% effective radius with altitude
nice_blue = [0 0.4470 0.741];
nice_orange = [0.8500, 0.3250, 0.0980];


figure;

% Optical depth is ALWAYS computed from cloud top to cloud base. Let's
% determine if the profile was sampled while the plane was ascending or
% descending
dz_dt = diff(vocalsRex.altitude)./diff(vocalsRex.time);

if mean(dz_dt)>0

    % The plane was ascending. So we should flip the effective radius
    % measurement, because it starts from cloud bottom. The optical depth
    % vector starts from cloud top.

    errorbar(fliplr(vocalsRex.re), vocalsRex.tau, fliplr(re_uncertainty), 'horizontal','-o','Color','black', 'MarkerSize',10,...
        'MarkerFaceColor','black','LineWidth',1);

elseif mean(dz_dt)<0

    % The plane was descending. Therefore, the effective radius
    % measurements start from cloud top. The optical depth
    % vector starts from cloud top as well, so we don't need any altering

    errorbar(vocalsRex.re, vocalsRex.tau, re_uncertainty, 'horizontal','-o','Color','black', 'MarkerSize',10,...
        'MarkerFaceColor','black','LineWidth',1);

else

    error([newline, 'Somethings up with the optical depth vector. Check this against the altitude and effective radius',...
        newline])

end

set(gca,'YDir','reverse')
ylabel('$\tau$','interpreter','latex','FontSize',35);
xlabel('$r_{e}$ $$(\mu m)$$','Interpreter','latex')
title('Comparison between in-situ and MODIS retrieved $r_e$', 'Interpreter','latex')
grid on; grid minor; hold on;

% Fit a curve to the in-situ data to show the capability we are interested
% in devloping

% curve_fit_linewidth = 6;
% curve_fit_color = mySavedColors(1,'fixed');                        % Bright pink
%
% f = fit(tau', fliplr(double(vocalsRex.re))', 'smoothingspline','SmoothingParam',0.9);
% hold on;
% plot(f(tau),tau,'Color',curve_fit_color,'LineStyle',':', 'LineWidth',curve_fit_linewidth);

% Plot the z-space in meters on the right axis
yyaxis right
ylim([0, abs(vocalsRex.altitude(end) - vocalsRex.altitude(1))])
set(gca,'YColor','black')
ylabel('Altitude within cloud $(m)$', 'Interpreter','latex','FontSize',30);
yyaxis left

% Label cloud top and cloud bottom
% Create textbox
annotation('textbox',[0.029,0.865079365079366,0.051,0.077777777777778],...
    'String',{'Cloud','Top'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',[0.029,0.096825396825397,0.051,0.077777777777778],...
    'String',{'Cloud','Bottom'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',22,...
    'FitBoxToText','off');



% Plot the modis droplet estimate as a constant vertical line

% grab the MODIS LWP to plot
modis_lwp_2plot = modis.cloud.lwp(vocalsRex.modisIndex_minDist(pixel_2Plot)); % g/m^2

xl0 = xline(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(pixel_2Plot)),':',...
    ['MODIS $$r_{2.1} = $$',num2str(modis.cloud.effRadius17(vocalsRex.modisIndex_minDist(pixel_2Plot))), '$$\mu m$$'], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_blue);
xl0.LabelVerticalAlignment = 'bottom';

% Plot the MODIS optical depth estiamte as a constant horizontal line
yl0 = yline(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(pixel_2Plot)),':',...
    ['MODIS $$\tau_{2.1} = $$',num2str(modis.cloud.optThickness17(vocalsRex.modisIndex_minDist(pixel_2Plot)))], 'Fontsize',22,...
    'Interpreter','latex','LineWidth',2,'Color',nice_orange);
yl0.LabelHorizontalAlignment = 'left';





% Let's compute the mean number concentration within this cloud and print
% it on our plot

mean_Nc = mean(vocalsRex.total_Nc);

dim = [.2 .5 .3 .3];
str = ['$$< N_c >_{in-situ} = \;$$',num2str(round(mean_Nc)),' $$cm^{-3}$$',newline,'$$LWP_{in-situ} = $$',num2str(round(LWP_vocals,1)),' $$g/m^{2}$$'];
annotation('textbox',dim,'String',str,'FitBoxToText','on','Interpreter','latex','FontSize',18,'FontWeight','bold');
set(gcf,'Position',[0 0 1000 630])

% Create a Legend with only the two black curves
%legend('Vocals Rex In-situ Measurement', 'Desired Retrieval Profile', 'Interpreter','latex', 'Location','best')
legend('Vocals Rex In-situ Measurement', 'Interpreter','latex', 'Location','best')


end