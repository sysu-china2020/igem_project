function [fitresult, gof] = createFit(x, y)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'DOX-MFI fit' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 28-Oct-2020 02:15:15 自动生成


%% Fit: 'DOX-MFI fit'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'a*log(x+b)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.915735525189067 0.792207329559554 0.959492426392903];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'DOX-MFI fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'DOX-MFI fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on


