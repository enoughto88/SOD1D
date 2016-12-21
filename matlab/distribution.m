%1D distribution visualization for SOD1D
%Directries are written in Linux format
%Input data file is distribution.******.dat

clear all
close all
%Preamble
set(0,'defaultAxesFontName', 'Times');
set(0,'defaultTextFontName', 'Times');
set(0,'defaultUicontrolFontSize',16);
set(0,'defaultAxesFontSize',16);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 5.5 4.5]);

%-----Parameters (Change here)-----
%Directories for input data and output figures
dirin  = '/Users/kawashima/Dropbox/SOD1D/output/distribution/';
dirout = '/Users/kawashima/Dropbox/SOD1D/figure/distribution/';
fname1 = 'distribution';
%File number
file = 20;
%Grid parameter
nx   = 100;
XL   = 1.2;
ANALYTIC = 1; %0:No,1:Yes
time = 2.0e-3*100; %If ANALYTIC = 1, time must be set
%Horizontal axis setting
xrange(1)   =-0.6; xrange(2)   = 0.6;
xTickInterval = 0.2;
%Vertical axis setting
rhorange(1) = 0.0; rhorange(2) = 1.1; 
urange(1)   = 0.0; urange(2)   = 1.1;
Prange(1)   = 0.0; Prange(2)   = 1.1;
Trange(1)   = 0.6; Trange(2)   = 1.4;
rhoTickInterval = 0.2;
uTickInterval   = 0.2;
PTickInterval   = 0.2;
TTickInterval   = 0.2;
%----------------------------------

xTick    = xrange(1):xTickInterval:xrange(2);
rhoTick  = rhorange(1):rhoTickInterval:rhorange(2);
uTick    = urange(1):uTickInterval:urange(2);
PTick    = Prange(1):PTickInterval:Prange(2);
TTick    = Trange(1):TTickInterval:Trange(2);
IntxTick = num2str(xTick','%2.1f\n');
IntrhoTick = num2str(rhoTick','%2.1f\n');
IntuTick = num2str(uTick','%2.1f\n');
IntPTick = num2str(PTick','%2.1f\n');
IntTTick = num2str(TTick','%2.1f\n');

dx = XL/nx;
xx = -XL/2+dx/2:dx:XL/2-dx/2;
%Preparation of arrays
rho = zeros(nx,1);
vel = zeros(nx,1);
pre = zeros(nx,1);
tem = zeros(nx,1);

%Average the data of fileini ~ fileend
if file<10
   fname2  = '.00000';
elseif file<100
   fname2  = '.0000';
elseif file<1000
   fname2  = '.000';
elseif file<10000
   fname2  = '.00';
else
   fname2  = '.0';
end
fname = strcat(dirin,fname1,fname2,num2str(file));
data  = dlmread([fname,'.dat']);
for i = 1:1:nx
   rho(i) = data(i,1);
   vel(i) = data(i,2);
   pre(i) = data(i,3);
end
tem = pre./rho;

data = sod_analytic(time,XL);

%2D contour map of space potential
plot(xx,rho,'ko','MarkerSize',5);
grid off
hold on
set(gca, 'XLim', xrange);
set(gca, 'YLim', rhorange);
set(gca,'XTick',xTick)
set(gca,'XTickLabel',IntxTick)
set(gca,'YTick',rhoTick)
set(gca,'YTickLabel',IntrhoTick)
xlabel('x-Position','FontSize',16)
ylabel('Density','FontSize',16)
if(ANALYTIC == 1)
   hold on
   plot(data.x,data.rho,'-k','LineWidth',2);
   legend('Numerical','Analytical')
end
saveas(figure(1),strcat(dirout,'rho.png'));
hold off

plot(xx,vel,'ko','MarkerSize',5);
grid off
hold on
set(gca, 'XLim', xrange);
set(gca, 'YLim', urange);
set(gca,'XTick',xTick)
set(gca,'XTickLabel',IntxTick)
set(gca,'YTick',uTick)
set(gca,'YTickLabel',IntuTick)
xlabel('x-Position','FontSize',16)
ylabel('Velocity','FontSize',16)
if(ANALYTIC == 1)
   hold on
   plot(data.x,data.u,'-k','LineWidth',2);
   legend('Numerical','Analytical')
end
saveas(figure(1),strcat(dirout,'vel.png'));
hold off

plot(xx,pre,'ko','MarkerSize',5);
grid off
hold on
set(gca, 'XLim', xrange);
set(gca, 'YLim', Prange);
set(gca,'XTick',xTick)
set(gca,'XTickLabel',IntxTick)
set(gca,'YTick',PTick)
set(gca,'YTickLabel',IntPTick)
xlabel('x-Position','FontSize',16)
ylabel('Pressure','FontSize',16)
if(ANALYTIC == 1)
   hold on
   plot(data.x,data.P,'k-','LineWidth',2);
   legend('Numerical','Analytical')
end
saveas(figure(1),strcat(dirout,'pre.png'));
hold off

plot(xx,tem,'ko','MarkerSize',5);
grid off
hold on
set(gca, 'XLim', xrange);
set(gca, 'YLim', Trange);
set(gca,'XTick',xTick)
set(gca,'XTickLabel',IntxTick)
set(gca,'YTick',TTick)
set(gca,'YTickLabel',IntTTick)
xlabel('x-Position','FontSize',16)
ylabel('Temperature','FontSize',16)
if(ANALYTIC == 1)
   hold on
   plot(data.x,data.T,'-k','LineWidth',2);
   legend('Numerical','Analytical')
end
saveas(figure(1),strcat(dirout,'tem.png'));
hold off
