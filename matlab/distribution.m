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
time = 2.0e-3*100;
%----------------------------------


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
hold on
plot(data.x,data.rho,'-k','LineWidth',2);
grid off
hold on
set(gca, 'XLim', [-XL/2,XL/2]);
set(gca, 'YLim', [0,1.1]);
set(gca,'XTick',[-0.6 -0.4 -0.2  0 0.2 0.4 0.6])
set(gca,'XTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2','0.4','0.6'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0 1.2])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('x-Position','FontSize',16)
ylabel('Density','FontSize',16)
saveas(figure(1),strcat(dirout,'rho.png'));
hold off

plot(xx,vel,'ko','MarkerSize',5);
hold on
plot(data.x,data.u,'-k','LineWidth',2);
grid off
hold on
set(gca, 'XLim', [-XL/2,XL/2]);
set(gca, 'YLim', [0,1.1]);
set(gca,'XTick',[-0.6 -0.4 -0.2  0 0.2 0.4 0.6])
set(gca,'XTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2','0.4','0.6'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0 1.2])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('x-Position','FontSize',16)
ylabel('Velocity','FontSize',16)
saveas(figure(1),strcat(dirout,'vel.png'));
hold off

plot(xx,pre,'ko','MarkerSize',5);
hold on
plot(data.x,data.P,'k-','LineWidth',2);
grid off
hold on
set(gca, 'XLim', [-XL/2,XL/2]);
set(gca, 'YLim', [0,1.1]);
set(gca,'XTick',[-0.6 -0.4 -0.2  0 0.2 0.4 0.6])
set(gca,'XTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2','0.4','0.6'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0 1.2])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('x-Position','FontSize',16)
ylabel('Pressure','FontSize',16)
saveas(figure(1),strcat(dirout,'pre.png'));
hold off

plot(xx,tem,'ko','MarkerSize',5);
hold on
plot(data.x,data.T,'-k','LineWidth',2);
grid off
hold on
set(gca, 'XLim', [-XL/2,XL/2]);
set(gca, 'YLim', [0.6,1.2]);
set(gca,'XTick',[-0.6 -0.4 -0.2  0 0.2 0.4 0.6])
set(gca,'XTickLabel',{'-0.6','-0.4','-0.2','0.0','0.2','0.4','0.6'})
set(gca,'YTick',[0 0.2 0.4 0.6 0.8 1.0 1.2])
set(gca,'YTickLabel',{'0.0','0.2','0.4','0.6','0.8','1.0','1.2'})
xlabel('x-Position','FontSize',16)
ylabel('Temperature','FontSize',16)
saveas(figure(1),strcat(dirout,'tem.png'));
hold off
