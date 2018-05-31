close all
clear 
clc

% add PACE path (if PACE folder is not in working directory, change this)
currDir = pwd;
addpath(genpath([currDir '\PACE']));

% fix seed for simulations
rand('state',0)

% to display plots, set this to 1
plots = 1;

%% READ DATA AND PREPARE INPUT
data = textread('cd4.txt','','delimiter',' ');
patientID = data(:,2);
nPatients = max(patientID);
cd4  = data(:,end);
time = data(:,3);

for i = 1:nPatients
    t{i} = time(patientID == i)';
    y{i} = cd4(patientID == i)';
    nMeasurements(i,1) = length(t{i});
end

% scheduled measurement days
measureDays = [0 2 7 10 14 28 56 84 168];

% compute min, max, median observations
[m1, i1] = min(nMeasurements);
[m2, i2] = max(nMeasurements);
p = prctile(nMeasurements,50);
display(['Minimum number of observations: ' num2str(m1) ' (subject ' num2str(i1) ')'])
display(['Maximum number of observations: ' num2str(m2) ' (subject ' num2str(i2) ')'])
display(['Median number of observations: '  num2str(p)])

% design plot
if plots 
    createDesignPlot(t, 0, 0, 1, 'CD4');
end

%% PLOT 4 RANDOM SUBJECTS
patientPlot = [5 15 29 40];
if plots  
    figure()
    for i = 1:length(patientPlot)
        j = patientPlot(i);
        subplot(2,2,i)    
        plot(t{j}/7,y{j},'ko','LineWidth',1.4, ...
            'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
        title(['Subject ' num2str(j)],'Fontweight','bold')
        xlabel('Time (Weeks)','Fontweight','bold')
        ylabel('CD4','Fontweight','bold')
        xlim([0 30])
        ylim([0 500])
        vline(measureDays/7,'r:')
    end
end

%% PLOT ALL MEASUREMENTS
if plots 
    figure()
    hold on
    for i = 1:nPatients    
        plot(t{i}/7,y{i},'ko','LineWidth',1.4, ...
            'MarkerSize',2,'MarkerEdgeColor','r','MarkerFaceColor','r')
    end
    title('All Measurements','Fontweight','bold')
    xlabel('Time (Weeks)','Fontweight','bold')
    ylabel('CD4','Fontweight','bold')

    vline(measureDays/7,'r:')
end

%% ESTIMATE MODELS UNDER EACH SPECIFICATION ({FVE90,FVE95,AIC,BIC} x {CV,GCV,GM})
tic           
for i = 0:2
    for j = 0:2
        p1 = setOptions('yname','CD4','regular',0,'selection_k','AIC', ...
               'bwmu_gcv',i,'bwxcov_gcv',j,'screePlot',0,'designPlot',0, ...
               'corrPlot',0,'numBins',0,'verbose','off');
        p2 = setOptions('yname','CD4','regular',0,'selection_k','BIC', ...
               'bwmu_gcv',i,'bwxcov_gcv',j,'screePlot',0,'designPlot',0, ...
               'corrPlot',0,'numBins',0,'verbose','off');
        p3 = setOptions('yname','CD4','regular',0,'selection_k','FVE', ...
               'bwmu_gcv',i,'bwxcov_gcv',j,'FVE_threshold',0.9,'screePlot',0,'designPlot',0, ...
               'corrPlot',0,'numBins',0,'verbose','off');
        p4 = setOptions('yname','CD4','regular',0,'selection_k','FVE', ...
               'bwmu_gcv',i,'bwxcov_gcv',j,'FVE_threshold',0.95,'screePlot',0,'designPlot',0, ...
               'corrPlot',0,'numBins',0,'verbose','off');
       
        for k = 1:4
            eval(['yy = FPCA(y,t,p' num2str(k) ');'])
            if i == 0
                estModels{1,j+1,k} = yy;
            elseif i == 1
                estModels{3,j+1,k} = yy;
            else
                estModels{2,j+1,k} = yy;
            end
        end
    end
end
toc

%% PLOT CORRELATION SURFACES UNDER EACH SPECIFICATION ({FVE90,FVE95,AIC,BIC} x {CV,GCV,GM})
% get number of figures
if plots 
    h =  findobj('type','figure');
    nf = length(h);
end    

% 4 figures, 3x3 subplots (for each specification)
for i = 1:3
for j = 1:3
for k = 1:4
    yy = estModels{i,j,k};
    nEigFn(i,j,k) = getVal(yy,'no_opt');
    t1 = getVal(yy,'out21');        
    xcorr = getVal(yy,'xcorr');
    
    if plots 
        f = figure(nf+k);
        subplot(3,3,(i-1)*3 + j)
        
        % if enough principal components, plot correlation surface
        if nEigFn(i,j,k) > 1        
            mesh(t1/7,t1/7,xcorr);
            view(0,90)  
            axis tight
        else % otherwise leave empty and write "Too Few PCs"
            text(55/7,100/7,'Too few PCs')
            xlim([0 200/7])
            ylim([0 200/7])            
        end
        xlabel('t (in weeks)','Fontweight','bold')        
        set(f, 'Position', [100, 100, 1049, 895]);            
        
        % write headings for top and left most plots
        subplot(3,3,1)
        title('CV','Fontweight','bold')
        ylabel('CV         ','Fontweight','bold','rot',0)
        subplot(3,3,2)
        title('GCV','Fontweight','bold')
        subplot(3,3,3)
        title('GM','Fontweight','bold')
        subplot(3,3,4)
        ylabel('CV         ','Fontweight','bold','rot',0)
        subplot(3,3,4)
        ylabel('GCV         ','Fontweight','bold','rot',0)
        subplot(3,3,7)
        ylabel('GM         ','Fontweight','bold','rot',0)
    end
end
end
end   

%% PLOT MEAN AND VARIANCE FUNCTIONS
for i = 1:3
    yy = estModels{i,1,1};
    t1 = getVal(yy,'out1');        
    mu(i,:) = getVal(yy,'mu');
    yy = estModels{1,i,1};
    t2 = getVal(yy,'out21');
    xcov(:,i) = diag(getVal(yy,'xcov'));   
end

if plots
    f = figure();
    set(f, 'Position', [100, 100, 900, 350]);
    
    % plot mean function
    subplot(1,2,1)
    plot(t1/7,mu,'LineWidth',1.5)
    axis tight
    grid on
    legend('CV','GCV','GM','Location','NorthWest')
    xlabel('Time (in Weeks)','Fontweight','bold')
    ylabel('CD4 Count Mean','Fontweight','bold')
    title('Mean Function','Fontweight','bold')

    % plot variance function
    subplot(1,2,2)
    plot(t2/7,xcov,'LineWidth',1.5)
    axis tight
    grid on
    legend('CV','GCV','GM','Location','NorthWest')
    xlabel('Time (in Weeks)','Fontweight','bold')
    ylabel('CD4 Count Variance','Fontweight','bold')
    title('Variance Function','Fontweight','bold')
end

%% COMPUTE CV-CV WITH 3 EIGENFUNCTIONS
p1 = setOptions('yname','CD4','regular',0,'selection_k',3, ...
                'bwmu_gcv',0,'bwxcov_gcv',0,'screePlot',0,'designPlot',0, ...
                'corrPlot',0,'numBins',0,'verbose','off','FVE_threshold',0.9);            
       
tempY = FPCA(y,t,p1);
t1 = getVal(tempY,'out1');
phi = getVal(tempY,'phi');

%% PLOT EIGFENFUNCTIONS
if plots 
    f = figure();
    set(f, 'Position', [100, 100, 1200, 500]);
    
    % first eigenfunction
    subplot(1,3,1)
    plot(t1/7,phi(:,1),'k','LineWidth',1.5)
    axis tight
    grid on
    xlabel('Time (in Weeks)','Fontweight','bold')
    ylabel('${\bf \phi_1(t)}$', 'Interpreter', 'latex');
    title('First Eigenfunction','Fontweight','bold')
    ylim([-0.15 0.17])

    % second eigenfunction
    subplot(1,3,2)
    plot(t1/7,phi(:,2),'k','LineWidth',1.5)
    axis tight
    grid on
    xlabel('Time (in Weeks)','Fontweight','bold')
    ylabel('${\bf \phi_2(t)}$', 'Interpreter', 'latex');
    title('Second Eigenfunction','Fontweight','bold')
    ylim([-0.15 0.17])

    % third eigenfunction
    subplot(1,3,3)
    plot(t1/7,phi(:,3),'k','LineWidth',1.5)
    axis tight
    grid on
    xlabel('Time (in Weeks)','Fontweight','bold')
    ylabel('${\bf \phi_3(t)}$', 'Interpreter', 'latex');
    title('Third Eigenfunction','Fontweight','bold')
    ylim([-0.15 0.17])
end

%% COMPUTE CV-CV WITH 2 EIGENFUNCTIONS, PREDICTION + 95% CI
p =  setOptions('yname','CD4','regular',0,'selection_k',2, ...
                'bwmu_gcv',0,'bwxcov_gcv',0,'screePlot',0,'designPlot',0, ...
                'corrPlot',0,'numBins',0,'verbose','off','FVE_threshold',0.9);            
       
model = FPCA(y,t,p);
t1 = getVal(model,'out1');
mu = getVal(model,'mu');
phi = getVal(model,'phi');
[ypred, xi_new, xi_var] = FPCApred(model, y, t, 0);
for i = 1:nPatients
    xHat{i} = mu + xi_new(i,:)*phi';
    se{i}   = sqrt(chi2inv(0.95,2)* diag(phi * xi_var{i} *phi'));
    se2{i}   = sqrt(norminv(0.975)* diag(phi * xi_var{i} *phi'));
end

%% PLOT PREDICTIONS
patientPlot = [5 15 29 40];

if plots
    figure()
    for i = 1:length(patientPlot)
        j = patientPlot(i);
        subplot(2,2,i)    
        hold on
        plot(t{j}/7,y{j},'ko','LineWidth',1.4, ...
            'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
        plot(t1/7,xHat{j},'b','LineWidth',1.5);
        plot(t1/7,xHat{j}+se{j}','k:','LineWidth',1.5);
        plot(t1/7,xHat{j}-se{j}','k:','LineWidth',1.5);
        plot(t1/7,xHat{j}+se2{j}','k--','LineWidth',1.5);
        plot(t1/7,xHat{j}-se2{j}','k--','LineWidth',1.5);
        title(['Subject ' num2str(j)],'Fontweight','bold')
        xlabel('Time (Weeks)','Fontweight','bold')
        ylabel('CD4','Fontweight','bold')       
        grid on
        axis tight
        ylim([-50 625])
    end
end

%% OUT-OF-SAMPLE PREDICTION
removeIdx = sort(mysample(1:nPatients,4,0));
newIdx = setdiff(1:nPatients,removeIdx);

% create new estimation set with removed subjects
k = 1;
for i = 1:length(newIdx)
    newy{k} = y{newIdx(i)};
    newt{k} = t{newIdx(i)};
    k = k+1;
end

model = FPCA(newy,newt,p);

% randomly remove measurements from prediction set
k = 1; 
for i = 1:length(removeIdx)    
    n = length(y{removeIdx(i)});
    nRemove = max(floor(rand(1,1)*n),1);
    sampleIdx = sort(mysample(1:n,nRemove,0));
    newpredy{k} = y{removeIdx(i)}(sampleIdx);
    newpredt{k} = t{removeIdx(i)}(sampleIdx);
    removeIdx2 = setdiff(1:n,sampleIdx);
    removedy{k} = y{removeIdx(i)}(removeIdx2);
    removedt{k} = t{removeIdx(i)}(removeIdx2);
    k = k+1;
end

% predict new trajectories with removed data
t1 = getVal(model,'out1');
mu = getVal(model,'mu');
phi = getVal(model,'phi');
[ypred, xi_new, xi_var] = FPCApred(model, newpredy, newpredt, 0);
for i = 1:length(removeIdx)   
    xHat_new{i} = mu + xi_new(i,:)*phi';
    se_new{i}   = sqrt(chi2inv(0.95,2)* diag(phi * xi_var{i} *phi'));
    se2_new{i}  = sqrt(norminv(0.975)* diag(phi * xi_var{i} *phi'));
end

%% PLOT NEW PREDICTIONS
if plots 
    figure()
    for i = 1:length(removeIdx)
        subplot(2,2,i)    
        hold on
        plot(newpredt{i}/7,newpredy{i},'ko','LineWidth',1.4, ...
            'MarkerSize',4,'MarkerEdgeColor','r','MarkerFaceColor','r')
        plot(removedt{i}/7,removedy{i},'ko','LineWidth',1.4, ...
            'MarkerSize',4,'MarkerEdgeColor','g','MarkerFaceColor','g')
        plot(t1/7,xHat_new{i},'b','LineWidth',1.5);
        plot(t1/7,xHat_new{i}+se_new{i}','k:','LineWidth',1.5);
        plot(t1/7,xHat_new{i}-se_new{i}','k:','LineWidth',1.5);
        plot(t1/7,xHat_new{i}+se2_new{i}','k--','LineWidth',1.5);
        plot(t1/7,xHat_new{i}-se2_new{i}','k--','LineWidth',1.5);
        title(['Subject ' num2str(removeIdx(i))],'Fontweight','bold')
        xlabel('Time (Weeks)','Fontweight','bold')
        ylabel('CD4','Fontweight','bold')       
        grid on
        axis tight        
    end
end