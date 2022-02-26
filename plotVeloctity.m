function plotVeloctity()
% Read all .dat files in directory
Files = dir('*-*.dat'); % find all files that have the - and .dat
for i = length(Files):-1:1 % go through all files
    if ~isempty(strfind(Files(i).name,'conv')) %if it matches 'Conv' 
        Files(i) = []; %remove line from matrix
    end
end
clc; % clear command window
disp('Detected Files: ');
for i = 1:length(Files) %print out remaining files 
    disp([ '   ' num2str(i) ': ' Files(i).name]);
end

%%

for i = 1:length(Files)
    tmp  = csvread(Files(i).name); %set the contents of the file to a matrix tmp
    disp(Files(i).name); % display the anme of file

    x    = tmp(:,1); % Pulls the 1st column as a vectors from the C printout
    y    = tmp(:,2); % Pulls the 2nd column as a vectors from the C printout
    Temp = tmp(:,3); % Pulls the 3rd column as a vectors from the C printout
    xvec = unique(x); %Gets values of x-axis (unique values in vector x)
    yvec = unique(y); %Gets values of y-axis (unique values in vector y)
    nR   = length(yvec); % Sets the grid size in Y
    nC   = length(xvec); % sets the grid size in X
    tgrd = reshape(Temp,[nR nC]); % Reshapes vector to size [nR nC]
    
    xMin = min(xvec);
    xMax = max(xvec);
    yMin = min(yvec);
    yMax = max(yvec);
    
    figure
    
    imagesc(xvec,yvec,tgrd); % Plots the three vectors
    grid on;
    grid minor;
    hold on; % In order to add contour lines to current plot
    [C,h]= contour(xvec,yvec,tgrd, 'LineColor', [0 0 0]); %adds the black lines
    clabel(C,h); %labels the lines with numbers
    hold off; % Tell matlab that we are done with the plot
    box on; %removes the top and right border lines
    axis xy image; %keeps the aspect ratio to 1:1
    axis([xMin xMax yMin yMax]);
    title(strrep(Files(i).name,'.dat','')); % adds Title
    xlabel(' x','FontSize',14); %adds x-axis label
    ylabel('y','FontSize',14); %adds y-axis label
    colormap('jet'); %sets color style 
    colorbar; %adds the color scale on the legend
    
    %set(gca,'LooseInset',get(gca,'TightInset'))
    set(gca,'OuterPosition',[0 0 1 1]);
    pause(0.1);


    %p = get(gca,'Position');
   %  width=p(3)-p(1);
    %height=p(4)-p(2);
    %set(gcf,'PaperUnits','Inches');
   % set(gcf,'PaperSize',[width,height]);
   %set(gcf,'position',[p(1) p(2) width height]);
    %set(gcf,'PaperPositionMode','manual');   
   % set(gcf,'PaperPosition',[0 0 width height]);

    % pulls apart the file dir/name and stores name to fname
    [~,fname,~] = fileparts(Files(i).name); 
    % adds _Results.png to the end of the file name and saves it
    saveas(gcf,[ fname '.png']);    
end
end