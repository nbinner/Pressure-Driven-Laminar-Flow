function plotVorticity()
% Read all .dat files in directory
Files = dir('*Vorticity*.dat'); % find all files that have the - and .dat
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
    vor  = csvread(Files(i).name); %set the contents of the file to a matrix vel
    disp(Files(i).name); % display the anme of file

    X = vor(:,1);
    Y = vor(:,2);
    U = vor(:,3);
    V = vor(:,4);
    xvec = unique(X); %Gets values of x-axis (unique values in vector x)
    yvec = unique(Y); %Gets values of y-axis (unique values in vector y)
    nR   = length(yvec); % Sets the grid size in Y
    nC   = length(xvec); % sets the grid size in X
    %%tgrd = reshape(plotVorticity,[nR nC]); % Reshapes vector to size [nR nC]
    
    xMin = min(xvec);
    xMax = max(xvec);
    yMin = min(yvec);
    yMax = max(yvec);
    
    
    quiver(X,Y,U,V,1.5) 
    grid on;
    box on;
    axis xy image; %keeps the aspect ratio to 1:1
    axis([xMin xMax yMin yMax]);
    title(strrep(Files(i).name,'.dat','')); % adds Title
    xlabel(' x','FontSize',14); %adds x-axis label
    ylabel('y','FontSize',14); %adds y-axis label
    colormap('jet');
    grid minor;
    axis equal
    hold off
     [~,fname,~] = fileparts(Files(i).name); 
    % adds _Results.png to the end of the file name and saves it
    saveas(gcf,[ fname '.png']);   
end
end