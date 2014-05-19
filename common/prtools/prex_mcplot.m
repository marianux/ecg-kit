%PREX-MCPLOT   PRTools example on a multi-class classifier plot
help prex_mcplot
echo on
gridsize(100)
            % Generate twice normally distributed 2-class data in 2D
a = +gendath([20,20]);    % data only
b = +gendath([20,20]);    % data only
            % Shift the second data by a vector [5,5]
            % and combine it with the first dataset in to A
A = [a; b+5];      
            % Generate 4-class labels
lab = genlab([20 20 20 20],[1 2 3 4]');
            % Construct a 4-class dataset A
A = prdataset(A,lab);    
A = setname(A,'4-class dataset')
            % Plot this 4-class dataset
figure
            % Make a scatter-plot of the right size
scatterd(A,'.'); drawnow; 
            % Compute normal densities based quadratic classifier  
w = qdc(A);     
            % Plot filled classification regions
plotc(w,'col'); drawnow;   
hold on;
            % Redraw the scatter-plot
scatterd(A);     
hold off
echo off
