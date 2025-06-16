

clc;
clear

D = 6;
D_1 = 2;

%27 degree PIF4-OE2 negative feedback
p_max = 85;      % day/max
p_min = 0.95;    % min
p_max0 = 100;    % night/max


% Define the x values
x_values = linspace(0, 24, 500);

% Calculate the corresponding y values using the negative_feedback function
y_values = arrayfun(@(x) negative_feedback(x,D_1,D,p_max,p_min,p_max0), x_values);




hold on;
plot(x_values,y_values,':','LineWidth', 2.5,'Color','red')
%plot(x_values,y_values,'LineWidth', 2.5,'Color','black')
hold on;


ax = gca;
% Set the font size of the x-axis and y-axis tick labels
ax.XAxis.FontSize = 20; % Set x-axis tick label font size to 14
ax.YAxis.FontSize = 20; % Set y-axis tick label font size to 14

pbaspect([1 1.2 1])
rectangle('position',[0,0,24,70],'EdgeColor','black','lineWidth',2)
ax = gca;
ax.TickLength = [0.02 0.035];
ax.LineWidth = 2;
set(gca,'TickDir','out');
axis tight;


xlim([0 24])
xticks([0 8 16 24])
yticks([0 25 50 70])
ylim([0 70])
% Save the x and y values to a text file
data = [x_values' y_values'];  % Combine x and y values into a matrix
save('negative_feedback_data_OE_27.txt', 'data', '-ascii', '-tabs');  % Save as a text file
writematrix(data, 'negative_feedback_data_OE_27.xlsx');

% Negative feedback function

function p_0 = negative_feedback(x,D_1,D,p_max,p_min,p_max0)

t0 = mod(x,24);

if t0 <= D_1 
    
    p_0 = -((p_max-p_min)/D_1)*t0+p_max ;
   
end

if t0 > D_1 && t0 <= D
    
    p_0 = p_min;
    
end

if t0 > D 


    t0 = t0-D;
    p_0 = ((p_max0-p_min)/(24-D))*t0+p_min ;
    
  
end


end
