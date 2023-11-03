
data=[];
data_points=1500;

for i=1:data_points

    F=1000 + 1000*rand(1);  % feed between 1000 and 2000
    S= 1000 + 500*rand(1);  % solvent between 1000 and 1500
    stages=2 + floor(6*rand(1)); %stages between 2 to 7

    percentage= percentage_crosscurrent(stages, S,F); 
    % function calling for creating the LLE model using the 
    % right angle triangle diagram and finding the final  percentage 
    % extracted 

    data(i,1)=percentage; %the first column of the data is the percentages
    data(i,2)=F;
    data(i,3)=S;
    data(i,4)=stages;

    i
end
data=sortrows(data); %sorting the data on the basis of percentages
csvwrite("crosscurrent_ML_data.txt",data) % the csv file will get saved in the same directory as the file
