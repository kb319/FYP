
data = readtable('home68.csv', 'Format', '%{dd/MM/yyyy HH:mm}D %f');

disp(data.Properties.VariableNames);
disp(data(1:5, :));

data.Properties.VariableNames = {'Timestamp', 'Value'};

startTime = min(data.Timestamp);
endTime = max(data.Timestamp);
fullTimeRange = (startTime:hours(1):endTime)';

fullTimeTable = table(fullTimeRange, 'VariableNames', {'Timestamp'});

mergedData = outerjoin(fullTimeTable, data, 'Keys', 'Timestamp', 'MergeKeys', true);

for i = 1:height(mergedData)
    if isnan(mergedData.Value(i))
        previousDayIndex = i - 24;
        if previousDayIndex > 0
            mergedData.Value(i) = mergedData.Value(previousDayIndex);
        end
    end
end

writetable(mergedData, 'home68_filled.csv');
%%
consecutiveNaNCount = 0;

moreThan24NaNs = false;

startNaNIndex = NaN;

for i = 1:height(mergedData)
    if isnan(mergedData.Value(i))
        consecutiveNaNCount = consecutiveNaNCount + 1;
        
        if consecutiveNaNCount == 1
            startNaNIndex = i;
        end
    else
        consecutiveNaNCount = 0;
    end
    
    if consecutiveNaNCount > 24
        moreThan24NaNs = true;
        break;
    end
end

if moreThan24NaNs
    disp('There are more than 24 consecutive NaNs in the dataset.');
    disp(['The sequence starts at: ', datestr(mergedData.Timestamp(startNaNIndex))]);
else
    disp('There are no more than 24 consecutive NaNs in the dataset.');
end
