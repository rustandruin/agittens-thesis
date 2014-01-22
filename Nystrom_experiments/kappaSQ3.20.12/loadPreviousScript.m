disp('Would you like to save the experimental data to plot again later?')
disp('1 - Yes')
disp('2 - No')
saveyn = condInput('input responce (default=2): ','saveyn',' isequal(saveyn,1) || isequal(saveyn,2)','','2');
if isequal(saveyn,1)
    disp('Please input a file name.')
    currentDateStr=datestr(now, 'HHMM_ddmmyyyy');
    go=0;
    while go==0
        try
            saveName = input(['input filename (default=',currentDateStr,'.mat): '],'s');
            if isempty(saveName);saveName=[currentDateStr,'.mat'];end;
            save(saveName);
            go=1;
        catch
            disp('Please input a valid file name.');
        end
    end
end