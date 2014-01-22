go=0;
while(go==0)% input various plot types
    bndpltcnt=0;%total number of bounds plotted so far
    disp(' ')
    disp('/\/\/\/\/\/\/\/')
    disp(' ')
    
    disp('1 - Numerical Experiment with Algorithm 2 sampling method, p_k=1/m')
    title1='Algorithm 2 Experiment';
    plt1LgndStr='Algorithm 2';
    
    disp('2 - Numerical Experiment with Algorithm 3 sampling method, gamma=c/m')
    title2='Algorithm 3 Experiment';
    plt2LgndStr='Bernoulli';
    
    disp('3 - Thm. 3.2, bound1:              Theoretical bound for Algorithm 2')
    title3='Thm 3.2';
    plt3LgndStr='3.2';
    
    disp('4 - Thm. 3.5, bound4:              Theoretical bound for Algorithm 2')
    title4='Thm 3.5';
    plt4LgndStr='3.5';
    
    disp('5 - Cor. 4.3, Bound6:              Theoretical bound for Algorithm 3 Sampling, gamma=c/m')
    title5='Cor. 4.3';
    plt5LgndStr='4.3';
    
    disp('6 - Cor. 3.10, Bound3:             Theoretical bound for Algorithm 2')
    title6='Cor. 3.10';
    plt6LgndStr='3.10';
    
    disp('7 - Thm. 3.13 and Cor. 4.6 Bound8: Theoretical bound for Algorithm 2 and 3 Sampling, gamma=c/m')
    title7='';%This is handled at the end of this script as this is a special case
    plt7LgndStrA='3.13';
    plt7LgndStrB='4.6';
    plt7LgndStr='';
    
    disp(' ')
    %     disp('6 - Theoretical bound for Blendenpik Sampling via UserFriendlyTailBoundsForSumsOfRandomMatrices(old thm4.2), gamma=c/m, THIS IS IN BETA')
    disp('Please input what you would like to plot.')
    disp('ie. 134 would plot the first, third and fourth options.')
    disp('You are not allowed to plot more than 4 bounds total.')
    plotTypes = input('Input desired plots: ','s');disp(' ');
    if ~isempty(plotTypes)
        go=1;
    end
    for i=1:length(plotTypes)
        if isequal(go,1) && (isequal(plotTypes(i),'1') || isequal(plotTypes(i),'2') || isequal(plotTypes(i),'3') || isequal(plotTypes(i),'4') || isequal(plotTypes(i),'5') || isequal(plotTypes(i),'6') || isequal(plotTypes(i),'7') || isequal(plotTypes(i),'8') || isequal(plotTypes(i),'9') || isequal(plotTypes(i),'0'))
            if str2double(plotTypes(i))>7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%put in new value
                disp('Invalid plot type.  Please only enter numbers 1-7.')
                go=0;
                break;
            end
            for j=i+1:length(plotTypes)
                if isequal(plotTypes(i),plotTypes(j)) && isequal(go,1)
                    disp('Please enter each plot type only once.');disp(' ')
                    go=0;
                end
            end
        elseif isequal(go,0)
        else
            go=0;
            disp('Please enter only numbers 0-9 with no commas or spaces.');disp(' ')
        end
        if isequal(go,1) && bndpltcnt<4 && (isequal(plotTypes(i),'3') || isequal(plotTypes(i),'4') || isequal(plotTypes(i),'5') || isequal(plotTypes(i),'6') || isequal(plotTypes(i),'7'))
            bndpltcnt=bndpltcnt+1;%Here we count the number of bounds plotted%%%%%%%%%%%%%%%%%%%%%%%%%UPDATE ABOVE
        elseif isequal(go,1) && bndpltcnt==4
            go=0;
            isequal(go,1);
            disp('Please only select 4 total bounds to plot.')
        end
    end
end
%Here we make sure that the legend string for plot 7 is the one we want.
for i=1:length(plotTypes)
    if isequal(plotTypes(i),'1') || isequal(plotTypes(i),'3') || isequal(plotTypes(i),'4') || isequal(plotTypes(i),'6')
        if isempty(plt7LgndStr)
            plt7LgndStr = plt7LgndStrA;
            title7='Thm. 3.13';
        elseif ~isequal(plt7LgndStr,plt7LgndStrA)
            plt7LgndStr = [plt7LgndStrA,', ',plt7LgndStrB];
            title7='Thm. 3.13, Cor. 4.6';
        end
    end
    if isequal(plotTypes(i),'2') || isequal(plotTypes(i),'5')
        if isempty(plt7LgndStr)
            plt7LgndStr = plt7LgndStrB;
            title7='Cor. 4.6';
        elseif ~isequal(plt7LgndStr,plt7LgndStrB)
            plt7LgndStr = [plt7LgndStrA,', ',plt7LgndStrB];
            title7='Thm. 3.13, Cor. 4.6';
        end
    end
end
if isempty(plt7LgndStr)
    plt7LgndStr = [plt7LgndStrA,', ',plt7LgndStrB];
    title7='Thm. 3.13, Cor. 4.6';
end