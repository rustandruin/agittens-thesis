%License
%     Copyright (C) 2012  Thomas Wentworth
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%------------- BEGIN CODE --------------

function [ outq4rqfq3f3f34f, outq4rqfq3f3f34fStringq4rqfq3f3f34f ] = condInput( questionStringq4rqfq3f3f34f, varNameq4rqfq3f3f34f, inputConditionsq4rqfq3f3f34f, extraVarsq4rqfq3f3f34f, defaultValueq4rqfq3f3f34f )
% %NOTE: SEE COMMENTS FOR READABLE CODE
% %CONDINPUT This function intakes an input from the user and then checks
% %that the input meets certain requierements.  If it does not then it will
% %ask again until the user wises up and inputs something acceptable.
% %
% %
% %We use funny variable names beause the code allows the user to input his own variables and variable names.
% %We do not want the user defined variables to interfere with the condInput's variables
% %
% %SEE CODE HERE:
% %
% %
% %function [ out, outString ] = condInput( questionString, varName, inputConditions, extraVars, defaultValue )
% %%out is the value of the variable you want to get from the user
% %%outString is the string the user gave to define out
% %%
% %%questionString (string) is the question we ask the user
% %%varName (string) is the name of the variable you are inputting
% %%inputConditions is a string of conditions to check.  If these are not met we ask the question again
% %%extraVars(string) is a string where you define extra variables to be used.  Please
% %%terminate statements with ';'.  We will execute this string
% %%defaultValue(string) this is the default value.  If the user inputs nothing we will use this value
% %go=1;
% %if nargin>3 && ~isempty(extraVars)
% %    eval(extraVars);
% %elseif nargin>5 || nargin<3
% %    error('you must have 3-5 string inputs (questionString, varName, inputConditions, extraVars, defaultValues)');
% %end
% %while isequal(go,1)
% %    try
% %    outString=input(questionString,'s');    
% %    if isempty(outString) && isequal(nargin,5)
% %        out = eval(defaultValue);
% %        outString=defaultValue;
% %        disp(['Using default value, ',varName,' = ',defaultValue])
% %        eval([varName,' = out;']);
% %        eval(['if ',inputConditions,'; go=0;else;disp(''You must Satisfy the conditions:'');disp(''',inputConditions,''');end;'])
% %    elseif isempty(outString) && nargin<5
% %        disp('Please input a value')
% %    else
% %        out = eval(outString);
% %        eval([varName,' = out;']);
% %        eval(['if ',inputConditions,'; go=0;else;disp(''You must Satisfy the conditions:'');disp(''',inputConditions,''');end;'])
% %    end
% %    catch
% %        disp('You must Satisfy the conditions:');
% %        disp(inputConditions);
% %    end
% %end
% %eval(['out=',varName,';']);
% %end
% %
% %


goq4rqfq3f3f34f=1;
if nargin>3 && ~isempty(extraVarsq4rqfq3f3f34f)
    eval(extraVarsq4rqfq3f3f34f);
elseif nargin>5 || nargin<3
    error('you must have 3-5 string inputs (questionStringq4rqfq3f3f34f, varNameq4rqfq3f3f34f, inputConditionsq4rqfq3f3f34f, extraVarsq4rqfq3f3f34f, defaultValueq4rqfq3f3f34fs)');
end
while isequal(goq4rqfq3f3f34f,1)
    try
    outq4rqfq3f3f34fStringq4rqfq3f3f34f=input(questionStringq4rqfq3f3f34f,'s');
    if isempty(outq4rqfq3f3f34fStringq4rqfq3f3f34f) && isequal(nargin,5)
        outq4rqfq3f3f34f = eval(defaultValueq4rqfq3f3f34f);
        outq4rqfq3f3f34fStringq4rqfq3f3f34f=defaultValueq4rqfq3f3f34f;
        disp(['Using default value, ',varNameq4rqfq3f3f34f,' = ',defaultValueq4rqfq3f3f34f])
        eval([varNameq4rqfq3f3f34f,' = outq4rqfq3f3f34f;']);
        
        eval(['if ',inputConditionsq4rqfq3f3f34f,'; goq4rqfq3f3f34f=0;else;disp(''You must Satisfy the conditions:'');disp(''',inputConditionsq4rqfq3f3f34f,''');end;'])
        
    elseif isempty(outq4rqfq3f3f34fStringq4rqfq3f3f34f) && nargin<5
        disp('Please input a value')
    else
        outq4rqfq3f3f34f = eval(outq4rqfq3f3f34fStringq4rqfq3f3f34f);
        eval([varNameq4rqfq3f3f34f,' = outq4rqfq3f3f34f;']);
        eval(['if ',inputConditionsq4rqfq3f3f34f,'; goq4rqfq3f3f34f=0;else;disp(''You must Satisfy the conditions:'');disp(''',inputConditionsq4rqfq3f3f34f,''');end;'])
    end
    catch
        disp('You must Satisfy the conditions:');
        disp(inputConditionsq4rqfq3f3f34f);
    end
end
eval(['outq4rqfq3f3f34f=',varNameq4rqfq3f3f34f,';']);
end

