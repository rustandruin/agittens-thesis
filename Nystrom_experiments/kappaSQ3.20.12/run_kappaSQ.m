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
function [ ] = kARiPlotter( )
%kARiPlotter - This function produces the experiments described in [our
%paper].  It can plot all of the bounds expressed in a form such that
%epsilon is given in terms of m,n,mu,c,delta.  It can also run numerical
%experiments with various test matrices an two sampling stratigies.
%
% Syntax:  [] = function_name()
%
% Inputs:
%    No inputs, simply follow the on-screen commands.
%
% Outputs:
%    No outputs, this function produces figures which can then be saved at
%    the users preference.
%
% Example: 
%
% Other m-files required: 
%
% initialInputScript.m - Part of kARiPlotter.m
% inputPlotSpecificVarsScript.m - Part of kARiPlotter.m
% inputPlotTypesScript.m - Part of kARiPlotter.m
% loadPreviousScript.m - Part of kARiPlotter.m
% plot1Script.m - Part of kARiPlotter.m
% plot2Script.m - Part of kARiPlotter.m
% plot3Script.m - Part of kARiPlotter.m
% runAllTestsScript.m - Part of kARiPlotter.m
%
% Subfunctions: 
%
% blendSample.m - Performs Algorithm 2 sampling
% coherence.m - This function returns coherence
% condInput.m - Input with conditions and default value
% dataCollector.m - Collects data from numerical experiments (opt 1 and 2)
% exactlyCSample.m - Performs Algorithm 3 sampling 
% generateT3.m - A matrix generation technique
% matricesWithHadamardStructure.m - A matrix generation technique
% mtxGenMethod1.m - A matrix generation technique
% mtxGenMethod2.m - A matrix generation technique
% mtxGenMethod3.m - A matrix generation technique
% runSampling.m - Perform the individual numerical test, get kappa
% specificCohIntRot.m - Used by mtxGenMethod1
% stacksOfDiagonalMatrices.m - A matrix generation technique
% theoBound1.m - A theoretical bound
% theoBound3.m - A theoretical bound
% theoBound4.m - A theoretical bound
% theoBound6.m - A theoretical bound
% theoBound8.m - A theoretical bound
% theoBound9.m - A theoretical bound
% theoDataCollector.m - Collect data from theoBound
% MAT-files required: Optional
%

% Author: Thomas Wentworth
% North Carolina State University
% email: twentworth@gmail.com
% Website: none
% March 2012; Last revision: 20-March-2012

clc
disp('Would you like to plot previously computed data?')
disp('1 - Yes')
disp('2 - No')
loadPrev = condInput('Input choice (default=2): ','loadPrev','loadPrev==1 || loadPrev==2','','2');
if loadPrev == 2
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    initialInputScript;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script inputs all initial variables
    %These include axis type, m,n,mu, and c.
    %We also save these variables to .kARiPlotterTemp.mat so that the user
    %does not need to enter them every time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
    %These are some initial variables that we must define
    ylbstr='\kappa(SQ)';
    titlestr='';
    titlestr2='';
    plotstr='semilogy(';
    lgndstr='legend(';
    plotstr2='plot(';
    lgndstr2='legend(';
    deltaStr='';
    titleAdd='';
    titleStr='';
    delta=-1;%we are just initializing delta here, we will ask for it later.
    yPlotmin=1;
    yPlotmax=1;
    plotStylesStore=cell(4,1);
    plotStylesStore{1}='-';
    plotStylesStore{2}=':';
    plotStylesStore{3}='-.';
    plotStylesStore{4}='--';%here we store the 4 types of plotting styles for bounds
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inputPlotTypesScript;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script inputs what we are going to plot (bounds, experiments,etc)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
    disp('Would you like to plot in Color or Black/White?')
    disp('1 - Color')
    disp('2 - Black/White')
    coloryn=condInput('Input choice (default=1): ','coloryn','coloryn==1 || coloryn==2','','1');
    bndpltcnt=0;%This variable must be set to 0 here, you can not move it up.
    disp(' ')
    disp('/\/\/\/\/\/\/\/')
    disp(' ')
    mtxgstr='';%stores how matrices are generated
    mtxMethod=-1;
    rns=-1;
    go=1;%Of course we want to go. Used in next script.
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inputPlotSpecificVarsScript;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script inputs all plot specific variables, ie, stuff that depends
    %on what the user selected in inputPlotTypesScript.
    %We also display warning messages about failure conditions.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
    disp(' ')
    disp('/\/\/\/\/\/\/\/')
    disp(' ')
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loadPreviousScript;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script allows the user to load previously computed data to plot
    %again.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
    disp(' ')
    disp('/\/\/\/\/\/\/\/')
    disp(' ')
    theoData=[];
    cohErr=0;
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    runAllTestsScript;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script runs all of the requested experiments and collects the
    %data for the plots.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%   
elseif loadPrev==1
    disp('Displaying .mat files in current directory')
    try
        ls *.mat
    catch
        disp('No .mat files found in current directory')
    end
    disp('Please input which file you would like to load')
    loadName='';
    while ~isequal(exist(loadName,'file'),2)
        loadName = input('fileName: ', 's');
    end
    load(loadName)
    loadPrev=1;%we must change this value back
end%end of load prev if statement
%DO PLOTTING
if isempty(titlestr)
    error('Nothing to plot')
end
figure(1);
%%%%%%
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot1Script;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script plots the main results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
%%%%%%
if ~isempty(titlestr2)
    figure(2);
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot2Script;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script plots the failure probabilities for experimental data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
    figure(3);
    %%%%%%
    %%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plot3Script;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This script plots figure1 and 2 in one figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%
    %%%%%%
end
if loadPrev==2 && saveyn==2%only offer to save if we aren't loading a file
    loadPreviousScript;
elseif saveyn==1 && loadPrev==2
    save(saveName);
end
end
function [A,coh] = ldMtxGen(A,coh)
%This is used in the case that the user wants to load a matrix from a file.
end