while(go==1)%HERE WE INPUT DATA FOR SPECIFIC PLOT TYPES
    for i=1:length(plotTypes)
        if (isequal(plotTypes(i),'1') || isequal(plotTypes(i),'2')) && isequal(mtxMethod,-1)
            disp('A new matrix is generated for each c or mu value.')
            disp('The same matrix is used for multiple runs.')
            disp('How many runs per c value would you like for your numerical experiments? ');
            [rns,~]=condInput('Input runs (default=3): ','rns','round(rns)==rns && rns>0','','3');
            go2=1;
            while(go2==1)
                disp('How would you like to generate your test matrices?')
                disp('1 - A = generateT3(m,n,mu)  (fast)')%generateT3
                disp('2 - A = Partially iterative method similar to 3 with randomness, not described in paper, working Beta (slower)')%mtxGenMethod1(m,n,mu)
                disp('3 - A = Algorithm 4 (slower)')%mtxGenMethod2(m,n,mu)
                disp('4 - A = Algorithm 5')%mtxGenMethod3(m,n,mu)
                disp('5 - A = matrices with hadamard structure, note: m must be a power of 2  (fast)')
                disp('6 - A = stacks of diagonal matrices, note: m/n must be an integer greater than 2  (fast)')
                disp('7 - Load A from a file, only for x-axis = c')
                disp('  note 1: Option 1 is an intelligent combination of options 4 and 5 and zero padding')
                disp('  which maximizes zero rows.  It will work for most matrix sizes and mu values.  The program')
                disp('  may crash if an invalid combination of m,n,mu is chosen (this is rare).')
                disp('  note 4: We pick the smallest m_s such that n/m_s>=coh and then use Alg 4 with zero padding to make an m by n matrix')
                
                mtxMethod = condInput('Input Method (default=1): ','mtxMethod','round(mtxMethod)==mtxMethod && 0<mtxMethod && 8>mtxMethod','','1');disp(' ');%%%%PICK MATRIX GENERATION TYPE
                go2=0;
                if mtxMethod==1
                    mtxM=@generateT3;
                    mtxgstr=', MatrixType=1';
                elseif mtxMethod==2
                    mtxM=@mtxGenMethod1;
                    mtxgstr=', MatrixType=2';
                elseif mtxMethod==3
                    mtxM=@mtxGenMethod2;
                    mtxgstr=', MatrixType=3';
                elseif mtxMethod==4
                    mtxM=@mtxGenMethod3;
                    mtxgstr=', MatrixType=4';
                elseif mtxMethod==5
                    mtxM=@matricesWithHadamardStructure;
                    mtxgstr=', MatrixType=4';
                    if (log2(m)-round(log2(m)))> 10^-14
                        disp('m must be a power of 2 for option 4')
                        go2=1;
                    end
                elseif mtxMethod==6
                    mtxM=@stacksOfDiagonalMatrices;
                    mtxgstr=', MatrixType=5';
                    if abs(round(m/n)-m/n)>10^-14 || round(m/n)<2
                        disp('m/n must be an integer greater than 2')
                        go2=1;
                    end
                elseif mtxMethod==7
                    disp('Displaying .mat files in current directory')
                    try
                        ls *.mat
                    catch
                        disp('No .mat files found in current directory')
                    end
                    disp('Please input a file name.')
                    disp('The file should have a variable named ''mtxIn'' which is an m by n matrix')
                    
                    go334455=1;
                    if ~isequal(xAxisType,2)
                        go334455=0;
                        go2=1;
                        disp('For this option you must have x-axis = c')
                    end
                    while go334455==1
                        try
                            tmpflnm = input(['input filename: '],'s');
                            if isempty(tmpflnm);error('inputflnm');end;
                            load(tmpflnm,'mtxIn');
                            if ~isequal(size(mtxIn),[m n])
                                disp('input matrix has wrong m and n')
                                error('wrong m and n')
                            end
                            if rank(mtxIn)<n
                                disp('input matrix rank deficient')
                                error('rank')
                            end
                            if max(max(abs(eye(n) - mtxIn'*mtxIn)))>10^-14
                                disp('Warning, input matrix does not have orthonormal columns')
                                disp('Ignore?')
                                disp('1 - Yes')
                                disp('2 - No')
                                [igErr,~]=condInput('Input answer (default=2): ','igErr','isequal(igErr,1) || isequal(igErr,2)','','2');
                                if isequal(igErr,2)
                                    error('oncols')
                                end
                            end
                            mtxInCoh=coherence(mtxIn);
                            mtxM = @(var1,var2,var3) ldMtxGen(mtxIn, mtxInCoh);
                            mtxgstr=', MatrixType=7';
                            go334455=0;
                        catch
                            disp('Please input a valid file name.');
                        end
                    end
                end
            end
        end
        if delta<0 && (isequal(plotTypes(i),'3') || isequal(plotTypes(i),'4') || isequal(plotTypes(i),'5') || isequal(plotTypes(i),'6')|| isequal(plotTypes(i),'7')|| isequal(plotTypes(i),'8'))
            [delta,~]=condInput('Input delta (default=.01): ','delta',' delta>0 && delta<1','','.01');
        end
        if isequal(plotTypes(i),'3')%DISPLAY MESSAGES
            disp(' ');
            disp('Please note that for plot type 3')
            disp('we must have 0<epsilon<1 and')
            disp('min(n, cc*log(cc/sqrt(delta))) <=c')
            disp('must hold, where')
            disp('cc=96*m*mu')
            disp('we will not plot any points where it does not hold')
        elseif isequal(plotTypes(i),'4')
            disp(' ');
            disp('Please note that for plot type 4')
            disp('we must have 0<epsilon<1')
            disp('we will not plot any points where it does not hold')
        elseif isequal(plotTypes(i),'5')
            disp(' ');
            disp('Please note that for plot type 5')
            disp('we must have 0<epsilon<1')
            disp('we will not plot any points where it does not hold')
        elseif isequal(plotTypes(i),'6')
            disp(' ');
            disp('Please note that for plot type 6')
            disp('we must have 0<epsilon<1')
            disp('we will not plot any points where it does not hold')
        elseif isequal(plotTypes(i),'7')
            disp(' ');
            disp('Please note that for plot type 7')
            disp('!!!epsilon IS SOLVED FOR NUMERICALLY!!!')
            disp('we must have 0<epsilon<1')
            disp(['we must have \delta=',num2str(delta),'+/-',num2str(delta*10^-12),'%'])
            disp('we will not plot any points where it does not hold')
            %         elseif isequal(plotTypes(i),'8')
            %             disp(' ');
            %             disp('Please note that for plot type 8')
            %             disp('!!!epsilon IS SOLVED FOR NUMERICALLY!!!')
            %             disp('we must have 0<epsilon<1')
            %             disp(['we must have \delta=',num2str(delta),'+/-',num2str(delta*10^-12),'%'])
            %             disp('we will not plot any points where it does not hold')
        end
        if ~isequal('m',xlbstr(end)) && isequal(xAxisType,2) && ((isequal(plotTypes(i),'2') || isequal(plotTypes(i),'5')) || (isequal(plotTypes(i),'7') && isequal(length(plotTypes),1))) %add info about gamma
            if isequal(xAxisType,2)
                xlbstr= 'c = \gamma m';
            end
            titleAdd=', c = \gamma m';
        end
    end
    go=go-1;
end
