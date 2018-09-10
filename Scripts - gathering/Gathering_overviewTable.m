%some preliminary data
input0  = '..\..\..\Cycle sensor data\2. Dati sperimentali\2. Dati finali\';
list0   = dir(input0); list0(1:2) = [];  list0 = list0([list0.isdir]);
exp_meth = [{'staining'},{'no stain'},{'staining'},{'no stain'},{'no stain'}];
acq_meth = [{'ImageStream'},{'ImageStream'},{'ImageStream'},{'Microscopy'},{'ImageStream'}];
ana_meth = [{'IDEAS'},{'IDEAS'},{'MATLAB'},{'MATLAB'},{'MATLAB'}];


%the real core
finalCells =  cell(0,11);
for aOUT=1:length(list0)
    input1  = [input0 list0(aOUT).name '/'];
    list1   = dir(input1);  list1(1:2) = [];
    fNumber = length(list1);

    tFname  = cell(fNumber,1);  tCond = tFname;     tRepl = tFname; tTPnumb = tFname;
    tGrRate = tFname;           tODextr = tFname;   tImNr = tFname;
    for bOUT=1:fNumber
        input2      = [input0 list0(aOUT).name '/' list1(bOUT).name];
        tFname(bOUT)  = {list1(bOUT).name};
        tCond(bOUT)   = {list1(bOUT).name(end-8:end-6)};
        tRepl(bOUT)   = {list1(bOUT).name(end-4)};
        
        x= readtable(input2);
        tGrRate(bOUT) = {round(unique(x.growth_rate),4)};
        tImNr(bOUT)   = {sum(x.nucleiNr>0)};
        if aOUT==1 || aOUT==3
            tTPnumb(bOUT) = {length(unique(x.tPoint))};
            tODextr(bOUT) = {'variable'};
        else                  
            tTPnumb(bOUT) = {NaN};
            tODextr(bOUT) = {unique(x.OD_extraction)};
        end
        disp([aOUT bOUT])
    end
    tExpMeth = repelem(exp_meth(aOUT),fNumber)';
    tAcqMeth = repelem(acq_meth(aOUT),fNumber)';
    tAnaMeth = repelem(ana_meth(aOUT),fNumber)';
    
    if aOUT==1 || aOUT==3 %dye experiment
        run('..\Scripts - auxiliary\aux_dye_rateCalc.m'); clearvars dyeOD dyeODturbi grVec16 nonExistant
        check1  = not(grFlag16);
        check2  = cell2mat(tImNr)>=10000;
    else
        run('..\Scripts - auxiliary\aux_noD_rateCalc.m'); clearVars grVecNoD %grFlagNoD
        if aOUT==2
            grFlagNoD(17)   = [];   grFlagNoD(11) = [];
            check1      = not(grFlagNoD);
            check2      = cell2mat(tImNr)>=10000;
        elseif aOUT==4
            check1      = not(grFlagNoD);
            check2      = cell2mat(tImNr)>=300;
        else
            grFlagNoD(17)   = [];
            check1      = not(grFlagNoD);
            check2      = cell2mat(tImNr)>=10000;            
        end
    end
    tCheck = num2cell(check1 & check2);
    
    finalRows   = [tFname, tCond, tRepl, tExpMeth, tAcqMeth, tAnaMeth, tTPnumb, tGrRate, tODextr, tImNr, tCheck];
    finalCells  = [finalCells;finalRows]; %#ok<AGROW>
    disp ([num2str(aOUT) ' completed'])
end

newTitles = [{'file_name'},{'condition'},{'replica'},{'exp_method'},{'acquision_tool'},{'analysis_method'},...
             {'N_tPoints'},{'growth_rate'},{'OD_extraction'},{'N_images'},{'usable_set'}];
tabEND = cell2table(finalCells,'VariableNames',newTitles);
writetable(tabEND,[input0 'overviewTable.txt']);
