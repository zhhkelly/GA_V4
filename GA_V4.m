%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Genetic Algorithm                                                       % 
%                                                                         %
% for optimizing cell-scale EP model parameters                           %
% based on experimental AP morphological features                         %
%                                                                         %
% Authors: Kelly Zhang, Chelsea E. Gibbs, Patrick M. Boyle                %
% Version: 4.0.0                                                          %
% Last updated: 1/10/2021                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize UI and GA
try
    %Close any open UIs from previous runs
    delete(findall(0));
    disp("Checkpoint passed: Old UIs closed. This is the latest run.")
catch
    %If no open UIs, print message
    disp("Checkpoint passed: No old UIs open. This is the latest run.")
end

% Open UI
%.mlapp version should match name here
myApp                   = GA_UI_V4; 
uiwait(myApp.UIFigure);

% Experiment-specific Parameters
model_name              = myApp.ModelName.Value;
experiment_name         = myApp.ExperimentName.Value;
population_size         = myApp.PopulationSize.Value;
ElitePop_size           = myApp.ElitePopulationSize.Value;
MatingPool_size         = population_size-ElitePop_size;
ChildrenPop_size        = MatingPool_size;
start_gen               = myApp.StartGeneration.Value;
end_gen                 = myApp.EndGeneration.Value;
sim_dur                 = myApp.SimDuration.Value; %ms
dt                      = myApp.dt.Value;
test_run_factor         = 5;
full_run_factor         = 1;
% test_run_factor is hard coded. The value is x times the mating pool size.
% Overfill number of variants to ensure full quiescent population.
% full_run_factor is 1 since we need to run only the quiescent full pop.

% Stiumulus Protocol
stim_curr               = myApp.StimCurrent.Value;
stim_dur                = myApp.StimDuration.Value;
stim_bcl                = myApp.StimBCL.Value;
stim_num                = myApp.StimNum.Value;

% Experimental Data
sample_size             = myApp.Exp_SS.Value;

% The empty quotes are array padding to allow for consistent dimensions
% even with blanks in the UI input.
target_features         = str2double([
    myApp.CL_mean.Value,        myApp.CL_SE.Value,          "";
    myApp.BPM_mean.Value,       myApp.BPM_SE.Value,         "";
    myApp.MDP_mean.Value,       myApp.MDP_SE.Value,         "";
    myApp.Peak_mean.Value,      myApp.Peak_SE.Value,        "";
    myApp.APA_mean.Value,       myApp.APA_SE.Value,         "";
    myApp.APD20_mean.Value,     myApp.APD20_SE.Value,       "";
    myApp.APD50_mean.Value,     myApp.APD50_SE.Value,       "";
    myApp.APD90_mean.Value,     myApp.APD90_SE.Value,       "";
    myApp.NotchDepth_mean.Value,myApp.NotchDepth_SE.Value,  "";
    myApp.vmaxUp_mean.Value,    myApp.vmaxUp_SE.Value,      "";
]);

Crossover_checkbox      = myApp.Crossover_checkbox.Value;
num_Xpoints             = myApp.NumCrossoverPoints.Value;

% Preallocate struct for storing variant data
CurrentGen_struct       = struct(...
    "time_data",[],"trace_data",[],"APfeature_data",[],...
    "adj_data",[],"percent_error",[],"fitness_score",[],"variant_id",[]);

%% Initialize openCARP compatible models
% Initialize model
MakeDynamicModel_command = sprintf(...
    "/Software/cme/openCARP/bin/make_dynamic_model.sh %s", model_name);
system(MakeDynamicModel_command)

%% Make copy of model files in experiment directory
MakeModelFileCopies_command = sprintf(...
    "cp ./%s* ./%s",...
    model_name, experiment_name);
system(MakeModelFileCopies_command)

%% Read in initial population data
% Check if a starting generation .txt file exists in experiment dir.
StartGenMat_path = sprintf(...
    "./%s/generation_%d/generation_%d_adjs.mat",...
    experiment_name,start_gen,start_gen);
if exist(StartGenMat_path, "file")
    disp("Checkpoint passed: starting generation .mat exists in dir")
    start_population=load(StartGenMat_path);
    start_population=struct2array(start_population);
    ParamOrder_array=string(start_population(:,1));
    start_data=double(start_population(:,2:end));
else
% If it doesn't exist, quit the run
    disp(...
        "Checkpoint FAILED: Use the LHS.m to make generation_x_adjs.mat")
    disp("Now exiting the GA");
    return
end

%% Make the .par files specific to the experiment
% Execute the GAHereDoc.sh to make a .par file for running FULL CARP sims
ParamOrder_string           = sprintf("%s ",ParamOrder_array);
GAHereDoc_command           = sprintf(...
    "bash -u ./GAHereDoc.sh %s %s %s %d %d %d %d %d %d %d %d %s",...
    "full_runs",experiment_name,model_name,start_gen,population_size,...
    stim_bcl,stim_num,stim_dur,stim_curr,sim_dur,...
    test_run_factor,ParamOrder_string);
system(GAHereDoc_command);

%% Make the full-population-sized mesh specific to the experiment
FullMesh_path               = sprintf(...
    "./%s", experiment_name);
MakeLineMesh_test_command   = sprintf(... 
    "bash -u ./MakeLineMesh.sh %s %d %s",...
FullMesh_path, population_size,"full_runs");
system(MakeLineMesh_test_command)

%% Make the run_type directory within the generation directory
MakeRuntypeDir_command      = sprintf(...
    "mkdir ./%s/generation_%d/%s", experiment_name, start_gen, "full_runs");
system(MakeRuntypeDir_command);

%% Make the .adj files for each parameter for starting generation
% The starting population requires a full-run initialization before
% iterating through further generations
MakeAdjlFile_func(experiment_name, start_gen,"full_runs",...
    length(ParamOrder_array),ParamOrder_array, population_size,start_data);

%% Execute the RunCarpSim.sh for starting generation
% Runs each variant simulation in full to generate the vm_exploded_xxx.dats
RunCarpSim_command          = sprintf(...
    "bash -u ./RunCarpSim_V4.sh %d %s %s %s %d %d",...
    start_gen,"full_runs",model_name,experiment_name,...
    population_size,sim_dur);
system(RunCarpSim_command);

%% Begin iterating through generations
for gen_counter             = start_gen:end_gen  
    %% Consolidate voltage data from vm_exploded_xxx.dats
    for pop_counter         = 1:population_size
        % Format the variant number in the form of xxx
        variant_id          = sprintf(...
            "%03d",pop_counter-1);
        % Open the vm_exploded_xxx.dat file
        VmExploded_file     = sprintf(...
            "./%s/generation_%d/full_runs/vm_exploded_%s.dat",...
            experiment_name,gen_counter,variant_id);
        file_ID             = fopen(VmExploded_file,'r');
        format_Spec         ="%f";
        array_size          =[2 Inf];
        % save vm_exploded_xxx.dat contents as array
        VmExploded_data     = fscanf(file_ID,format_Spec,array_size)';
        fclose(file_ID);
        % save time and voltage data to tablecl
        CurrentGen_struct(pop_counter).time_data        =... 
            VmExploded_data(:,1);
        CurrentGen_struct(pop_counter).trace_data       =... 
            VmExploded_data(:,2);
    end

    %% Consolidate AP morphological feature data
    for pop_counter         = 1:population_size
        CurrentGen_struct(pop_counter).APfeature_data   =...
            APFeatureCheck_func(...
                CurrentGen_struct(pop_counter).time_data,...
                CurrentGen_struct(pop_counter).trace_data,...
                size(target_features,1));
    end
    
    %% Consolidate adjustment factors for each parameter 
    CurrGenMat          = sprintf(...
        "./%s/generation_%d/generation_%d_adjs.mat",...
        experiment_name, gen_counter, gen_counter);  
    PopAdj_data         = load(CurrGenMat);
    PopAdj_data         = struct2array(PopAdj_data);
    ParamOrder_array    = string(PopAdj_data(:,1));
    ParamAdj_array      = PopAdj_data(:,2:end);
    num_params          = length(ParamOrder_array);
    
    for pop_counter     = 1:population_size
        CurrentGen_struct(pop_counter).adj_data =...
            ParamAdj_array(:,pop_counter);
    end
    
    %% Calculate Fitness Scores
    PercentError_array  = FitnessScore_func(...
    {CurrentGen_struct.APfeature_data},target_features(:,1));
    % Sum the percent errors for each variant
    % Get the magnitude of sum(percent_errors) since any direction away
    % from an error of zero is unfavorable
    fitness_scores      = sum(PercentError_array,1,'omitnan');
    % Store percent errors and fitness scores in CurrentGen_struct
    for variant_counter = 1 : population_size
        CurrentGen_struct(variant_counter).fitness_score =...
            fitness_scores(variant_counter);
        CurrentGen_struct(variant_counter).percent_error =...
            PercentError_array(:,variant_counter);
    end
    %% Assign Variant IDs.
    % IDs are necessary to keep track of which vm_exploded_xxx.dat goes
    % with each parameter set after struct reordering in the next step
    for variant_counter = 1:population_size
        variant_id = sprintf(...
            "%03d", variant_counter);
        CurrentGen_struct(variant_counter).variant_id = variant_id;
    end
    
    %% Extract elite population
    % The elite population is the X% of high scoring variants from the full
    % population. These variants' parameters will be passed to the next
    % generation without modifications.
    
    % Sort the current generation from good to bad fitness scores
    CurrentGen_table            = struct2table(CurrentGen_struct);
    CurrentGen_struct_sorted    = table2struct(...
                                  sortrows(...
                                  CurrentGen_table,'fitness_score'));
                              
    % Preallocate array for elite population. Array contains just parameter
    % adjustment factors. Each row is a parameter, each column a variant
    EliteParams_array           = zeros (num_params,ElitePop_size);
    
    % Extract elite population and store in the array.
    % Each column is a variant. Each row is a parameter adj
    for variant_counter         = 1:ElitePop_size
        EliteParams_array(:,variant_counter) =...
            CurrentGen_struct_sorted(variant_counter).adj_data';
    end
    
%     % Output elite population scores of each generation into a text file
%     elite_scores                = [CurrentGen_struct_sorted(...
%                                    1:ElitePop_size).fitness_score];
%     ScoreEvo_path               = sprintf(...
%         "./%s/score_evolution.txt",experiment_name);
%     dlmwrite(ScoreEvo_path,elite_scores,'-append')
%  
%% Plot Elite population traces
    f1 = figure;
    for variant_counter         = 1:ElitePop_size
        plot(CurrentGen_struct_sorted(variant_counter).time_data,...
             CurrentGen_struct_sorted(variant_counter).trace_data,...
             "DisplayName",...
             num2str(CurrentGen_struct_sorted(variant_counter).variant_id));
        hold on
    end
    hold off
    legend('location','eastoutside')
    SavePlot_path = sprintf(...
        "./%s/generation_%d/PlotElite_%d",...
        experiment_name,gen_counter,gen_counter);
    saveas(gcf, SavePlot_path, 'jpeg');
    close(f1)

     %% Save .mat for Current Generation
    SaveMat_path                = sprintf(...
        "./%s/generation_%d/SortedGeneration_%d.mat",...
        experiment_name, gen_counter, gen_counter);
    save(SaveMat_path, "CurrentGen_struct_sorted")
     
    %% Plot Feature Distrubution
%     f2 = figure;
%     PlotFeatures_func(... 
%         {CurrentGen_struct.APfeature_data},target_features(:,1),...
%         target_features(:,2));
%     SavePlot_path = sprintf(...
%         "./%s/generation_%d/PlotGen_%d",...
%         experiment_name,gen_counter,gen_counter);
%     saveas(gcf, SavePlot_path, 'fig');
%     close(f2)
    %% Termination check
    % The algorithm will stop running if there is a variant with all
    % features within the tolerance of the target data.
    
    % Get maximum allowed percent error for each feature.
    TargetError_array   = abs(target_features(:,2)...
                          ./target_features(:,1))*100;
    
    % Check whether the percent error for each parameter for each variant
    % is within the allowed error.
    for variant_counter = 1:population_size
        TolCheck_array  =...
            abs(CurrentGen_struct_sorted(variant_counter).percent_error')...
            -TargetError_array; 
        max_error       = max(TolCheck_array);
        if max_error    <=0
            disp(CurrentGen_struct_sorted(variant_counter).adj_data)
            fprintf(...
                "Variant number %d is the best fitting. Quitting GA."...
                ,variant_counter);
            return
        end
    end
    disp("Checkpoint passed: no best fit yet")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Paced runs
    % External stimulus protocol drives model beating rate
    if stim_num                 ~= 0
        % Create a overfill population of variants. Whole, even number.
        TestPool_size           = round(MatingPool_size*test_run_factor);
        if mod(TestPool_size,2) ~= 0
           TestPool_size        = TestPool_size+1;
        end
        
        % Preallocate array for this overfill population.
        % Each row is a parameter adj, each column is a variant.
        TestPoolParams_array    = zeros(num_params, TestPool_size);
        
%Tournament selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Two variants are chosen from the full population. Their fitness
        % scores are compared. The better fitness score variant will be
        % placed in the test pool.
        
        for tourney_counter     = 1:TestPool_size
            % Contestants is an array that contains two numbers.
            contestants         = datasample(...
                1:population_size,2,'Replace',true);
            % The smaller number is closer to the top of the
            % CurrentGen_struct_sorted sruct, and therefore corresponds to
            % a better fit variant
            tourney_winner      = min(contestants);
            
            % Place parameters of better scoring variant in test pool.
            TestPoolParams_array(:,tourney_counter) =...
                CurrentGen_struct_sorted(tourney_winner).adj_data;
        end
        
%Pair Up Variants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Randomly order the variants numbers in the test pool
        RandomOrder_array       = datasample(...
            1:TestPool_size,TestPool_size,"replace",false);
        % Pair up the variants by reshaping RandomOrder_array
        % Each column of PairedVariants_array has the 2 parent variant nums
        PairedVariants_array    = reshape(...
                                  RandomOrder_array,[2,TestPool_size/2]);
        
%Crossover Operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single or multipoint crossover. Explanation within function.
        % This is an optional operator. If skipped, the child population
        % will be a duplicate of the test pool.
        if Crossover_checkbox                == true
            TestChildParams_array            = CrossOver_func(...
                TestPoolParams_array,PairedVariants_array,num_Xpoints);
        else
            TestChildParams_array            = TestPoolParams_array;
        end

%Mutation Operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Randomly alter parameters to introduce more parameter diversity.
        % Explanation within function.
        TestChildParams_array_mut = mutation_func(TestChildParams_array);
        
%Make next generation directory%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The mutated test child variants belong to the next generation.
        % Create new generation directory
        MakeNewGen_func(experiment_name, gen_counter+1);
        
%Run the GAHereDoc.sh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        % Execute the GAHereDoc.sh to make a .par file for TEST CARP sims
        % Test sims are non-stimulated. We want to pick out variants where
        % the spontaneous, intrinsic bcl is longer than the stim_bcl.
        % sim_dur = stim_bcl because we are looking within one stim_bcl,
        % no point of running it longer
        
       GAHereDoc_command            = sprintf(...
            "bash -u ./GAHereDoc.sh %s %s %s %d %d %d %d %d %d %d %d  %s",...
            "test_runs",experiment_name,model_name,gen_counter+1,...
            population_size,stim_bcl,stim_num,stim_dur,stim_curr,...
            stim_bcl,test_run_factor,ParamOrder_string);
        system(GAHereDoc_command);
 
%Make the TEST mesh%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TestMesh_path               = sprintf(...
            "./%s/generation_%d", experiment_name, gen_counter+1);
         MakeLineMesh_test_command = sprintf(... 
            "bash -u ./MakeLineMesh.sh %s %d %s",...
            TestMesh_path, TestPool_size,"test_runs");
        system(MakeLineMesh_test_command)
        
%Make the run_type directory within the generation directory%%%%%%%%%%%%%%%
        MakeRuntypeDir_command       = sprintf(...
            "mkdir ./%s/generation_%d/%s",experiment_name,...
            gen_counter+1,"test_runs");
        system(MakeRuntypeDir_command);

%Run the MakeAdjFiles_func for TEST RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MakeAdjlFile_func(...
            experiment_name,gen_counter+1,"test_runs",...
            length(ParamOrder_array), ParamOrder_array,TestPool_size,...
            TestChildParams_array_mut)

% Run the RunCarpSim.sh for TEST RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The shell script runs the short simulations for the test pool.
        % the final input parameter of the .sh is also sim_dur=stim_bcl
        % Further explanations of the script are commented in its .sh file
        RunCarpSim_command          = sprintf(...
            "bash -u ./RunCarpSim_V4.sh %d %s %s %s %d %d",gen_counter+1,...
            "test_runs",model_name,experiment_name,TestPool_size,stim_bcl);
        system(RunCarpSim_command)
        
%Read in the quiescent variants file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        QuiescentVariants_path      = sprintf(...
            "./%s/generation_%d/%s/quiescent_variants.txt",....
            experiment_name,gen_counter+1,"test_runs");
        QV_ID                       = fopen(QuiescentVariants_path);
        format_spec                 = "%f";
        array_size                  = [1 Inf];
        QuiescentVariants_array     = fscanf(QV_ID,format_spec,array_size);
        fclose(QV_ID);
        
%Calculate new test_run_factor%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CURRENTLY KEPT STATIC.
        % Calculate proportion of TestPool_size that was quiescent
        % Add 20% of population as buffer
        %test_run_factor = (nnz(QuiescentVariants_array)+...
            %(0.2*TestPool_size))/TestPool_size;
        
        % In case factor goes below 1.2, re-initialize
        %if test_run_factor < 1.2
            % test_run_factor = 1.2;
        %end
        
%Trim population%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preallocate array for quiescent variants' parameters
        % Each row is a parameter, each column is a variant
        TrimmedChildParams_array    = zeros(...
            num_params,ChildrenPop_size);
 
        for variant_counter         = 1:ChildrenPop_size
            % Save parameters of variants listed in QuiescentVariants_aray
            % in the TrimmedChildParams_array.
            TrimmedChildParams_array(:,variant_counter)=...
                TestChildParams_array_mut(...
                :,QuiescentVariants_array(variant_counter));
        end
        
% Make the .par file for FULL runs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    GAHereDoc_command           = sprintf(...
        "bash -u ./GAHereDoc.sh %s %s %s %d %d %d %d %d %d %d %d %s",...
        "full_runs",experiment_name,model_name,gen_counter+1,...
        population_size,stim_bcl,stim_num,stim_dur,stim_curr,sim_dur,...
        full_run_factor,ParamOrder_string);
    system(GAHereDoc_command)

%Make the run_type directory within the generation directory%%%%%%%%%%%%%%%
        MakeRuntypeDir_command      = sprintf(...
            "mkdir ./%s/generation_%d/%s",...
            experiment_name,gen_counter+1,"full_runs");
        system(MakeRuntypeDir_command);

%Run the MakeAdjFiles_func for FULL RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CombinedPop_array           = [...
            EliteParams_array,TrimmedChildParams_array];
        MakeAdjlFile_func(...
            experiment_name,gen_counter+1,"full_runs",...
            length(ParamOrder_array), ParamOrder_array,population_size,...
            CombinedPop_array)
        
%Save the .mat for parameter adj%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GenXAdjs_array              = [ParamOrder_array,CombinedPop_array];
        SaveMat_path                = sprintf(...
            "./%s/generation_%d/generation_%d_adjs.mat",...
            experiment_name, gen_counter+1, gen_counter+1);
        save(SaveMat_path, "GenXAdjs_array")
        
%Run the RunCarpSim.sh for FULL RUNS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulate the variants of the next generation for the full sim_dur
         RunCarpSim_command         = sprintf(...
            "bash -u ./RunCarpSim_V4.sh %d %s %s %s %d %d",gen_counter+1,...
            "full_runs",model_name,experiment_name,population_size,sim_dur);
        system(RunCarpSim_command)
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
   elseif stim_num==0 && stim_bcl==0
%% Spontaneous Runs
        % Models beat due to intrinsic automaticity. 
        % No need for quiescence check.
        % Preallocate array for mating pool population.
        % Each row is a parameter adj, each column is a variant.
        MatingPoolParams_array      = zeros(num_params, ChildrenPop_size);
        
%Tournament selection%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Two variants are chosen from the full population. Their fitness
        % scores are compared. The better fitness score variant will be
        % placed in the mating pool.
        
        for tourney_counter         = 1:ChildrenPop_size
            % Contestants is an array that contains two numbers.
            contestants             = datasample(...
                1:population_size,2,'Replace',true);
            
            % The smaller number is closer to the top of the
            % CurrentGen_struct_sorted sruct, and therefore corresponds to
            % a better fit variant
            tourney_winner          = min(contestants);
            
            % Place parameters of better scoring variant in test pool.
            MatingPoolParams_array(:,tourney_counter) =...
                CurrentGen_struct_sorted(tourney_winner).adj_data;
        end
        
%Pair Up Variants%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Randomly order the variants numbers in the test pool
        RandomOrder_array           = datasample(...
            1:MatingPool_size,MatingPool_size,"replace",false);
        
        % Pair up the variants by reshaping RandomOrder_array
        % Each column of PairedVariants_array has the 2 parent variant nums
        PairedVariants_array        = reshape(...
            RandomOrder_array,[2,MatingPool_size/2]);
        
%Crossover Operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Single or multipoint crossover. Explanation within function.
        % This is an optional operator. If skipped, the child population
        % will be a duplicate of the mating pool.
        if myApp.Crossover_checkbox.Value   == true
            ChildParams_array               = CrossOver_func(...
                MatingPoolParams_array,PairedVariants_array,num_Xpoints);
        else
            ChildParams_array               = MatingPoolParams_array;
        end

%Mutation Operator%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Randomly alter parameters to introduce more parameter diversity.
        % Explanation within function.
        ChildParams_array_mut               = mutation_func(...
            ChildParams_array);
        
%Make next generation directory%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % The mutated test child variants belong to the next generation.
        % Create new generation directory
        MakeNewGen_func(experiment_name, gen_counter+1);

%Make the run_type directory within the generation directory%%%%%%%%%%%%%%%
        MakeRuntypeDir_command              = sprintf(...
            "mkdir ./%s/generation_%d/%s",experiment_name,...
            gen_counter+1,"full_runs");
        system(MakeRuntypeDir_command);
        
%Make .par files for FULL runs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        GAHereDoc_command           = sprintf(...
            "bash -u ./GAHereDoc.sh %s %s %s %d %d %d %d %d %d %d %d %s",...
            "full_runs",experiment_name,model_name,gen_counter+1,...
            population_size,stim_bcl,stim_num,stim_dur,stim_curr,sim_dur,...
            full_run_factor,ParamOrder_string);
        system(GAHereDoc_command)
   
% Create generation_x_adjs.mat%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CombinedPop_array   = [EliteParams_array, ChildParams_array_mut];
        generation_x_adjs   = [ParamOrder_array,CombinedPop_array];
        SaveMat_path        = sprintf(...
            "./%s/generation_%d/generation_%d_adjs.mat",...
            experiment_name, gen_counter+1, gen_counter+1);
        save(SaveMat_path, "generation_x_adjs")

%Run the MakeAdjFiles_func for FULL RUNS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MakeAdjlFile_func(...
            experiment_name,gen_counter+1,"full_runs",...
            length(ParamOrder_array), ParamOrder_array,population_size,...
            CombinedPop_array)
        
% Run the RunCarpSim.sh for FULL RUNS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulate the variants of the next generation for the full sim_dur
        RunCarpSim_command = sprintf(...
            "bash -u ./RunCarpSim_V4.sh %d %s %s %s %d %d",gen_counter+1,...
            "full_runs",model_name,experiment_name,population_size,sim_dur);
        system(RunCarpSim_command)
    else
        disp(...
        "Checkpoint FAILED: If non-paced stim BCL and Num must equal 0");
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALL FUNCTIONS (Must stay at end of file)                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Make new generation directories
function GenDir_path    = MakeNewGen_func(experiment_name,gen_num)
    GenDir_path         = sprintf(...
        "./%s/generation_%d", experiment_name,gen_num);
    MakeNewGen_command  = sprintf(...
        "mkdir %s", GenDir_path);
    system(MakeNewGen_command);
    fprintf(...
        "Checkpoint passed: generation %d dir has been made", gen_num);
end

%% New adjs File Function
function MakeAdjlFile_func(...
    experiment,gen_counter,run_type,num_params,...
    ParamOrder_array,num_variants,pop)
    
    for NumParams_counter   =1:num_params
        adj_file=[num_variants;"intra"];
        adj_file_name       = sprintf(...
            "./%s/generation_%d/%s/generation_%d_%s.adj",...
            experiment,gen_counter,run_type,gen_counter,...
            ParamOrder_array(NumParams_counter));
        writematrix(adj_file,adj_file_name,'FileType',"text");
        nodes               = (0:num_variants-1)';
        adj_file            = [nodes, pop(NumParams_counter,:)'];
        dlmwrite(adj_file_name,adj_file,'-append','delimiter',' ')
    end
    fprintf("Checkpoint passed: .adj files have been created.");
end

%% AP Morphological Feature Check Function
function APFeature_array    = APFeatureCheck_func(time_data,AP_data,num_features)
% flag is a multiplier to feature values. Features with physiologically
% irrelevant values will be flagged. This causes the feature value to be
% 100 times the target, drastically worsening the fitness score.
flag                = 1000;
% Preallocate APFeature_array
APFeature_array = zeros(num_features,1);
%TROUGH SEARCH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the values and indices of most negative Vms 
[local_min,min_ind] = findpeaks(-AP_data);
% Get only trough Vms. Neccessary because local_mins might contain notches.
trough_vms          = -local_min(-local_min<0);
% Translate trough indices to actual time in ms.
trough_ind          = min_ind(-local_min<0);
trough_times        = time_data(trough_ind);
% Combine trough_times and trough_vms into one array.
trough_array        = [trough_times,trough_vms];
num_troughs         = size(trough_array,1);

% If less than 3 troughs, there is repolarization/depolarization faiure.
% Return feature check with all inf to mark as irrelevant AP.
% This will prevent errors thrown in following checks.
if num_troughs                              < 3
    for feature_counter                     = 1:num_features
        APFeature_array(feature_counter,1)  = inf;
    end
    return
end

% Extract last trough Vm as MDP. We assume last is at steady state.
MDP                 = trough_array(end,2);
% Calculate CL as time between the last two troughs
CL                  = trough_array(end,1)-trough_array((end-1),1);
% Calculate beats per minute
CL_in_sec           = CL/1000;
BPM                 = 60/CL_in_sec;

%NOTCH SEARCH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get notch Vms.
notch_vms           = local_min(-local_min>0);
% Translate notch indices to actual time in ms.
notch_ind           = min_ind(-local_min>0);
notch_times         = time_data(notch_ind);
% Combine notch_times and notch_vms into one array
notch_array         = [notch_times,notch_vms];
len_notches          = length(notch_array);


%EXTRACT LAST AP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract trace of the last full AP, which we assume is at steady state.
% We will define a full AP as the trace between two troughs.
% Find left and right indicies of troughs to isolate snapshot of last AP.
l_trough_ind        = find(time_data==trough_array((end-1),1));
r_trough_ind        = find(time_data==trough_array(end,1));
last_AP_time_data   = time_data(l_trough_ind:r_trough_ind);
last_AP_data        = AP_data(l_trough_ind:r_trough_ind);

%PEAK SEARCH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the values and indices of most positive Vms in the last full AP
[local_max,max_ind] = findpeaks(last_AP_data);
% Need to pick out true peak from post-notch peaks. True peak will
% always occur before post-notch peak, so get first val from local_max.
Peak = max(local_max(:));

% Calculate notch depth
if len_notches~=0
    notch           = Peak-notch_array(end,2);
else
    notch           = 0;
end

% Calculate AP amplitude
APA                 = Peak-MDP;

%CALCULATE MAX DV/DT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dv                  = diff(last_AP_data);
dt                  = diff(last_AP_time_data);
dvdt                = dv./dt;
% Get values and indicies of max dv/dts
[dvdt_max,dvdt_max_ind]     = findpeaks(dvdt);
% Extract max dv/dt. Findpeaks could have picked out multiple
% depolarizations depending on notch existence or abnormal morphologies.
% Find the max of the dvdt_max to get true upstroke velocity.
dvdtMax             = max(dvdt_max);
% Find time of dvdtMax
dvdtMax_ind         = dvdt_max_ind(dvdt_max==dvdtMax);
dvdtMax_time        = last_AP_time_data(dvdtMax_ind);


%CALCULATE APD20,50,90%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APDx is defined as time between dvdtMax and x% repolarization.
% V 4.0.0 only checks 3 APDs. In the future, make this user defined.
APDx                = [20,0;
                       50,0;
                       90,0];
for x               = 1:length(APDx)
    APDx_vm         = MDP+(100-APDx(x,1))*APA/100;
    APDx_ind        = find(last_AP_data>=APDx_vm,1,"last");
    APDx_time       = last_AP_time_data(APDx_ind);
    try
        APDx(x,2)   = APDx_time-dvdtMax_time;
    catch
        APDx(x,2)   = APDx_vm*flag;     
    end
    
end
APD20               = APDx(1,2);
APD50               = APDx(2,2);
APD90               = APDx(3,2);

APFeature_array     = [CL;BPM;MDP;Peak;APA;...
                       APD20;APD50;APD90;notch;dvdtMax];   
end

%% Fitness Score Calculation Function
function PercentError_array = FitnessScore_func(...
    APfeature_data,target_features)
% The fitness function is the sum(percent error for each feature).
% APfeature_data should be a cell containing arrays. Each array is the
% feature value set for a variant.

% Preallocate array for percent errors
   PercentError_array = zeros(length(target_features),length(APfeature_data));
% Each column of PercentError_array is for a variant
% Each row  is for a feature
   for variant_counter  = 1:length(APfeature_data)
       PercentError_array(:,variant_counter) =...
       (abs((APfeature_data{variant_counter}-target_features))...
       ./target_features)*100;
   end
   
end

%% Crossover Function
function ChildParams_array = CrossOver_func(...
    Pool,PairedVariants_array,num_Xpoints)

    % Preallocate ChildParams_array
        ChildParams_array  = zeros(...
            size(Pool,1),numel(PairedVariants_array));
        
    for pair_counter        = 1:size(PairedVariants_array,2)
        % Create arrays containing param adjs of paired parent variants
        parent1_params      = Pool(:,PairedVariants_array(1,pair_counter));
        parent2_params      = Pool(:,PairedVariants_array(2,pair_counter));
        
        % Randomly generate an array of crossover points
        Xpoints_array       = datasample(...
            1:length(parent1_params),num_Xpoints,"replace",false);
       
        % Children params begin as identical to respective parents 
        child1_params       = parent1_params;
        child2_params       = parent2_params;
        
        % Create placeholder param arrays for holding swapped values
        temp1_params        = parent1_params;
        temp2_params        = parent2_params;
        
        
        % Children params recombine to generate new children "genotype"
        % At every Xpoint, every param after that point swaps between the
        % two children.
        for Xpoint                      = Xpoints_array
            temp1_params(Xpoint:end)    = child2_params(Xpoint:end);
            temp2_params(Xpoint:end)    = child1_params(Xpoint:end);
            child1_params               = temp1_params;
            child2_params               = temp2_params;
        end
        
        % Insert resulting children params in the consolidated array
        ChildParams_array(:,2*pair_counter-1)    = child1_params;
        ChildParams_array(:,2*pair_counter)      = child2_params;
    end
end
    
%% Mutation Function
function MutChildParams_array = mutation_func(ChildParams_array)

    num_variants            = size(ChildParams_array, 2);
    num_params              = size(ChildParams_array, 1);

    % Preallocate MutChildParams_array
    MutChildParams_array = zeros (num_params,num_variants);
        
        % Each child variant will undergo mutation
        for variant_counter     = 1:num_variants
            % Each parameter will have a 50% chance of getting mutated
            Params2Mut_array    = datasample([0,1],num_params)';
            % The amount the chosen parameters will be mutated by is chosen
            % from a normal distribution curve: mu = 0 sigma = 0.25
            Amt2Mut = normrnd(0,0.25);
            % All of the chose pameters increase/decrease by this amount
            % The absolute value is to preven a negative adj factor
            MutChildParams_array(:,variant_counter) =...
                abs(ChildParams_array(:,variant_counter)...
                + Params2Mut_array .* Amt2Mut);
        end
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Visualization Functions                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotFeatures_func(...
    APfeature_data,target_features,target_SE)
    % APfeature_data should be a cell containing arrays. Each array is the
    % feature value set for a variant.
    % Combine all the variant arrays into one big array. Preallocate:
    num_features = length(target_features);
    num_variants = length(APfeature_data);
    APfeature_array = zeros(num_features,num_variants);
    for variant_counter = 1:num_variants
        APfeature_array(:,variant_counter) =...
            APfeature_data{variant_counter};
    end
    
    % Initialize a tiled plot layout
    tile_dim = ceil(sqrt(num_features));
    tiledlayout(tile_dim,tile_dim);
    
    % Target features:
    feature_labels = [...
        "Cycle Length (ms)";
        "Beats/Minute";
        "Diastolic Potential (mV)";
        "Peak (mV)";
        "AP Amplitude (mV)";
        "APD20 (ms)";
        "APD50 (ms)";
        "APD90 (ms)";
        "Notch Depth (mV)";
        "vmaxUp (mV/ms)";
        ];
    
    % In each tile, plot the scatter of feature value distribution
    for tile_counter = 1:num_features
        nexttile
        
        % Upper tolerance
        top_line = target_features(tile_counter)+target_SE(tile_counter);
        hold on
        % Lower tolerance
        bot_line = target_features(tile_counter)-target_SE(tile_counter);
        hold on
        % x-coordinates is the variant number
        x = 1:num_variants;
        % y-coordinate is the value for a feature
        y = APfeature_array(tile_counter,:);
        % Define graph limits
        boundaries = [0, num_variants];
        
        % Shade the target feature range
        area(boundaries,[top_line,top_line],bot_line);
        area_color = [194, 224, 164]./255;
        colororder(area_color)
        hold on
        
        % Plot the feature values for each variant as scatter
        % scatter_color = [191, 106, 8]./255;
        scatter(x,y, 'filled', 'SizeData', 30, 'MarkerFaceAlpha', 0)
       
        % Aesthetics 
        text(x,y,split(num2str(1:num_variants)),'FontSize',6)
        xlim(boundaries);title(feature_labels(tile_counter))
        set(gca,'XColor', 'none')
        set(gcf,'position',[100,100,1000,1000])
    end 
    disp("Checkpoint passed: figure has been generated and saved")
end