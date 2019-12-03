%% Simulations.m
% Inputs - none
% Outputs - lifetimes
% [lifetimes] = Simulations()
function [lifetimes] = Simulations()
global SIMOPTS;

%% Initialize generation variables
GENVARS = struct('babies',[],'basic_map',[],'land',[],'run',[],'exp_name',[]);

%% Initialize output
if SIMOPTS.only_lt,  
  global lifetimes; 
  global populations; 
end %#ok<*TLEV>h
lifetimes = zeros(length(SIMOPTS.SIMS),1); %number of generations each simulation lasts

%% Naming scheme
[base_name,dir_name] = NameAndCD(SIMOPTS.make_dir,SIMOPTS.do_cd); 

%% Determine whether to do lifetimes or relaxations


if   ~SIMOPTS.only_lt || SIMOPTS.write_over,  
 
%% Simulation loop
i = 0;
for run = SIMOPTS.SIMS 
  go_pop = 1; go_tx = 1;  go_ty = 1;
  
  new_dir_name = split_cd(dir_name,SIMOPTS.split,run,SIMOPTS.make_dir,SIMOPTS.do_cd);
  
  if ~SIMOPTS.write_over, 
    if ~SIMOPTS.only_lt, 
      fprintf(['Attempting ' base_name int2str(run) '\n']);
    elseif SIMOPTS.only_lt, 
      fprintf(['Attempting lifetimes of ' base_name '\n']);
    end
  end
  %% Determine whether to do simulation for the run value
  if ~SIMOPTS.only_lt, 
    go_pop = mat_exist([new_dir_name 'population_' base_name int2str(run)]);
    go_tx = mat_exist([new_dir_name 'trace_x_' base_name int2str(run)]);
    go_ty = mat_exist([new_dir_name 'trace_y_' base_name int2str(run)]);
  end %if not only_lt and only_relax
  if (~go_pop || ~go_tx || ~go_ty) || SIMOPTS.write_over || SIMOPTS.only_lt,  
  run_name = int2str(run);
  GENVARS.run = run;
  i = i +1;

  if isempty(SIMOPTS.load_name) && SIMOPTS.loaded, SIMOPTS.load_name = base_name(1:(end-1)); end

  %% FITNESS LANDSCAPE
  if ~SIMOPTS.loaded, 
    if SIMOPTS.landscape_heights(1)~=SIMOPTS.landscape_heights(2), 
      GENVARS.basic_map = rand(SIMOPTS.basic_map_size)*SIMOPTS.landscape_heights(2);
    elseif SIMOPTS.landscape_heights(1)==SIMOPTS.landscape_heights(2), 
      GENVARS.basic_map = SIMOPTS.landscape_heights(1)*ones(SIMOPTS.basic_map_size);
    end
  else, 
    if SIMOPTS.landscape_heights(1)==SIMOPTS.landscape_heights(2), 
      GENVARS.basic_map = SIMOPTS.landscape_heights(1)*ones(SIMOPTS.basic_map_size);
    else, 
      load([new_dir_name 'basic_map_' SIMOPTS.load_name]);
      GENVARS.basic_map = basic_map;  clear basic_map
    end
  end
  % for linear landscapes
  if SIMOPTS.linear,  
    GENVARS.land = ceil(interp1(GENVARS.basic_map,2));
  else, 
    GENVARS.land = ceil(interp2(GENVARS.basic_map,2)); %round up the landscape values
  end

  %% ORIGINAL POPULATION
  if ~SIMOPTS.loaded, 
    %Set initial mutability(ies) for population
    [IBN] = setInitialMutabilities();
    %Set initial locations for population
    location = [rand(SIMOPTS.IPOP,1)*(size(GENVARS.land,1)-1)+1 ...
      rand(SIMOPTS.IPOP,1)*(size(GENVARS.land,1)-1)+1]; 
    %Build the first population
    GENVARS.babies = [location IBN];  clear IBN location
  else, 
    %Load the first population
    load([new_dir_name 'babies_' SIMOPTS.load_name]);
    GENVARS.babies = babies;  clear babies
  end

  %just to let you know in the Command Window what you're running
  if ~SIMOPTS.only_lt, 
    fprintf(['normal sim of ' base_name run_name '\n']); 
    GENVARS.exp_name = [base_name run_name];
  elseif SIMOPTS.only_lt, 
    simstr = length(int2str(SIMOPTS.SIMS(1))) +1 ...
      +length(int2str(SIMOPTS.SIMS(end)));
    fprintf(['lifetimes of ' base_name(1:(end-simstr)) run_name '\n']);
    GENVARS.exp_name = [base_name];
  end

  %% generation loop
  [GENOUTS] = Generations(GENVARS);

  %% RECORD INFORMATION
  if ~SIMOPTS.only_lt, 
    RecordData(base_name,GENVARS,GENOUTS,dir_name);  % save data (returns a 1 for saved if successful)
  else, 
    populations(i,1:numel(GENOUTS.population)) = GENOUTS.population;
  end

  lifetimes(i) = GENOUTS.finished;

  %% RESET & CLEAN-UP
  clear GENVARS GENOUTS;
  end %population, trace_x, or trace_y does not exist OR write_over is true
end %SIMS
end %only_lt & only_relax & write_over
%% Record lifetimes
if SIMOPTS.only_lt, 
  save_lifetimes(base_name,dir_name);
  save_populations(base_name,dir_name);
end
end