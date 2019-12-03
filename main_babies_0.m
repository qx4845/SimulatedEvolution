%% main_babies_0.m *******************************************************
% This is the primary m-file to run the evolution simulation, babies. For reference,
% the original defaults are: overpop = 0.25; death_max = 0.70; NGEN = 2000; range = 7; 
% IPOP = 300; landscape_movement = 2; landscape_heights = [1 4]; basic_map_size = [12 12];
% Cluster information must be run post simulation for RAM resource reasons on most
% computers. However, enough information may be output so that every indiv may be traced
% by generation, cluster, phylogeny, death, etc. Furthermore, landscape scenarios to model 
% Neutral Theory (#_Flatscape), static Natural Selection (Frozenscape), rolling Natural 
% Selection (#_Shifting), and drastic Natural Selection (#_Shock) are included. There is
% no feedback mechanism as an option for now. Reproduction options include Assortative
% Mating, Bacterial Splitting, and Random Mating. Assortative and Random are sexual
% reproduction types, whereas Bacterial is asexual. Finally, one may choose simulations of
% competition of many mutabilities in the population, competition of two mutability 
% populations, or single mutability values for the population.

%% START CLEAN & RANDOM
clear; clc; close all;
rng('shuffle','twister');

%% PARAMETERS
% Data creation options
  % Save options
  save_coords = 0; %save trace_x & trace_y? no = 0, yes = 1
  save_seeds = 0;
  save_parents = 0; %save parents? no = 0, yes = 1
  save_kills = 1; %save kills? no = 0, yes = 1
    save_rivalries = 1; %save rivalries? no = 0, yes = 1
  source = 'E:\Wedekind\babies_root\';
  for_external = 2; %output to Data if 0, to Data\external if 1, to both Data & Data\external if 2
  split = 0; %output to Data if 0, output to Data\#_#+split if nonzero
  make_dir = 1; %create new directory for sims if it doesn't exist? no = 0, yes = 1
  do_cd = 0;  %change directories to wherever data will be saved? no = 0, yes = 1
  OS = 0; %windows = 0, linux = 1
  write_over = 0;
  pool_size = 0;
  record = 0; % not sure what the variable will do yet. need to find out
  
  % What data to create
  no_bio = 0; %model biology (1 allows pops to go to < limit and continue running)
  only_lt = 0; %record only times to fixation if 1 and disregard all other data
  only_relax = 1;
  do_simulations = 1; %generates population, trace_x, trace_y, trace_noise, trace_cluster_seed*
  do_clustering = 0;  %generates trace_cluster_seed* & allows building and locating clusters
    do_cluster_seeds = 1;
    do_build_clusters = 1;  %generates num_clusters, trace_cluster, orgsnclusters, 
    do_locate_clusters = 1; %generates centroid_x, centroid_y, cluster_diversity
    do_build_gyration_radii = 1;
    do_build_diameters = 0;
    do_correlation_lengths = 1;
    fix = 0;
  do_genealogies = 0; %generates lineage info of original population & of the clusters
    do_indiv_lineage = 1; %generates num_descendants
    do_indiv_cluster_lineage = 1; %generates num_descendant_clusters, descendant_clusters
    do_cluster_lineage = 1; %generates num_clusters_produced, clusters produced, 
                            %num_clusters_fused, clusters_fused
                            

% Simulation numbers
SIMS = [1:10]; %simulation identifiers used in simulation looping

% Initial parameter settings
NGEN = 2000; %number of generations to run
IPOP = 8196; %number of initial population
  do_default_density = 0;
limit = 3; %minimum cluster size/extinction population +1
loaded = 0; %choose to load predefined variables, 0 = no load, 1 = load load_name variables
load_name = ['']; %string which identifies predefined variables babies and basic_map
append = 0;

% Birth settings
  % Reproduction options
  distribution = 0; %offspring distrubtion: 0=uniform, 1=normal
  reproduction = 1; %assortative mating = 0, bacterial cleaving = 1, random mating = 2
  % Mutability option
  exp_type = 0; %same mutability = 0; competition = 1; duel = 2
  % Single Mutability settings
  dmu = 0.0;  mutability = [0.3];%[0.38:dmu:0.52];%[0.26:dmu:0.45];%0.350:dmu:0.400];%:dmu:0.420];
  % Two mu competition settings
  bi = [1.31 2.31; 150 150]; %[mu1 mu2; IPOP1 IPOP2];
  % Mutability competition settings
  range = 7; % 0 to range possible mutabilities

% Death settings
  % Local options
  dop = 0.25; overpop = 0.25; %if closer than this distance, overpopulated, and baby dies
  random_walk = 0; %0 = coalescing, 1 = annihilating
  % Global options
  ddm = 0.01;  death_max = [0.21]; %percent of random babies dying varies from 0 to this value.
  indiv_death = 0;  %random percentage of entire pop dies = 0; individual probability = 1

% Landscape options
shock = 0; %Set this to 1 to generate a new random map every landscape_movement generations
  shock_heights = [1 1]; %min and max landscape heights during shock
  shock_duration = 2; %duration of a shock (generations)
  shock_repeat = 0; %repeated shocks, 0 = one shock, 1 = shocks every landscape_movement generations
  shock_over = 0; %flag for when single shock has occurred, 0 = hasn't finished, 1 = finished shock
landscape_movement = 1; %land moves every "landscape_movement" generations
  %There is only a default shift of 1 basic_map row, we could add this in if 
  %we really want to, but it may be unnecessary.
landscape_heights = [1 4]; %min and max of landscape; for flat landscapes, only min is taken
  %If both values are the same, then a flatscape will be generated.
basic_map_sizes = [6];%[[12]' [12]']; %X and Y lengths for the basic map size
linear = 0;
  %May just put in a check on whether one basic_map_size values is 1
periodic = [0 0]; %periodic boundary conditions for the x & y coords, respectively
  %if [1 0] or [0 1], then cylindrical boundaries wrapping x or y edges, respectively
  %if [1 1], then the landscape becomes toroidal
  %does not yet work with varied landscape iterpolations
  
% Initialize any global variables
if only_lt, global lifetimes; global populations; end %#<NUSED,TLEV>

global SIMOPTS;
%% Simulations loop
%if matlabpool('size')==0 && pool_size>0, matlabpool open;  end
for bms = basic_map_sizes,  
  basic_map_size = [bms bms]; 
  if do_default_density,  %IPOP above assumed as if for 12x12
    IPOP = scale_IPOP(IPOP,basic_map_size); %IPOP = IPOP/37544*(max_population(bms))
  end
for op = overpop, 
for dm = death_max, 
for mu = mutability, 
  % Bundle simulation options into SIMOPTS
  SIMOPTS = struct('record',record,'save_seeds',save_seeds,'save_parents',save_parents,...
    'save_kills',save_kills,'save_rivalries',save_rivalries,...
    'source',source,'OS',OS,'pool_size',pool_size,...
    'no_bio',no_bio,'only_lt',only_lt,'append',append,...
    'write_over',write_over,'split',split,...
    'make_dir',make_dir,'do_cd',do_cd,'for_external',for_external,...
    'do_simulations',do_simulations,...
    'do_clustering',do_clustering,'do_cluster_seeds',do_cluster_seeds,...
    'do_build_clusters',do_build_clusters,...
    'do_locate_clusters',do_locate_clusters,'do_build_gyration_radii',do_build_gyration_radii,...
    'do_build_diameters',do_build_diameters,'do_correlation_lengths',do_correlation_lengths,...
    'fix',fix,'only_relax',only_relax,...
    'do_genealogies',do_genealogies,'do_indiv_lineage',do_indiv_lineage,...
    'do_cluster_lineage',do_cluster_lineage,'do_indiv_cluster_lineage',do_indiv_cluster_lineage,...
    'SIMS',SIMS,'NGEN',NGEN,'IPOP',IPOP,'limit',limit,'loaded',loaded,'load_name',load_name,...
    'op',op,'random_walk',random_walk,'dm',dm,'indiv_death',indiv_death,...
    'shock',shock,'shock_heights',shock_heights,'shock_duration',shock_duration,...
    'shock_repeat',shock_repeat,'shock_over',shock_over,...
    'landscape_movement',landscape_movement,'landscape_heights',landscape_heights,...
    'basic_map_size',basic_map_size,'linear',linear,'periodic',periodic,...
    'distribution',distribution,'reproduction',reproduction,'exp_type',exp_type,'mu',mu,...
    'bi',bi,'range',range,'run',[]);
  if do_simulations==1, sim_times = Simulations(); end  
  if do_clustering==1,  clu_times = Clustering();  end  
  if do_genealogies==1, lin_times = Genealogies(); end  
end
end
end
end
% if matlabpool('size')>0, 
%   matlabpool close;  
% %   save('CPAR_time','asdf');
% else, 
% %   save('REG_time','asdf');
% end
%% Update information
% updates as of Apr 2013 ADS & DMK
% 1) Fixed a bug in Overpopulation Death around February 2013 which caused a mixup in the organisms
% that survive overpopulation/competitive death. Any data before February/March 2013 cannot be
% used reliably for individual genealogical information. -ADS
% 2) Periodic boundaries have been added for cylinders with absorbing end cap boundaries & toroid
% topology. AdjustLandscape.m has not yet been updated, nor have many analysis functions which rely
% on distances and spatial measures. At this point, take care with periodic data until further
% notice. -ADS
% 3) Can now measure cluster diameters by using a Floyd-Warshall algorithm in
% build_cluster_diameters.m. This function was added to the Clustering.m function which runs with
% simulations instead of a previous version which ran on the analysis side
% (get_cluster_diameters.m).
% 4) Parallelization has been added to some simulation and clustering algorithms. Use of CPU, GPU,
% and the campus cluster Bortas are available. At this point Bortas is running a trial version of
% Matlab's Distributed Computing Toolbox, so cluster use may not be available past April. Some new
% functions and scripts include: CPUPAR_get_cluster_diameters.m (now obsolete),
% CPUPAR_build_cluster_diameters.m, GPUPAR_build_cluster_diameters.m. There is some work on
% converting Simulations to work with Bortas by putting runs in parallel, but this has not been
% debugged sufficiently yet. Since Bortas is a unix based system, there is now an operating system
% option "OS" to work with the different file system. The main changes are within NameAndCD.m, but
% the details have not been fully sorted out. -ADS
% 5) Added get_path_length.m which measures the "time" steps of cluster formation. There is a
% complementary script which saves results to an excel file called write_path.m. -DMK
% 6) Can run purely dynamical simulations with "no_bio", which means simulations can run
% until every single organism has died out rather than stopping for populations dropping below
% "limit". -ADS
% 7) Can run relaxation simulations which will end when the population reaches "relaxed". -ADS
% 8) Population data is now ALWAYS saved even for lifetimes & relaxations. This can be dangerous for
% RAM resources since long time simulations O(NGEN)>~10^6 coupled with many simulations
% O(SIMS)>~10^2 may use up too much memory. -ADS
% 9) There is a new folder set of functions labeled "helper". Functions within are support
% functions and may, with a little alteration, be used in more general codes. Functions included:
% cat_row.m, check_variable.m, indiv_update.m, make_data_name.m, mat_exist.m, NameAndCD.m,
% par_script_gen_update.m, print_check_result.m, proper_name.m, script_gen_update.m,
% try_catch_load.m. -ADS
% 10) There is a new analysis script which checks aspects of data files to make sure that they are
% compatible and exist. The main script for this is check_data_set.m, and it uses data_check_A.m,
% data_check_B.m, data_check_C.m, open_all_data.m, get_all_data.m, initialize_check_strings.m,
% initialize_check_values.m, data_check_initial.m, data_check_time.m, data_check_dimA.m,
% data_check_dimB.m, data_check_length.m, data_check_min.m, data_check_max.m, and print_report.m.
% This script is only partly done, so more data_check_*.m functions will almost certainly be created
% soon. The reports determine whether the given data set measures pass or fail basic checks which
% are determined by the bounds of how simulations are allowed to run and the measured values are
% listed with the corresponding expected values. Using this on a data set will be recommended since
% it will provide a sturdy level of confidence for data about to be published or presented. -ADS
% 11) There are updates to ALL simulation, clustering, and genealogies functions. Mostly, the
% updates reflect changes corresponding to simplified code, the addition of global populations,
% addition of newer/updated helper functions, more parallelized code, etc. -ADS
% 12) There are updates to several analysis scripts including: write_lifetimes.m, write_clusters.m,
% write_populations.m, write_R.m, get_characteristic_length_distribution.m,
% get_correlation_length_distribution.m, generalize_base_name.m, get_coalescent_probabilties.m,
% get_indiv_lineage_broken.m, get_lineage_stats.m get_cluster_lineage_hists.m, get_R.m, get_kills.m,
% kill_v_pop.m, get_num_clusters.m, write_path.m, get_cluster_tree.m, get_cluster_diameters.m,
% CPUPAR_get_cluster_diameters.m, check_seed_distances.m, get_finite_cluster_of_gen.m,
% get_finite_clusters.m, get_cluster_dimension.m,
% get_fractal_cluster_dimension_from_number_density.m, get_path_length.m, get_populations.m,
% lifetimes_relaxations_distributions.m, get_correlation_length_growth_rate.m,
% get_lifetimes_to_survival_probabilities.m, plot_indivs.m, data_check_.m data_check_initial.m,
% data_check_time.m, data_check_dimA.m, data_check_dimB.m, get_fisher_exponent.m,
% get_transition_probabilities.m, get_cluster_density_dimension.m, get_all_data.m, data_check_B.m,
% data_check_A.m, data_check_length.m, data_check_min.m, data_check_max.m,
% initialize_check_strings.m, print_report.m, data_check_C.m, open_all_data.m, check_data_set.m,
% data_check_notes.m, plot_lineage.m, get_spatial_correlation_lengths.m,
% get_characteristic_length_distributions.m. Finally, main_analysis.m has been updated to reflect
% the changes to all of the mentioned updated (some new) scripts and functions above. -ADS
%
% updates as of Feb 2012 ADS & DMK
% 1) Capability to change random death from population size based to 
% individually based was added. -DMK
% 2) There is a new analysis function which can look at percolation of
% largest clusters. The new script is percolation_lengths.m. -ADS
% 3) Writing Excel files for use in Sigma Plot has been modified and
% expanded. The script, write_sigma_plot_files.m, still exists, but it
% includes new sub-scripts: write_populations.m, write_clusters.m,
% write_time_to_fixation.m. This rearrangement results in three separate
% xls files for each type of data. New capability to write out a similar
% file for the nearest neighbor measure, R, was also added. This
% sub-script, write_R.m, may be enacted by toggling record_R in
% main_analysis.m. -ADS
% 4) Added mkdir capability to NameAndCD.m. Any new simulations do not
% require manually making the appropriate directories. Instead, NameAndCD
% will automatically create new directories to house a new set of
% simulation data. -ADS
%
% updates as of Jan 2012 ADS -
% 1) There are many new analysis functions added. See main_analysis.m for
% all the fancy new analyses.
% 2) There is basic functionality for loaded runs. Use main_loader to make load
% files from pre-existing data. There is no functionality yet for self made
% data starts. It may be while until that's included, since it's not yet
% needed.
% 3) Included directory handling to NamingScheme.m which is now
% NameAndCD.m. See NameAndCD.m documentation for how to customize the use
% of the source parameter in main_babies_#.m for the location of your data.
% 4) Took out the flat option. For flatscapes, just set both
% landscape_heights values to the same value.
%
% updates as of July 2011 ADS - 
% 1) An option to load predefined variables for babies and basic_map is now available. Note 
% that predefined shiftedscapes are not yet possible with this update. Simulations will run
% akin to a historical contingency idea based on a moment of landscape in time and not an 
% underlying landscape throughout time. Functions updated include: main_babies, Simulations, 
% NamingScheme.  
% 2) The function, Record (now on 4th iteration), no longer saves a trace_noise file if the 
% exp_type is 0 (for single mutability). Functions updated include: Simulations, Record.
% 3) Fixed a bug with populations less than limit being saved, especially in bacterial
% simulations. The function updated to remove this bug is Simulations. Two new functions
% were created to fix population data that includes the excess population, modify_data.m
% and fix_population.m. See their help info for more information. The final updates for
% this bug include modification of build_clusters, locate_clusters, and
% build_trace_cluster_seed which had been temporarily changed to account for the excess
% population data value. They are now back to running as normal with the exception that
% build_trace_cluster_seed and build_clusters now take limit as an input.
%
% updates as of June 2011 ADS - 
% 1) Reproduction options for assortative mating (default), random mating and bacterial
% splitting included. Functions updated for reproduction toggling include: main_babies, 
% Simulations, Generations, MakeBabies, FindMates, NamingScheme.
% 2) Re-added in option to save kill counts. Functions updated include: main_babies,
% Simulations, Generations, OverpopulationLimit, RandomDeath, CliffJumpers, Record.
% 3) New option to determine minimum population needed to run (Limit). Since bacterial
% species is very controversial, the mimumum population needed to determine a species may
% be different sizes, so change Limit to accomodate for different potential minimum
% cluster sizes. Functions updated include: main_babies, Simulations, Generations.
% 4) Some function names end with a 2 or 3 to indicate these latest changes.
% 5) Changed the two competiting mutability handling so that exp_type = 2 determines
% competition between two mutabilities. Functions updated include: main_babies,
% Simulations, setInitialMutabilities, NamingSchemes.
%
% updates as of May 2011 ADS -
% 1) Created two new primary functions, Generations and Simulations.  These
% are meant to help break up babies so that multiple experiments may be run
% consecutively.  Generations is just the generations loop of babies_noBias
% and Simulations is just the simulations loop of babies_noBias.
%
% updates as of March 2011 ADS - 
% 1) Revamped functions throughout the program into explicit function calls. 
% 2) Included naming scheme according to competition or not, distribution 
%   of offspring type, and landscape scenarios.
% 3) Can now save parents for absolute phylogeny determination.  Some 
%   function names end with a 2.  This indicates that they are set up to 
%   save the parents.  Without the 2, the functions will not save data for 
%   parents.
% 4) Extended range of possible mu values within competition up from 1 to 7. 
% 5) Eased changes needed between scenarios so only the beginning values 
%   need changing. 
% 6) Altered the order of Overpopulation Limit check so that there is no 
%   longer bias for the early organisms listed in variable babies. 
% 7) Some simulations require more resources than available on most
%   lab computers, so I've set this up so that only 1st and 2nd nearest
%   neighbors are saved.  This allows for cluster data to be generated from
%   the trace_cluster_seed file.  Changes associated with this are 
%   designated on lines with a long string of the percent signs.  The 
%   original cluster variables are still in place, but are commented out.  
%   Use the script, determine_cluster_info, to determine full cluster data 
%   set.
% 8) Defaults are noted under the Parameters section.  If running defaults,
%   there will be no indication of those settings within the
%   simulation/filename.  However, the naming scheme does not yet account
%   for every parameter change you may make, so that will need to be
%   included in order to make the filenames describe the simulation
%   scenario explicitly.
%
%
% evol model ND & SB 2009 - Babies 30  --   Shifting landscape with tracers
% new version summer 2009 to have babies with different noises competing
% against each other
%
% Overpopulation derived from small circles surrounding each individual.
% Mating using the matching hypothesis - mate with the one closest to you
% (in terms of phenotype = assortative mating).
%
% A fitness landscape lies below the phenotype placements.  It is 
% based on a random matrix, with values in the landscape (rounded)
% equaling the number of babies an organism with that fitness would produce.
%
% SO - a certain phenotype puts you in a place on the map, where a certain
% fitness lies below you, and determines how many babies you produce.
%
% A Clustering mechanism is included, with statistics recorded.
%
% Still needing a feedback mechanism for landscape.. do have a changing
% landscape option based on randomness.
%**************************************************************************

%% old stuff

%cluster_perc_cutoff = 4.293; %when grouping organisms, percent cutoff represents 
%the percent change in the distance between the organisms when the next baby 
%is included.  A way to limit clustering based on closeness -if too many
%are clustered together, go lower. if clusters are too sparse, go higher...
%NOTE:  I've never used cluster_perc_cutoff. This was included in Nate's
%original design of the program.  I think we've abandoned this concept for
%various reasons.