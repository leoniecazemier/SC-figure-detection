%% sorted_units
% 
%  'ses' = Mouse and date, links to session_info
%  'depthR' = unit depth relative to SC surface
%  'trialStop' = last trial before mouse stopped performing task (trials after this were excluded)
%  'inclStableCell' = inclusion criterion, indicates whether cell was stable during recording sesssion
%  'inclRFavailable' = inclusion criterion, indicates whether RF was detectable
%  'RFPOS' = RF position 1 = Inside the figure, 2 = Outside the figure, 3 = On the edge, 
%       4 = Unclear, 5 = not analysed, unit already excluded 
%  'perfOrientation' = task performance on orientation task
%  'perfPhase'  =  task performance on phase task
%  'inclOrientation' = inclusion criterion, indicates whether performance on orientation task was sufficient during this session
%  'inclPhase' = inclusion criterion, indicates whether performance on phase task was sufficient during this session
%  'spikeRaster' = Spikes during each trial
%  'inclResponsiveCell' = inclusion criterion, indicates whether unit was visually responsive

%% session_info
 
%  'ses' = Mouse and date, links to sorted_units
%  'recordingSide' = hemisphere that was recorded from
%  'trialGroups' = trial behavioural result (see below for translation)
%  'itidurations' = Inter-trial-interval durations
%  'FigureOrientation' = Figure grating orientations, either 0 or 90 degrees
%  'RTs' = Reaction times

% Trial group translation:
% 1 = Orientation, Fig on RF, Hit
% 2 = Orientation, Fig on RF, Error
% 3 = Orientation, Fig on RF, Miss
%
% 4 = Orientation, Ground on RF, Hit
% 5 = Orientation, Ground on RF, Error
% 6 = Orientation, Ground on RF, Miss
%
% 7 = Phase, Fig on RF, Hit
% 8 = Phase, Fig on RF, Error
% 9 = Phase, Fig on RF, Miss
%
% 10 = Phase, Ground on RF, Hit
% 11 = Phase, Ground on RF, Error
% 12 = Phase, Ground on RF, Miss
%
% 13 = Contrast, Fig on RF, Hit
% 14 = Contrast, Fig on RF, Error
% 15 = Contrast, Fig on RF, Miss
%
% 16 = Contrast, Ground on RF, Hit
% 17 = Contrast, Ground on RF, Error
% 18 = Contrast, Ground on RF, Miss

%% 
