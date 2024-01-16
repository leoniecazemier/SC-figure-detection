%% Data descriptions

% The struct 'behaviour' has a separate field containing data for each
% mouse. For each mouse, properties of trials are concatenated across
% sessions. There is also some summarized information. 

% The following properties are present: 

% optolats: Latency of optogenetic interference relative to stimulus onset (NaN if no interference)
% oris: Orientation of figure (0 or 90 deg)
% RTs: Reaction times in ms
% contrasts: Figure grating contrast
% bgcontrasts: Background grating constrast
% gavepas: Whether a passive reward (not based on behaviour) was given during trial
% drum: whether drumming was active (drumming means immediate repetition of error trialuntil correct response is given)
% semidrum: whether semidrum was active (semidrum means repetition of error trial at the end of a sequence of 4 trials)
% OOP: whether trial is phase task trial
% repeated: Whether trial is a repetition due to (semi) drumming
% bglum: background grating luminance
% OptoOn: Whether optogenetic interference trials were interleaved at this
        % time of the session (interleaved opto was only turned on when mice performed well)
% reactions: Whether trial was a hit, error, or missed trial by the mouse
% sides: On which side the figure was displayed
% useBase: Whether the session this trial belonged to had good enough performance on contrast task to include in analysis
% useOri: Whether the session this trial belonged to had good enough performance on orientation task to include in analysis
% useOOP: Whether the session this trial belonged to had good enough performance on phase task to include in analysis
% propR: Proportion of right-side trials (usually 50%, sometimes if mice
        % had severe bias we would tailor the proportions to try and get rid of bias
% perfBase: Performance of mouse on contrast task during the session this trial belonged to      
% perfOri: Performance of mouse on orientation task during the session this trial belonged to   
% perfOOP: Performance of mouse on phase task during the session this trial belonged to  
% RTrightVec: Lick times on right spout
% RTleftVec: Lick times on left spout

% Summary statistics: 

% ---- Contrast task ----
% baseOptoH: For each latency, number of hits during optogenetic interference
% baseOptoE: For each latency, number of errors during optogenetic interference
% baseOptoP: For each latency, accuracy of task performance
% baseOptoLats: optogenetic interference latencies for the baseOpto summarized statistics
% basePperf: For each latency, p-value for test whether performance during optogenetic interference
        %  is lower than without optogenetic interference
% basePchance: For each latency, p-value for test whether performance during optogenetic interference
        %  is above chance level
 
% ---- Orientation task -----
% As base task, all variables named full...Ori

% ---- Phase task -----
% As base task, all variables named full...OOP
      



        
        


