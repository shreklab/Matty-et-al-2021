% This code analyzes tracking videos for 
% Figure 2 of Matty et al (2021)
%
% Created by Jessica Haley
% 8/9/2021
% MATLAB 2021a

%% Prepare Workspace (Clear and Add Paths)

clear
close all

addpath(genpath('Z:\jhaley\analysis'))
% addpath(genpath('Z:\jhaley\experiments\behavior'))
% addpath(genpath('F:\EuniceData'))

dataPath = 'Z:\jhaley\analysis\Eunice\';
figurePath = 'Z:\jhaley\analysis\Eunice\Figures\';
analysisPath = 'Z:\jhaley\analysis\Eunice\Analysis\';

%% Get File Names for Analysis

hungerFolders = {'WF';'FD'}; % Well-Fed, Food Deprived
odorFolders = {'Full';'NoCu';'NoD'}; % Copper + Diacetyl, Diacetyl Only, Copper Only

videoFile = {};
wormLabFile = {};
hunger = {};
odor = {};

for i = 1:length(hungerFolders)
    for j = 1:length(odorFolders)
        % Get file info for the current directory and extract file names
        theseVideos = dir([dataPath,hungerFolders{i},'\',...
            odorFolders{j},'\*.avi']);
        videoFile = [videoFile;{theseVideos.name}'];
        
        % Get file info for the current directory and extract file names
        theseWormLab = dir([dataPath,hungerFolders{i},'\',...
            odorFolders{j},'\Position\*.csv']);
        wormLabFile = [wormLabFile;{theseWormLab.name}'];
        
        % Store satiety and odor conditions
        hunger = [hunger;repmat({hungerFolders{i}},length(theseVideos),1)];
        odor = [odor;repmat({odorFolders{j}},length(theseVideos),1)];
    end
end

%% Define Experiment Info

numExp = length(videoFile); % # experiments
frameRate = 3; % frame rate
wormLabScale = ones(numExp,1); wormLabScale(1:9) = 32.37; % arbitrary scale used in wormlab
xPix = 1024; yPix = 1024;

% 2 cameras were used in these experiments. Files are labeled #f or #s to 
% denote which setup was used for each video. Scale bars will be calculated
% seperately for each camera.
camera = zeros(numExp,1); % 0 = s, 1 = f
for i = 1:numExp
    camera(i) = isempty(strfind(videoFile{i},'f_0'));
end

%% Get Video Info, Mask Sharpie/Shadow Marks, and Get Scale

% Start parallel processing
if isempty(gcp('nocreate'))
    pool = parpool(24); % start parellel pool
end

numFrames = zeros(numExp,1);
firstFrames = zeros(numExp,xPix,yPix);
firstFramesBW = zeros(numExp,xPix,yPix);
videoMasks = zeros(numExp,xPix,yPix); % all sharpie/shadow marks
lineMasks = zeros(numExp,xPix,yPix); % line where copper was dripped
arcMasks = zeros(numExp,xPix,yPix); % edge(s) of plate to be used as scale reference
lineOrients = zeros(numExp,1);
lineExtrema = zeros(numExp,8,2);
plateRadii = zeros(numExp,1);
plateCenters = zeros(numExp,2);

parfor i = 1:numExp
    % Get video object, # of frames, and first frame
    vidObject = VideoReader(videoFile{i});
    numFrames(i) = vidObject.NumFrames;
    firstFrame = rgb2gray(read(vidObject,1));
    
    % Use histogram to get threshold for image
    [counts,~] = imhist(firstFrame,256);
    thresh = otsuthresh(counts);
    
    % Use threshold with high tolerance to get mask for valid tracks (i.e.
    % any tracks in a shadow or under sharpie will not be counted)
    firstFrameBW = imbinarize(firstFrame,thresh*1.15);
    objects = regionprops(~bwlabel(firstFrameBW),{'Area',...
        'MajorAxisLength','BoundingBox','Image','Orientation',...
        'Extrema','FilledImage','MinorAxisLength'});
    objects([objects.Area] < 250) = []; % only get objects with 250 pixels+
    videoMask = zeros(xPix,yPix);
    for o = 1:length(objects)
        objectMask = zeros(xPix,yPix);
        ind = ceil(objects(o).BoundingBox);
        objectMask(ind(2)+[0:ind(4)-1],ind(1)+[0:ind(3)-1]) = ...
            objects(o).FilledImage;
        videoMask = videoMask + objectMask;
    end
    firstFrames(i,:,:) = firstFrame; % image of first frame
    firstFramesBW(i,:,:) = firstFrameBW; % BW image
    videoMasks(i,:,:) = double(videoMask > 0); % mask of all sharpie objects
    
    % Get line object which marks where Copper was placed. This will serve
    % as the reference for all distance measurements
    m = find([objects.MajorAxisLength] > 950 & ...
        [objects.MinorAxisLength] < 50,1);
    if isempty(m)
        % If the line is not dark enough, increase the threshold to get a 
        % thicker line that spans the whole plate
        if strcmp(videoFile{i},'fd_cu_nod_5s_053117.avi')
            firstFrameBW = imbinarize(firstFrame,thresh*1.25);
        elseif strcmp(videoFile{i},'fd_cu_d_2f_052517.avi')
            firstFrameBW = imbinarize(firstFrame,thresh*1.35);
        else
            firstFrameBW = imbinarize(firstFrame,thresh*1.5);
        end
        objects = regionprops(~bwlabel(firstFrameBW),{'Area',...
            'MajorAxisLength','BoundingBox','Image','Orientation',...
            'Extrema','FilledImage','MinorAxisLength'});
        m = find([objects.MajorAxisLength] > 950 & ...
            [objects.MinorAxisLength] < 50,1);
    end
    ind = ceil(objects(m).BoundingBox);
    lineMask = zeros(xPix,yPix);
    lineMask(ind(2)+[0:ind(4)-1],ind(1)+[0:ind(3)-1]) = ...
        objects(m).FilledImage;
    lineMasks(i,:,:) = lineMask; % mask of sharpie line
    lineOrients(i,:) = objects(m).Orientation; % angle of line (90 is perfect)
    lineExtrema(i,:,:) = objects(m).Extrema; % outside points of the line
    
    % Get arc object(s) for the edge of the plate
    firstFrameBW = imbinarize(firstFrame,thresh*0.5);
    objects = regionprops(~bwlabel(firstFrameBW),{'Area',...
        'MajorAxisLength','BoundingBox','Image','Orientation',...
        'Extrema','FilledImage','MinorAxisLength','Eccentricity'});
    if strcmp(videoFile{i},'wf_nocu_d_3f_053117.avi')
        objects = objects([1,3]);
    else
        bounds = [objects.BoundingBox];
        objects = objects([objects.MajorAxisLength] > 70 & ...
            abs([objects.Orientation]) > 59.6 & ...
            abs([objects.Orientation]) < 70.1 & ...
            bounds(1:4:end) < 1);
    end
    arcMask = zeros(xPix,yPix);
    for o = 1:length(objects)
        ind = ceil(objects(o).BoundingBox);
        arcMask(ind(2)+[0:ind(4)-1],ind(1)+[0:ind(3)-1]) = ...
            objects(o).FilledImage;
    end
    arcMasks(i,:,:) = arcMask; % mask of edge of plate
    
    % Fit a circle to arcMask and get radius and centerpoint.
    [x,y] = find(arcMask);
    n=length(x);  xx=x.*x; yy=y.*y; xy=x.*y;
    A=[sum(x) sum(y) n;sum(xy) sum(yy) sum(y);sum(xx) sum(xy) sum(x)];
    B=[-sum(xx+yy) ; -sum(xx.*y+yy.*y) ; -sum(xx.*x+xy.*y)];
    a=A\B;
    xc = -.5*a(1);
    yc = -.5*a(2);
    plateRadii(i)  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
    plateCenters(i,:) = -0.5*a(1:2);
end
numFrame = max(numFrames);

%% Import ScaleBar Video

vidObject = VideoReader([dataPath,'scaleBar.avi']);
scaleFrame = read(vidObject,1);

% Use histogram to get threshold for image
[counts,~] = imhist(scaleFrame,256);
thresh = otsuthresh(counts);

% Get arc object(s) for the edge of the plate
scaleFrameBW = imbinarize(scaleFrame,thresh*1.65);
objects = regionprops(~bwlabel(scaleFrameBW),{'Area',...
    'MajorAxisLength','BoundingBox','Image','Orientation',...
    'Extrema','FilledImage','MinorAxisLength','Eccentricity'});
objects = objects(54);
arcMask = zeros(xPix,yPix); % mask of edge of plate
for o = 1:length(objects)
    ind = ceil(objects(o).BoundingBox);
    arcMask(ind(2)+[0:ind(4)-1],ind(1)+[0:ind(3)-1]) = ...
        objects(o).FilledImage;
end

% Fit a circle to arcMask and get radius and centerpoint.
[x,y] = find(arcMask);
n=length(x);  xx=x.*x; yy=y.*y; xy=x.*y;
A=[sum(x) sum(y) n;sum(xy) sum(yy) sum(y);sum(xx) sum(xy) sum(x)];
B=[-sum(xx+yy) ; -sum(xx.*y+yy.*y) ; -sum(xx.*x+xy.*y)];
a=A\B;
xc = -.5*a(1);
yc = -.5*a(2);
scaleRadius  =  sqrt((a(1)^2+a(2)^2)/4-a(3));
scaleCenter = -0.5*a(1:2);

% Get scale bar edge
scaleFrameBW = imbinarize(scaleFrame,thresh);
objects = regionprops(~bwlabel(scaleFrameBW),'Extrema');
endPoints = objects.Extrema([1,6],:);
edgePosition = round(linspace(endPoints(1,1),endPoints(2,1),xPix));

% Get scale bar dashes
scaleTicks = scaleFrame(sub2ind([xPix,yPix],[1:xPix],edgePosition+90));
[pks,locs] = findpeaks(-double(scaleTicks),'MinPeakHeight',-30,...
    'MinPeakProminence',15,'MinPeakDistance',20);
scalePxPerMm = mean(diff(locs)); % px/mm

% Calculate plate radius in mm
scaleRadiusMm = scaleRadius/scalePxPerMm;
scaleRadiusCm = scaleRadiusMm/10;

%% Get scaling for each plate

% Some plates did not have enough contrast in the line to get a good fit
% for the plate. These will be excluded for determining the scaling
badFit = true(size(camera)); 
badFit([7,17,19,25,28,32,34,42]) = 0;

% Calculate scale of each camera (cm/px)
scalePx = [median(plateRadii(badFit & camera==0)),...
    median(plateRadii(badFit & camera==1))];
scaleCm = scaleRadiusCm./scalePx; % cm/px
scaleExpPx = camera; 
scaleExpPx(camera==0) = scalePx(1); scaleExpPx(camera==1) = scalePx(2);
scaleExpCm = scaleRadiusCm./scaleExpPx; % cm/px

%% Save experiment info in a table

expNum = [1:numExp]';
info = table(expNum,hunger,odor,videoFile,camera,wormLabFile,wormLabScale,...
    numFrames,plateRadii,plateCenters,lineOrients,scaleExpPx,scaleExpCm);

%% Load X,Y Positions from Exported .xls Files from WormLab

% Load all tracks into an array and store indices of tracks corresponding
% to each video
wormLabRaw = []; % all X,Y positions
wormLabNumTracks = zeros(numExp,1); % indices of tracks per video
for i = 1:numExp
    temp = readmatrix(wormLabFile{i},'Range',[6,1]);
    missingTime = numFrame - size(temp,1);
    if missingTime > 0
        temp = [temp;nan(missingTime,size(temp,2))];
    end
    wormLabRaw = [wormLabRaw,temp(1:numFrame,3:end)./wormLabScale(i)];
    wormLabNumTracks(i) = size(wormLabRaw,2);
end
wormLabNumTracks = [0;wormLabNumTracks./2]; % edges of tracks in each file
wormLabRaw(wormLabRaw > 1024) = NaN;
wormLabRaw(wormLabRaw <= 0) = NaN;

% Extract X and Y positions from loaded data
positions(1,:,:) = wormLabRaw(:,1:2:end); % X position
positions(2,:,:) = wormLabRaw(:,2:2:end); % Y position

% Convert X,Y into 1D index
positionsInd = sub2ind([1024 1024],ceil(reshape(wormLabRaw(:,1:2:end),1,[])),...
    ceil(reshape(wormLabRaw(:,2:2:end),1,[]))); % index of X,Y positions
positionsInd = reshape(positionsInd,numFrame,[]); % # frames x # tracks
positionsNaN = isnan(positionsInd); % 1 = frames w/out tracks, 0 = w/ tracks
positionsOnes = positionsInd; positionsOnes(positionsNaN) = 1;

%% Clean data based on sharpie marks, track length, area travelled

positionsClean = positions;

% Remove tracks that overlap with mask of Sharpie/plate
j = 1;
for i = 1:numExp
    % Get wormlab file and corresponding video file
    [~,expName,~] = fileparts(wormLabFile{i});
    vidNum = find(contains(videoFile,expName));
    
    % Get mask and dilate it to include the nearest 10 pixels
    mask = flip(squeeze((videoMasks(i,:,:) + lineMasks(i,:,:))>0))';
    mask = imdilate(mask,strel('square',10));
    
    % Remove overlapping tracks
    while j > wormLabNumTracks(i) && j <= wormLabNumTracks(i+1)
        overlap = mask(positionsOnes(:,j)); % 1 = overlap, 2 = no overlap
        positionsClean(:,overlap>0,j) = NaN; % remove points w/ overlap
        j = j+1;
    end
end

% Delete tracks that are too short (time) or small (area/length)
minFrames = 10*frameRate; % minimum track length = 10 seconds
minArea = 30; % minimum area traveled (pixels^2)
min1D = 10; % minimum distance traveled in X or Y (pixels)

% Get number of frames that each track spans
numTracks = sum(~positionsNaN,1); % # of frames in each track

% Get area traveled (bounding box around track)
minXY = squeeze(min(positions,[],2,'omitnan'));
maxXY = squeeze(max(positions,[],2,'omitnan'));
dX = maxXY(1,:)-minXY(1,:);
dY = maxXY(2,:)-minXY(2,:);
dArea = dX.*dY;

% Get indices of tracks and remove them
removeTracks = find(dX < minArea | dX < min1D | dY < min1D | ...
    numTracks < minFrames);
positionsClean(:,:,removeTracks) = NaN;
xPositions = squeeze(positionsClean(1,:,:));
yPositions = squeeze(positionsClean(2,:,:));

%% Extract frame, experiment #, and positions into a table

[frame,track] = find(~isnan(xPositions));
expNum = nan(size(frame));
for i = 1:numExp
    ind = find(track >= wormLabNumTracks(i)+1 & ...
        track <= wormLabNumTracks(i+1));
    expNum(ind) = i;
end
data = table(expNum,frame,track);
ind = find(~isnan(xPositions));
data.xPosition = xPositions(ind);
data.yPosition = yPositions(ind);

% If jump in track, replace one value with NaN
skip = find([0;~diff(data.frame)==1]);
data.xPosition(skip) = NaN; data.yPosition(skip) = NaN;

%% Calculate # of worms in each experiment

% Find number of worm tracks in each frame
numWorms = zeros(numExp,numFrame); % # of worms in each frame
parfor i = 1:numExp
    for j = 1:numFrame
        numWorms(i,j) = sum(data.expNum == i & data.frame == j);
    end
end

% Find the max, min, mean number of possible worms in each experiment 
[maxWorms,bestFrame] = max(numWorms,[],2);
meanWorms = mean(numWorms,2);
[minWorms,worstFrame] = min(numWorms,[],2);

% Add worm info to table
info.meanWorms = meanWorms; % mean number of tracks in a frame
info.minWorms = minWorms; % minimum number of tracks in a frame
info.maxWorms = maxWorms; % maximum number of tracks in a frame

%% Get distance to line (pixels)

lineEdges = [mean(lineExtrema(:,[1,2],:),2),... % top sharpie mark
    mean(lineExtrema(:,[5,6],:),2)]; % bottom sharpie mark
info.lineX = squeeze(lineEdges(:,:,1));
info.lineY = squeeze(lineEdges(:,:,2));
data.lineDistancePx = nan(size(data,1),1);
data.lineDistanceCm = nan(size(data,1),1);

for i = 1:numExp
    ind = find(data.expNum == i);
    data.scaleCmPerPx(ind) = info.scaleExpCm(i);
    
    % Use linear interpolation to get X position of line for each track Y position
    linePosition = interp1q(info.lineY(i,:)',info.lineX(i,:)',...
        data.yPosition(ind));
    
    % Calculate distance from each track position to the line
    data.lineDistancePx(ind) = data.xPosition(ind) - linePosition;
    
    % Convert to real units (Cm)
    data.lineDistanceCm(ind) = data.lineDistancePx(ind) * info.scaleExpCm(i);
end

%% Calculate Path Length and Velocity (pixels/frame)

% Use 2 second moving time window to calculate Euclidean distance
winSize = 2*frameRate; % 2 seconds
dX = data.xPosition(1+winSize:end) - data.xPosition(1:end-winSize);
dY = data.yPosition(1+winSize:end) - data.yPosition(1:end-winSize);
data.pathLengthPx = [nan(ceil(winSize/2),1);...
    sqrt(dX.^2 + dY.^2);...
    nan(floor(winSize/2),1)]./winSize;

% Replace first and last second of each track with NaN
data.pathLengthPx(find(diff(data.track))+[-2:3]) = NaN;
data.velocityMicronPerS = 10000*frameRate*data.pathLengthPx.*data.scaleCmPerPx;

%% Calculate probability of location and mean velocity as a function of distance from line

% Start parallel processing
if isempty(gcp('nocreate'))
    pool = parpool(24); % start parellel pool
end

% Get indices of all data in each 1 mm bin
binSize = 0.1; % 1 mm (cm)
distRange = [-4.4:binSize:1.1]-(binSize/2); % 1 mm bins (cm)
numDistBins = length(distRange)-1;
locationSum = zeros(numExp,numDistBins);
velocityMean = zeros(numExp,numDistBins);
velocitySD = zeros(numExp,numDistBins);
parfor i = 1:numExp
    for j = 1:numDistBins
        ind = find(data.expNum == i & ...
            data.lineDistanceCm >= distRange(j) & ...
            data.lineDistanceCm < distRange(j+1));
        
        % Get total number of worms in each mm bin
        locationSum(i,j) = length(ind);
        
        % Get mean and standard deviation of path length of tracks in each mm bin
        velocityMean(i,j) = mean(data.velocityMicronPerS(ind),'omitnan');
        velocitySD(i,j) = std(data.velocityMicronPerS(ind),'omitnan');
    end
end
% Get probability of worms residing in each mm bin
locationProb = locationSum./sum(locationSum,2);

% Add data to table
analysis = info(:,1:3);
analysis.locationSum = locationSum;
analysis.locationProb = locationProb;
analysis.velocityMean = velocityMean;
analysis.velocitySD = velocitySD;

%% Calculate cummulative tracks past the line as a function of time (15 min bins)

timePoints = [15;30;45];
numTimeBins = length(timePoints);
crossingsCum = zeros(numExp,numTimeBins);
tracksCum = zeros(numExp,numTimeBins);
for i = 1:numExp
    for j = 1:numTimeBins
        % Get cumulative sum of unique tracks past the line
        ind = find(data.expNum == i & ...
            data.frame <= timePoints(j)*60*frameRate & ...
            data.lineDistanceCm >= 0.05);
        crossingsCum(i,j) = length(unique(data.track(ind)));
        
        % Get cumulative sum of unique tracks at each time point
        ind = find(data.expNum == i & ...
            data.frame <= timePoints(j)*60*frameRate);
        tracksCum(i,j) = length(unique(data.track(ind)));
    end
end

% Get fraction of cumulative unique tracks past the line
crossingsCumFrac = crossingsCum./tracksCum;

% Add data to table
analysis.tracksCum = tracksCum;
analysis.crossingsCum = crossingsCum;
analysis.crossingsCumFrac = crossingsCumFrac;

%% Calculate summary statistics (mean + std) across conditions

numExpPerCondition = nan(length(hungerFolders),length(odorFolders));
trackStats = nan(length(hungerFolders),length(odorFolders),numFrame,3);
locationStats = nan(length(hungerFolders),length(odorFolders),numDistBins,3);
velocityStats = nan(length(hungerFolders),length(odorFolders),numDistBins,3);
crossingStats = nan(length(hungerFolders),length(odorFolders),numTimeBins,3);

for i = 1:length(hungerFolders)
    for j = 1:length(odorFolders)
        % Get indices and number of experiments matching hunger & odor condition
        ind = find(strcmp(info.hunger,hungerFolders{i}) & ...
            strcmp(info.odor,odorFolders{j}));
        numExpPerCondition(i,j) = length(ind);
        
        % Probability of worms residing in each mm bin
        locationStats(i,j,:,1) = mean(locationProb(ind,:),1,'omitnan');
        locationStats(i,j,:,2) = std(locationProb(ind,:),1,'omitnan');
        
        % Velocity (um/second) of worms in each experiment as a
        % function of distance from the line
        velocityStats(i,j,:,1) = mean(velocityMean(ind,:),1,'omitnan');
        velocityStats(i,j,:,2) = std(velocityMean(ind,:),1,'omitnan');

        % Number of tracks (worms) in each frame
        trackStats(i,j,:,1) = mean(numWorms(ind,:),1,'omitnan');
        trackStats(i,j,:,2) = std(numWorms(ind,:),1,'omitnan');
    end
end

% Calculate Standard Error of the Mean (SEM)
locationStats(:,:,:,3) = locationStats(:,:,:,2)./sqrt(numExpPerCondition);
velocityStats(isnan(velocityStats)) = 0;
velocityStats(:,:,:,3) = velocityStats(:,:,:,2)./sqrt(numExpPerCondition);
trackStats(:,:,:,3) = trackStats(:,:,:,2)./sqrt(numExpPerCondition);

%% Export crossings for Molly to graph in Prism

for j = 1:length(odorFolders)
    prismTable = table(timePoints);
    for i = 1:length(hungerFolders)
        % Get indices and number of experiments matching hunger & odor condition
        ind = find(strcmp(info.hunger,hungerFolders{i}) & ...
            strcmp(info.odor,odorFolders{j}));
        prismTable.(hungerFolders{i}) = crossingsCumFrac(ind,:)';
    end
    writetable(prismTable,[figurePath,'CumalitiveCrossings_',odorFolders{j},'.xlsx']);
end

%% Export p(location) and velocity for Molly to graph in Prism

for j = 1:length(odorFolders)
    prismTable = table((distRange(1:end-1)+binSize/2)');
    for i = 1:length(hungerFolders)
        % Get indices and number of experiments matching hunger & odor condition
        ind = find(strcmp(info.hunger,hungerFolders{i}) & ...
            strcmp(info.odor,odorFolders{j}));
        prismTable.(hungerFolders{i}) = locationProb(ind,:)';
    end
    writetable(prismTable,[figurePath,'ProbabilityLocation_',odorFolders{j},'.xlsx']);
end

for j = 1:length(odorFolders)
    prismTable = table((distRange(1:end-1)+binSize/2)');
    for i = 1:length(hungerFolders)
        % Get indices and number of experiments matching hunger & odor condition
        ind = find(strcmp(info.hunger,hungerFolders{i}) & ...
            strcmp(info.odor,odorFolders{j}));
        prismTable.(hungerFolders{i}) = velocityMean(ind,:)';
    end
    writetable(prismTable,[figurePath,'Velocity_',odorFolders{j},'.xlsx']);
end

%% Export info and analysis tables

writetable(info,[figurePath,'info.xlsx']);
writetable(analysis,[figurePath,'analysis.xlsx']);
        
%% Plot Colors

blue = [0,144,178]/255;
gray = [0.8,0.8,0.8];
cyan = [0,1,1];

%% Plot p(Location) and Velocity vs. distance from line
xData =  distRange(1:end-1)+binSize/2;

for i = 1:2
    switch i
        case 1 
            currentVar = locationStats;
            yLimits = [0 0.08];
            yLabel = 'p(Location)';
            fileLabel = 'Location';
        case 2
            currentVar = velocityStats;
            yLimits = [0 250];
            yLabel = 'Velocity (\mum/s)';
            fileLabel = 'Velocity';
    end
    
    for j = 1:3 % odor state
        yData = permute(currentVar(:,j,:,:),[3,1,4,2]);
        figure;
        hold on
        
        % Plot SEM for each condition
        patch([xData,flip(xData)],[yData(:,1,1)'+yData(:,1,3)',...
            flip(yData(:,1,1)'-yData(:,1,3)')],'k','FaceAlpha',0.3,...
            'LineStyle','none');
        patch([xData,flip(xData)],[yData(:,2,1)'+yData(:,2,3)',...
            flip(yData(:,2,1)'-yData(:,2,3)')],blue,'FaceAlpha',0.3,...
            'LineStyle','none');
        
        % Plot Mean for each condition
        p1 = plot(xData,yData(:,1,1),'k','LineWidth',1); % WF
        p2 = plot(xData,yData(:,2,1),'Color',blue,'LineWidth',1); % FD
        
        % Cover the copper barrier +/- 0.5 mm and plot a line
        patch([-0.05 -0.05 0.05 0.05],[yLimits,flip(yLimits)],...
            'w','FaceAlpha',1,'LineStyle','none');
        plot([0 0],yLimits,'Color',cyan,'LineWidth',3);
        
        xlim([-4 0.5])
        ylim(yLimits)
        xlabel('Distance from Barrier (cm)')
        ylabel(yLabel)
        title(odorFolders{j})
        legend([p1,p2],{['WF (n = ',num2str(numExpPerCondition(1,j)),')'],...
            ['FD (n = ',num2str(numExpPerCondition(2,j)),')']},...
            'box','off','location','NorthWest')
        set(gca,'FontSize',14,'FontName','Arial')
        
        % Save figure as .jpg and .svg then close
        saveas(gcf,[figurePath,fileLabel,'_',odorFolders{j},'.jpg'])
        saveas(gcf,[figurePath,fileLabel,'_',odorFolders{j},'.svg'])
        exportgraphics(gca,[figurePath,fileLabel,'_',odorFolders{j},'.pdf'],...
            'BackgroundColor','none','ContentType','vector')
        close(gcf);
    end
end

%% Plot paths

c = jet(45);
frameBins = round(linspace(1,numFrame,45));

parfor expNum = 1:numExp
    [~,expName,~] = fileparts(wormLabFile{expNum});
    vidNum = find(contains(videoFile,expName));
    expTracks = unique(data.track(data.expNum == expNum));
    
    figure('Position',[0 0 1000 1000])
    
    % Show sharpie marks and shadows in gray
    imshow(squeeze(flip((~videoMasks(vidNum,:,:)*0.3+0.7))));
    hold on
    
    % Plot copper barrier in cyan
    plot(lineEdges(expNum,[2,1],1),lineEdges(expNum,:,2),...
        'Color',cyan,'LineWidth',5)
    
    % Plot edge of the plate in black
    p = nsidedpoly(1000, 'Center', [plateCenters(expNum,2),xPix-plateCenters(expNum,1)],...
        'Radius', plateRadii(expNum));
    plot(p,'FaceAlpha',0,'LineWidth',2)
    
    % Plot tracks colored by time
    for j = 1:length(expTracks)
        for k = 1:length(frameBins)-1
            ind = find(data.track == expTracks(j) & ...
                data.frame >= frameBins(k) & data.frame < frameBins(k+1));
            plot(data.xPosition(ind),data.yPosition(ind),...
                'Color',c(k,:),'LineWidth',0.7);
        end
    end
    
    xlim([1 1024])
    ylim([1 1024])
    xticks([])
    yticks([])
    set(gcf,'Toolbar','none')
    
    % Plot 1 cm scale bar
    scaleBar = plot(xPix*16/25 + [0 scaleExpPx(expNum)/5],...
        yPix*24/25 + [0 0],'Color','k','LineWidth',3);
    scaleText = text(xPix*16/25 + scaleExpPx(expNum)/10,...
        yPix*23/25,'1 cm','HorizontalAlignment','center',...
        'Color','k','FontName','Arial','FontSize',10);
    
    % Label figure with max number of worms
    text(xPix*24.5/25,yPix*0.5/25,['n = ',num2str(maxWorms(expNum))],...
        'HorizontalAlignment','right','VerticalAlignment','top','FontName','Arial')
    set(gca,'FontName','Arial','Visible','off');
    
    % Save figure as .jpg and .svg then close
    saveas(gcf,[figurePath,'Tracks_Exp',num2str(expNum,'%.2d'),'.jpg'])
    saveas(gcf,[figurePath,'Tracks_Exp',num2str(expNum,'%.2d'),'.svg'])
    exportgraphics(gca,[figurePath,'Tracks_Exp',num2str(expNum,'%.2d'),'.pdf'],...
            'BackgroundColor','none','ContentType','vector')
    close(gcf);
end

%% Export color map

figure('Position',[0 0 1000 1000])
colormap(jet)
c = colorbar('west','Ticks',0:5:45);
c.Label.String = 'Time (min)';
c.Label.Rotation = -90;
c.Label.VerticalAlignment = 'bottom';
caxis([0 45])
set(gca,'FontName','Arial','Visible','off');
saveas(gcf,[figurePath,'Colorbar.jpg'])
saveas(gcf,[figurePath,'Colorbar.svg'])
exportgraphics(gca,[figurePath,'Colorbar.pdf'],...
    'BackgroundColor','none','ContentType','vector')
close(gcf);

%% Create Figure Supplement for P(Location)

expNum = 27;

[~,expName,~] = fileparts(wormLabFile{expNum});
vidNum = find(contains(videoFile,expName));
expTracks = unique(data.track(data.expNum == expNum));
c = jet(length(expTracks));
rng(5)
colorOrder = randperm(length(c));

figure('Position',[0 0 500 500])
hold on

% Plot mm bin boundary in gray
mmScale = 0.1/scaleExpCm(expNum);
xBounds = [1 4]+5.5;
yBounds = [1 4]+27;
plot(mmScale*(xBounds(1)+1).*[1,1],mmScale.*yBounds,...
    'Color',gray,'LineWidth',2)
plot(mmScale*(xBounds(2)-1).*[1,1],mmScale.*yBounds,...
    'Color',gray,'LineWidth',2)

% Plot tracks in unique colors
for j = 1:length(expTracks)
    ind = find(data.track == expTracks(j));
    plot(data.xPosition(ind),data.yPosition(ind),...
        'Marker','o','Color',c(colorOrder(j),:),'LineWidth',1,...
        'MarkerSize',4,'MarkerFaceColor','w',...
        'MarkerEdgeColor',c(colorOrder(j),:));
    ind = find(data.track == expTracks(j) & ...
        data.xPosition >= mmScale.*(xBounds(1)+1) & ...
        data.xPosition <= mmScale.*(xBounds(2)-1));
    plot(data.xPosition(ind),data.yPosition(ind),...
        'Marker','o','MarkerFaceColor',c(colorOrder(j),:),...
        'MarkerEdgeColor',c(colorOrder(j),:),...
        'MarkerSize',4,'LineStyle','none');
end
bin1 = sum(data.expNum == expNum & ...
    data.xPosition >= mmScale.*(xBounds(1)) & ...
    data.xPosition < mmScale.*(xBounds(1)+1) & ...
    data.yPosition >= mmScale.*yBounds(1) & ...
    data.yPosition < mmScale.*yBounds(2));
bin2 = sum(data.expNum == expNum & ...
    data.xPosition >= mmScale.*(xBounds(1)+1) & ...
    data.xPosition <= mmScale.*(xBounds(2)-1) & ...
    data.yPosition >= mmScale.*yBounds(1) & ...
    data.yPosition < mmScale.*yBounds(2));
bin3 = sum(data.expNum == expNum & ...
    data.xPosition > mmScale.*(xBounds(2)-1) & ...
    data.xPosition <= mmScale.*xBounds(2) & ...
    data.yPosition >= mmScale.*yBounds(1) & ...
    data.yPosition < mmScale.*yBounds(2));

xlim(mmScale.*xBounds)
ylim(mmScale.*yBounds)
xticks([])
yticks([])
set(gcf,'Toolbar','none')
set(gca,'FontName','Arial','Visible','off');

% Save figure as .jpg and .svg then close
saveas(gcf,[figurePath,'ProbabilityLocation.jpg'])
saveas(gcf,[figurePath,'ProbabilityLocation.svg'])
exportgraphics(gca,[figurePath,'ProbabilityLocation.pdf'],...
    'BackgroundColor','none','ContentType','vector')
close(gcf);

%% Create Figure Supplement for Velocity

expNum = 27;

[~,expName,~] = fileparts(wormLabFile{expNum});
vidNum = find(contains(videoFile,expName));
expTracks = unique(data.track(data.expNum == expNum));
c = jet(length(expTracks));
rng(5)
colorOrder = randperm(length(c));

figure('Position',[0 0 500 500])
hold on

% Plot mm bin boundary in gray
mmScale = 0.1/scaleExpCm(expNum);
xBounds = [1 3]+9;
yBounds = [1 3]+36;

% Plot tracks in unique colors
ind = find(data.track == expTracks(168));
plot(data.xPosition(ind),data.yPosition(ind),...
    'Marker','o','Color','k','LineWidth',1,...
    'MarkerSize',4,'MarkerFaceColor','k',...
    'MarkerEdgeColor','k');
xlim(mmScale.*xBounds)
ylim(mmScale.*yBounds)
xticks([])
yticks([])
set(gcf,'Toolbar','none')
set(gca,'FontName','Arial','Visible','off');

% Plot 1 cm scale bar
scaleBar = plot(mmScale.*(xBounds(1) + [1.4 1.8]),...
    mmScale.*(yBounds(1) + [0.9 0.9]),'Color','k','LineWidth',3);
scaleText = text(mmScale.*(xBounds(1) + 1.6),...
    mmScale.*(yBounds(1) + 0.8),'400 \mum','HorizontalAlignment','center',...
    'Color','k','FontName','Arial','FontSize',10);

% Save figure as .jpg and .svg then close
saveas(gcf,[figurePath,'Velocity.jpg'])
saveas(gcf,[figurePath,'Velocity.svg'])
exportgraphics(gca,[figurePath,'Velocity.pdf'],...
    'BackgroundColor','none','ContentType','vector')
close(gcf);

%% Create Figure Supplement for Barrier Crossings

expNum = 27;

[~,expName,~] = fileparts(wormLabFile{expNum});
vidNum = find(contains(videoFile,expName));
expTracks = unique(data.track(data.expNum == expNum));
c = jet(length(expTracks));
rng(5)
colorOrder = randperm(length(c));

figure('Position',[0 0 500 500])

% Show sharpie marks and shadows in gray
imshow(squeeze(flip((~videoMasks(vidNum,:,:)*0.3+0.7))));
hold on

% Plot copper barrier in cyan
plot(lineEdges(expNum,[2,1],1),lineEdges(expNum,:,2),...
    'Color',cyan,'LineWidth',5)

% Plot edge of the plate in black
p = nsidedpoly(1000, 'Center', [plateCenters(expNum,2),xPix-plateCenters(expNum,1)],...
    'Radius', plateRadii(expNum));
plot(p,'FaceAlpha',0,'LineWidth',2)

% Plot tracks with unique color
for j = 1:length(expTracks)
    ind = find(data.track == expTracks(j) & ...
        data.frame >= 0*frameRate*60 & ...
        data.frame <= 15*frameRate*60);
    if ~isempty(ind)
        plot(data.xPosition(ind),data.yPosition(ind),...
            'Color',c(colorOrder(j),:),'LineWidth',0.7);
        scatter(data.xPosition(ind(1)),data.yPosition(ind(1)),100,...
            'MarkerFaceColor',c(colorOrder(j),:),...
            'MarkerEdgeColor',c(colorOrder(j),:));
        text(data.xPosition(ind(1)),data.yPosition(ind(1)),num2str(j),...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'Color','w','FontSize',6,'FontWeight','bold');
    end
end
% Label tracks
for j = 1:length(expTracks)
    ind = find(data.track == expTracks(j) & ...
        data.frame >= 0*frameRate*60 & ...
        data.frame <= 15*frameRate*60);
    if ~isempty(ind)
        scatter(data.xPosition(ind(1)),data.yPosition(ind(1)),100,...
            'MarkerFaceColor',c(colorOrder(j),:),...
            'MarkerEdgeColor',c(colorOrder(j),:));
        text(data.xPosition(ind(1)),data.yPosition(ind(1)),num2str(j),...
            'HorizontalAlignment','center','VerticalAlignment','middle',...
            'Color','w','FontSize',6,'FontWeight','bold','FontName','Arial');
    end
end

xlim([1 1024])
ylim([1 1024])
xticks([])
yticks([])
set(gcf,'Toolbar','none')

% Plot 1 cm scale bar
scaleBar = plot(xPix*16/25 + [0 scaleExpPx(expNum)/5],...
    yPix*24/25 + [0 0],'Color','k','LineWidth',3);
scaleText = text(xPix*16/25 + scaleExpPx(expNum)/10,...
    yPix*23/25,'1 cm','HorizontalAlignment','center',...
    'Color','k','FontName','Arial','FontSize',10);

% Label figure with max number of worms
text(xPix*24.5/25,yPix*0.5/25,['n = ',num2str(maxWorms(expNum))],...
    'HorizontalAlignment','right','VerticalAlignment','top','FontName','Arial')
set(gca,'FontName','Arial','Visible','off');

% Save figure as .jpg and .svg then close
saveas(gcf,[figurePath,'BarrierCrossings.jpg'])
saveas(gcf,[figurePath,'BarrierCrossings.svg'])
exportgraphics(gca,[figurePath,'BarrierCrossings.pdf'],...
    'BackgroundColor','none','ContentType','vector')
close(gcf);