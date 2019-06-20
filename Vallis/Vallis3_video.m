% Demo to create a movie file from a Gaussian and then optionally save it to disk as an avi video file.

%==============================================================================================
% Initialization code
clear all;
clc;
workspace;
numberOfFrames = 61;
hFigure = figure;

% Set up the movie structure.
% Preallocate movie, which will be an array of structures.
% First get a cell array with all the frames.
allTheFrames = cell(numberOfFrames,1);
vidHeight = 344;
vidWidth = 446;
allTheFrames(:) = {zeros(vidHeight, vidWidth, 3, 'uint8')};
% Next get a cell array with all the colormaps.
allTheColorMaps = cell(numberOfFrames,1);
allTheColorMaps(:) = {zeros(256, 3)};
% Now combine these to make the array of structures.
myMovie = struct('cdata', allTheFrames, 'colormap', allTheColorMaps);
% Create a VideoWriter object to write the video out to a new, different file.
% writerObj = VideoWriter('problem_3.avi');
% open(writerObj);
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
% set(gcf, 'renderer', 'zbuffer');

%==============================================================================================
% Create the movie.
% After this loop starts, BE SURE NOT TO RESIZE THE WINDOW AS IT'S SHOWING THE FRAMES, or else you won't be able to save it.
TT = linspace(-20,8,numberOfFrames);

transi = 8*12*30*24*60*60;  % in seconds
% ttot   = 1000*12*30;  % run length in days 
ttot = 42*12*30*24*60*60;
h     = 12*60*60;          % in seconds

%% initial conditions
u_init = -0.5;
T_w_init = 18;
T_e_init = 12;

for j = 1:numberOfFrames
%% Choice of flag
flag = 'none';




% parameters from Vallis
A = 1/12/30/24/60/60;
B = 2;
C = 1/4/30/24/60/60;
U = -0.45;
T_deep = TT(j);
T = 12;
dx = 7500*1000;
omega = 2*pi/(12*30*24*60*60);


%% derived parameters
no_steps_transi = transi/h;
no_steps = ttot/h;

%% initialise variables transient
y_vec_transi = zeros(3,no_steps_transi+no_steps);
t_vec_transi = zeros(1,no_steps_transi+no_steps);
y_vec = zeros(3,no_steps);
t_vec = zeros(1,no_steps);
y_init = [u_init;T_w_init;T_e_init];
y_vec_transi(:,1) = y_init;


%% integrate transient and real in one step
f1 = @(t,x) [B*(x(3)-x(2))/2/dx-C*(x(1)-U);...
           x(1)*(T_deep-x(3))/2/dx-A*(x(2)-T);...
           x(1)*(x(2)-T_deep)/2/dx-A*(x(3)-T)];
if strcmp(flag,'none')
    for i = 2:no_steps_transi+no_steps
        k1 = f1(t_vec_transi(i-1),y_vec_transi(:,i-1));
        k2 = f1(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k1/2);
        k3 = f1(t_vec_transi(i-1)+h/2,y_vec_transi(:,i-1)+h*k2/2);
        k4 = f1(t_vec_transi(i-1)+h,y_vec_transi(:,i-1)+h*k3);
        y_vec_transi(:,i) = y_vec_transi(:,i-1)+h*(k1+2*k2+2*k3+k4)/6;
        t_vec_transi(i) = h*(i-1);
    end
end


y_vec = y_vec_transi(:,no_steps_transi+1:end);
t_vec = t_vec_transi(no_steps_transi+1:end);

plot(y_vec(1,:),y_vec(3,:)-y_vec(2,:))
xlabel('u')
ylabel('difference')
axis([-15 15 -20 20])
title(['T is', num2str(T_deep)])

% plot(t_vec/60/60/24/30/12,y_vec(1,:))
% xlabel('Time in year')
% ylabel('u')
% title(['T is', num2str(T_deep)])

drawnow;
	thisFrame = getframe(gcf);
	% Write this frame out to a new video file.
%  	writeVideo(writerObj, thisFrame);
	myMovie(j) = thisFrame;
end

% for frameIndex = 1 : numberOfFrames
% 	z = exp(-(x-t(frameIndex)).^2-(y-t(frameIndex)).^2);
% 	cla reset;
% 	% Enlarge figure to full screen.
% % 	set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
% 	surf(x,y,z);
% 	axis('tight')
% 	zlim([0, 1]);
% 	caption = sprintf('Frame #%d of %d, t = %.1f', frameIndex, numberOfFrames, t(frameIndex));
% 	title(caption, 'FontSize', 15);
% end
% close(writerObj);

%==============================================================================================
% See if they want to replay the movie.
message = sprintf('Done creating movie\nDo you want to play it?');
button = questdlg(message, 'Continue?', 'Yes', 'No', 'Yes');
drawnow;	% Refresh screen to get rid of dialog box remnants.
close(hFigure);
if strcmpi(button, 'Yes')
	hFigure = figure;
	% Enlarge figure to full screen.
	% set(gcf, 'Units', 'Normalized', 'Outerposition', [0, 0, 1, 1]);
	title('Varying T deep', 'FontSize', 15);
	% Play the movie.
	movie(myMovie,1,4);
	close(hFigure);
end

%==============================================================================================
% See if they want to save the movie to an avi file on disk.
promptMessage = sprintf('Do you want to save this movie to disk?');
titleBarCaption = 'Continue?';
button = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
if strcmpi(button, 'yes')
	% Get the name of the file that the user wants to save.
	% Note, if you're saving an image you can use imsave() instead of uiputfile().
	startingFolder = pwd;
	defaultFileName = {'*.avi';'*.mp4';'*.mj2'}; %fullfile(startingFolder, '*.avi');
	[baseFileName, folder] = uiputfile(defaultFileName, 'Specify a file');
	if baseFileName == 0
		% User clicked the Cancel button.
		return;
	end
	fullFileName = fullfile(folder, baseFileName);
	% Create a video writer object with that file name.
	% The VideoWriter object must have a profile input argument, otherwise you get jpg.
	% Determine the format the user specified:
	[folder, baseFileName, ext] = fileparts(fullFileName);
	switch lower(ext)
		case '.jp2'
			profile = 'Archival';
		case '.mp4'
			profile = 'MPEG-4';
		otherwise
			% Either avi or some other invalid extension.
			profile = 'Uncompressed AVI';
	end
	writerObj = VideoWriter(fullFileName, profile);
    writerObj.FrameRate = 4;
	open(writerObj);
	% Write out all the frames.
	numberOfFrames = length(myMovie);
	for frameNumber = 1 : numberOfFrames 
	   writeVideo(writerObj, myMovie(frameNumber));
	end
	close(writerObj);
	% Display the current folder panel so they can see their newly created file.
	cd(folder);
	filebrowser;
	message = sprintf('Finished creating movie file\n      %s.\n\nDone with demo!', fullFileName);
	uiwait(helpdlg(message));
else
	uiwait(helpdlg('Done with demo!'));
end







       