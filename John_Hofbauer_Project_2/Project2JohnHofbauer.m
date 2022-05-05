%% Big picture of Project 2
    % Consider same points from diffrent camera
    % Project from 3D to 2D
    % Use triangluuation for reconstruction
    close all
    clc
    %profile on
    

% Initialization of VideoReader for the vue video.
%YOU ONLY NEED TO DO THIS ONCE AT THE BEGINNING
    filenamevue2mp4 = 'Subject4-Session3-24form-Full-Take4-Vue2.mp4';
    vue2video = VideoReader(filenamevue2mp4);

    filenamevue4mp4 = 'Subject4-Session3-24form-Full-Take4-Vue4.mp4';
    vue4video = VideoReader(filenamevue4mp4);

% Load Mocap Data Sets
    load('Subject4-Session3-Take4_mocapJoints.mat');

%%  2.1 Find the good frames (conidence = 1) 
    counter = 1;
    for mocapFnum = 1:length(mocapJoints) % Loop through the frames.
        if sum(mocapJoints(mocapFnum,:,4)) == 12 % Check that all 12 values are 1
            goodFrames(counter).frame = mocapFnum; %#ok<SAGROW> 
            counter = counter + 1;
        end
    end

%% 2.2 Load Camera Parameters AND check P matrix
    load('vue2CalibInfo.mat'); % prenamed vue2
    load('vue4CalibInfo.mat'); % prenamed vue4

    % Check that the Pmat and Kmat were created corrcetly.
    %vue2.orientation
    %vue2.position
    %[U,S,V] = svd(vue2.Pmat)
    %[U,S,V] = svd(vue2.orientation)

    vue2.myp = myp(vue2);
    if sum(vue2.myp,'all') - sum(vue2.Kmat * vue2.Pmat, 'all') < 10
        fprintf("The camra P matrix was calculated corectly\n");
    end
    vue4.myp = myp(vue4);

%% 2.3 Projecting 3D points into 2D pixel locations
    %find teh pseudeo in verso fo the p matris
    % multiply by the in verse of teh px
    starttime = datetime('now');
    
    for frame = 1:1:length(goodFrames)
        mocapFnum = goodFrames(frame).frame;
        xpoints = mocapJoints(mocapFnum,:,1); %array of 12 X coordinates
        ypoints = mocapJoints(mocapFnum,:,2); % Y coordinates
        zpoints = mocapJoints(mocapFnum,:,3); % Z coordinates

        
        % convert to 2d point on left image
        for i = 1:1:12
            %jointi = [xpoints(i); ypoints(i); zpoints(i)];
            goodFrames(frame).mocapJoints(i).X = xpoints(i);
            goodFrames(frame).mocapJoints(i).Y = ypoints(i);
            goodFrames(frame).mocapJoints(i).Z = zpoints(i);

            BigX = [xpoints(i); ypoints(i); zpoints(i); 1];
            littlex = vue2.Kmat * vue2.Pmat * BigX;
            vec2 = homo2cart2d(littlex);
            goodFrames(frame).left.x(i) = vec2(1);
            goodFrames(frame).left.y(i) = vec2(2);
            


        % convert to 2d point on right image


            BigX = [xpoints(i); ypoints(i); zpoints(i); 1];
            littlex = vue4.Kmat * vue4.Pmat * BigX;
            vec2 = homo2cart2d(littlex);
            goodFrames(frame).right.x(i) = vec2(1);
            goodFrames(frame).right.y(i) = vec2(2);

        end

    end
    fprintf("Projection from 3d points X to 2d points x for each frame on both cameras Complete Successfully\n");
    fprintf("delta Time = %s\n", string(datetime('now') -starttime));
    

%% 2.4 Triangulation back into a set of 3D scene points
    starttime = datetime('now');
    p1 = vue2.myp(1,:);
    p2 = vue2.myp(2,:);
    p3 = vue2.myp(3,:);
    p1prime = vue4.myp(1,:);
    p2prime = vue4.myp(2,:);
    p3prime = vue4.myp(3,:);

    M = 8;
    parfor (frame = 1:length(goodFrames), M)
        for i = 1:1:12 % for each point

            x = goodFrames(frame).left.x(i);
            y = goodFrames(frame).left.y(i);
            xprime = goodFrames(frame).right.x(i);
            yprime = goodFrames(frame).right.y(i);
            
            A = [y*p3 - p2;
                p1 - x*p3;
                yprime*p3prime - p2prime
                p1prime - xprime*p3prime];

            % solve AX = 0 using S V D

            [u, d] = eigs(A' * A);
            
            uu = u(:,4);
            uu = uu/uu(4);
            uu = homo2cart3d(uu);

            goodFrames(frame).calculatedjoints(i).x = uu(1);
            goodFrames(frame).calculatedjoints(i).y = uu(2);
            goodFrames(frame).calculatedjoints(i).z = uu(3);

        end
    end
    fprintf("Trianglation from 2d points x to 3d points X for each frame Complete Successfully\n")
    fprintf("delta Time = %s\n", string(datetime('now') -starttime));

%% 2.5 Measure error betweeen trianglulated and original 3D points
    starttime = datetime('now');
    ittorator = 1;
    for frame = 1:length(goodFrames)
        mocapFnum = goodFrames(frame).frame; % for each frame
        dxSum = 0;
        
        for i = 1:1:12 % for each point

            % Compute the L^2 (Euclidean) distance.
            % ((X-X')^2 + (Y-Y')^2 + (Z-Z')^2)
            actualJoint = [goodFrames(frame).mocapJoints(i).X;
                goodFrames(frame).mocapJoints(i).Y;
                goodFrames(frame).mocapJoints(i).Z];
            calculatedJoint = [goodFrames(frame).calculatedjoints(i).x ;
                goodFrames(frame).calculatedjoints(i).y;
                goodFrames(frame).calculatedjoints(i).z];
            dx = Euclideandistance(actualJoint, calculatedJoint);
            dxSum = dxSum + dx;

            % 1 joint all framse error 
             jointerror(i).list(frame) = dx;
            % all the frames all the joints
            error(ittorator) = dx;
            ittorator = ittorator + 1;
        end
        % every frame all the joints
        frameerror(frame) = dxSum;
    end
    fprintf("The error between pairs of 3d points on each frame Complete Successfully\n")
    fprintf("delta Time = %s\n", string(datetime('now') -starttime));
%% Computing epipolar lines between the two views. 
    %profile viewer
    for frame = 1:60:60%length(goodFrames)
        % create lines
        %{
        f = figure;
        f.Position = [100 400 1800 500];
        hold on;
        for i = 1:11
        X(i) = goodFrames(frame).mocapJoints(i).X;
        Y(i) = goodFrames(frame).mocapJoints(i).Y;
        Z(i) = goodFrames(frame).mocapJoints(i).Z;
        end
        plot3(X, Y, Z,'-x', 'MarkerSize', 5, 'LineWidth', 6, 'Color','blue');
        for i = 1:11
        X(i) = goodFrames(frame).calculatedjoints(i).x;
        Y(i) = goodFrames(frame).calculatedjoints(i).y;
        Z(i) = goodFrames(frame).calculatedjoints(i).z;
        end
        plot3(X, Y, Z,'-x', 'MarkerSize', 5, 'LineWidth', 1, 'Color','red');
        xlim([0, 1920])
        ylim([0, 1080])


       
        vue2video.CurrentTime = (mocapFnum-1)*(50/100)/vue2video.FrameRate;
        vid2Frame = readFrame(vue2video);
        f = figure;
        f.Position = [100 400 1800 500];
        f.Visible = true;
        subplot(1,2,2)
        imshow(vid2Frame);
        axis on
        hold on;
        plot(goodFrames(frame).left.x, goodFrames(frame).left.y, '-x', 'MarkerSize', 5, 'LineWidth', 1);
        xlim([0, 1920])
        ylim([0, 1080])
        vue4video.CurrentTime = (mocapFnum-1)*(50/100)/vue4video.FrameRate;
        vid4Frame = readFrame(vue4video);
        subplot(1,2,1)
        imshow(vid4Frame);
        axis on
        hold on;
        % Plot cross at row 100, column 50
        plot(goodFrames(frame).right.x, goodFrames(frame).right.y, '-x', 'MarkerSize', 5, 'LineWidth', 1);
        xlim([0, 1920])
        ylim([0, 1080])
        %}
        %{
        % Calucation the Fundamental matrix F
        % This is commeted out after talking to the TA and not having
        % enought time to debug, Also it is slower to calcualte the F
        % matrix for every frame compared to projecting the camera positon
        % onto each frame and connecting the dots.

        A = eye(9);
        for i = 1:1:9
            x = goodFrames(1).right.x(i);
            y = goodFrames(1).right.y(i);
            xprime = goodFrames(1).left.x(i);
            yprime = goodFrames(1).left.y(i);
            A(i,:) = [x*xprime x*yprime x y*xprime y*yprime y xprime yprime 1];
        end
        % use S V D to find f1 - f9
        [U,D,V] = svd(A);
        F = reshape(V(:,9),3,3)';
        [U,D,V] = svd(F);
        F = U*diag([D(1,1) D(2,2) 0])*V';
        
        % change the vector matrix into a 3x4 matrix for F
        %F = [uu(1) uu(2) uu(3);
        %    uu(4) uu(5) uu(6);
        %    uu(7) uu(8) uu(9)];
        % Calculate the epiplor lines for the left image. 
        % l = F^T * x' (line on the left image = fundamental matrix transpose
        % Multiplied by a point in the right image
        x = [goodFrames(1).right.x(i);
            goodFrames(1).right.y(i);
            1];
        l = F * x
        % convert l into a line eqution
        %line  = x*l(1) + y*l(2) + 1*(3) = 0
        line = @(x) (x*l(1) + l(3))/l(2);
        y1 = line(0);
        y2 = line(1920);
        %}
        mocapFnum = goodFrames(frame).frame; % for each frame
        %camera position
        BigX = [vue4.position(1); vue4.position(2); vue4.position(3); 1];
        littlex = vue2.Kmat * vue2.Pmat * BigX;
        vec2 = homo2cart2d(littlex);
        goodFrames(frame).left.camerax = vec2(1);
        goodFrames(frame).left.cameray = vec2(2);
    
        BigX = [vue2.position(1); vue2.position(2); vue2.position(3); 1];
        littlex = vue4.Kmat * vue4.Pmat * BigX;
        vec2 = homo2cart2d(littlex);
        goodFrames(frame).right.camerax = vec2(1);
        goodFrames(frame).right.cameray = vec2(2);
    
        
        % create lines
        line1 = [[goodFrames(frame).left.x(1) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(1) goodFrames(1).left.cameray]]; 
        line2 = [[goodFrames(frame).left.x(2) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(2) goodFrames(1).left.cameray]];
        line3 = [[goodFrames(frame).left.x(3) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(3) goodFrames(1).left.cameray]];
        line4 = [[goodFrames(frame).left.x(4) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(4) goodFrames(1).left.cameray]];
        line5 = [[goodFrames(frame).left.x(5) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(5) goodFrames(1).left.cameray]];
        line6 = [[goodFrames(frame).left.x(6) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(6) goodFrames(1).left.cameray]];
        line7 = [[goodFrames(frame).left.x(7) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(7) goodFrames(1).left.cameray]];
        line8 = [[goodFrames(frame).left.x(8) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(8) goodFrames(1).left.cameray]];
        line9 = [[goodFrames(frame).left.x(9) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(9) goodFrames(1).left.cameray]];
        line10 = [[goodFrames(frame).left.x(10) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(10) goodFrames(1).left.cameray]];
        line11 = [[goodFrames(frame).left.x(11) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(11) goodFrames(1).left.cameray]];
        line12 = [[goodFrames(frame).left.x(12) goodFrames(1).left.camerax]; [goodFrames(frame).left.y(12) goodFrames(1).left.cameray]];
        [line1(1, 1), line1(2, 1)] = increaselinesegment(line1, 1500);
        [line2(1, 1), line2(2, 1)] = increaselinesegment(line2, 1500);
        [line3(1, 1), line3(2, 1)] = increaselinesegment(line3, 1500);
        [line4(1, 1), line4(2, 1)] = increaselinesegment(line4, 1500);
        [line5(1, 1), line5(2, 1)] = increaselinesegment(line5, 1500);
        [line6(1, 1), line6(2, 1)] = increaselinesegment(line6, 1500);
        [line7(1, 1), line7(2, 1)] = increaselinesegment(line7, 1500);
        [line8(1, 1), line8(2, 1)] = increaselinesegment(line8, 1500);
        [line9(1, 1), line9(2, 1)] = increaselinesegment(line9, 1500);
        [line10(1, 1), line10(2, 1)] = increaselinesegment(line10, 1500);
        [line11(1, 1), line11(2, 1)] = increaselinesegment(line11, 1500);
        [line12(1, 1), line12(2, 1)] = increaselinesegment(line12, 1500);
    
        
        vue2video.CurrentTime = (mocapFnum-1)*(50/100)/vue2video.FrameRate;
        vid2Frame = readFrame(vue2video);
        f = figure;
        f.Position = [100 400 1800 500];
        f.Visible = true;
        subplot(1,2,2)
        imshow(vid2Frame);
        axis on
        hold on;
        % Plot cross at row 100, column 50
        plot([line1(1, 1) line1(1, 2)], [line1(2, 1) line1(2, 2)],  'cyan',...
            [line2(1, 1) line2(1, 2)], [line2(2, 1) line2(2, 2)],  'blue', ...
            [line3(1, 1) line3(1, 2)], [line3(2, 1) line3(2, 2)],  ...
            [line4(1, 1) line4(1, 2)], [line4(2, 1) line4(2, 2)],  ...
            [line5(1, 1) line5(1, 2)], [line5(2, 1) line5(2, 2)],  ...
            [line6(1, 1) line6(1, 2)], [line6(2, 1) line6(2, 2)],  ...
            [line7(1, 1) line7(1, 2)], [line7(2, 1) line7(2, 2)],  ...
            [line8(1, 1) line8(1, 2)], [line8(2, 1) line8(2, 2)],  ...
            [line9(1, 1) line9(1, 2)], [line9(2, 1) line9(2, 2)],  ...
            [line10(1, 1) line10(1, 2)], [line10(2, 1) line10(2, 2)],  ...
            [line11(1, 1) line11(1, 2)], [line11(2, 1) line11(2, 2)],  ...
            [line12(1, 1) line12(1, 2)], [line12(2, 1) line12(2, 2)], 'yellow', ...
            goodFrames(frame).left.x, goodFrames(frame).left.y, 'o', 'MarkerSize', 2, 'LineWidth', 1);
        
        vue4video.CurrentTime = (mocapFnum-1)*(50/100)/vue4video.FrameRate;
        vid4Frame = readFrame(vue4video);
    
        subplot(1,2,1)
    
        line1 = [[goodFrames(frame).right.x(1) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(1) goodFrames(1).right.cameray]];
        line2 = [[goodFrames(frame).right.x(2) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(2) goodFrames(1).right.cameray]];
        line3 = [[goodFrames(frame).right.x(3) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(3) goodFrames(1).right.cameray]];
        line4 = [[goodFrames(frame).right.x(4) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(4) goodFrames(1).right.cameray]];
        line5 = [[goodFrames(frame).right.x(5) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(5) goodFrames(1).right.cameray]];
        line6 = [[goodFrames(frame).right.x(6) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(6) goodFrames(1).right.cameray]];
        line7 = [[goodFrames(frame).right.x(7) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(7) goodFrames(1).right.cameray]];
        line8 = [[goodFrames(frame).right.x(8) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(8) goodFrames(1).right.cameray]];
        line9 = [[goodFrames(frame).right.x(9) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(9) goodFrames(1).right.cameray]];
        line10 = [[goodFrames(frame).right.x(10) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(10) goodFrames(1).right.cameray]];
        line11 = [[goodFrames(frame).right.x(11) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(11) goodFrames(1).right.cameray]];
        line12 = [[goodFrames(frame).right.x(12) goodFrames(1).right.camerax]; [goodFrames(frame).right.y(12) goodFrames(1).right.cameray]];
        [line1(1, 1), line1(2, 1)] = increaselinesegment2(line1, 1500);
        [line2(1, 1), line2(2, 1)] = increaselinesegment2(line2, 1500);
        [line3(1, 1), line3(2, 1)] = increaselinesegment2(line3, 1500);
        [line4(1, 1), line4(2, 1)] = increaselinesegment2(line4, 1500);
        [line5(1, 1), line5(2, 1)] = increaselinesegment2(line5, 1500);
        [line6(1, 1), line6(2, 1)] = increaselinesegment2(line6, 1500);
        [line7(1, 1), line7(2, 1)] = increaselinesegment2(line7, 1500);
        [line8(1, 1), line8(2, 1)] = increaselinesegment2(line8, 1500);
        [line9(1, 1), line9(2, 1)] = increaselinesegment2(line9, 1500);
        [line10(1, 1), line10(2, 1)] = increaselinesegment2(line10, 1500);
        [line11(1, 1), line11(2, 1)] = increaselinesegment2(line11, 1500);
        [line12(1, 1), line12(2, 1)] = increaselinesegment2(line12, 1500);
    
        imshow(vid4Frame);
        axis on
        hold on;
        % Plot cross at row 100, column 50
        plot([line1(1, 1) line1(1, 2)], [line1(2, 1) line1(2, 2)], 'cyan', ...
            [line2(1, 1) line2(1, 2)], [line2(2, 1) line2(2, 2)],  'blue', ...
            [line3(1, 1) line3(1, 2)], [line3(2, 1) line3(2, 2)],  ...
            [line4(1, 1) line4(1, 2)], [line4(2, 1) line4(2, 2)],  ...
            [line5(1, 1) line5(1, 2)], [line5(2, 1) line5(2, 2)],  ...
            [line6(1, 1) line6(1, 2)], [line6(2, 1) line6(2, 2)], ...
            [line7(1, 1) line7(1, 2)], [line7(2, 1) line7(2, 2)], ...
            [line8(1, 1) line8(1, 2)], [line8(2, 1) line8(2, 2)], ...
            [line9(1, 1) line9(1, 2)], [line9(2, 1) line9(2, 2)],  ...
            [line10(1, 1) line10(1, 2)], [line10(2, 1) line10(2, 2)],  ...
            [line11(1, 1) line11(1, 2)], [line11(2, 1) line11(2, 2)],  ...
            [line12(1, 1) line12(1, 2)], [line12(2, 1) line12(2, 2)], 'yellow', ...
            goodFrames(frame).right.x, goodFrames(frame).right.y, 'o', 'MarkerSize', 2, 'LineWidth', 1);
        
        %path = strcat('Output/epipolarLines', string(frame), '.png');
        %saveas(gcf,path)
        %close(f)
    end


%% Evaluation
fprintf("--- Error For each Frame ---\n")
fprintf("Mean %s\n", string(mean(frameerror,'all')));
fprintf("Standard deviation %s\n", string(std(frameerror)));
fprintf("Minimum %s\n", string(min(frameerror)));
fprintf("Median %s\n", string(median(frameerror)));
fprintf("Max %s\n", string(max(frameerror)));

fprintf("--- Total Error ---\n")
fprintf("Mean %s\n", string(mean(error,'all')));
fprintf("Standard deviation %s\n", string(std(error)));
fprintf("Minimum %s\n", string(min(error)));
fprintf("Median %s\n", string(median(error)));
fprintf("Max %s\n", string(max(error)));
fprintf("Sum %s\n", string(sum(error)));

fprintf("--- Joint Error ---\n")
fprintf("joint 1,\t Mean,\t\t\t STD,\t\t\t Min,\t\t\t Midian,\t\t\t Max \n")
for i = 1:12
    fprintf("joint %s,\t %s,\t %s,\t %s,\t %s,\t %s \n", string(i), string(mean(jointerror(i).list)), ...
        string(std(jointerror(i).list)), string(min(jointerror(i).list)), ...
        string(median(jointerror(i).list)), string(max(jointerror(i).list)) )
end
%plot(frameerror)

%% Fucntions 

function vec2 = homo2cart2d(vec3)
    x = vec3(1)/vec3(3);
    y = vec3(2)/vec3(3);
    vec2 = [x; y];
end

function delta = Euclideandistance(vec3point1, vec3point2)
    delta = sqrt((vec3point1(1)-vec3point2(1))^2 + (vec3point1(2)-vec3point2(2))^2 + (vec3point1(3)-vec3point2(3))^2);
end

function vec3 = homo2cart3d(vec4)
    x = vec4(1)/vec4(4);
    y = vec4(2)/vec4(4);
    z = vec4(3)/vec4(4);
    vec3 = [x; y; z];
end

function p = myp(vue)
    % K = Intrinsics
    K = padarray(vue.Kmat,[0 1],0,'post');

    % R = rotation 3x3 padded 
    R = padarray(vue.Rmat,[1 1],0,'post');
    R(4, 4) = 1;

    % IC = the position matix
    IC = padarray(vue.position'*-1,[0 3],0,'pre');
    IC(1, 1) = 1;
    IC(2, 2) = 1;
    IC(3, 3) = 1;
    IC = padarray(IC,[1 0],0,'post');
    IC(4, 4) = 1;

    % P = KR{I|-C}
    p = K * R * IC;
end

function [newx, newy] = increaselinesegment(line, amount)
    slope = abs((line(2,1)-line(2,2))/(line(1,1)-line(1,2)));
    newx = line(1, 1) + amount;
    newy = line(2,1) + amount*slope;
end

function [newx, newy] = increaselinesegment2(line, amount)
    slope = abs((line(2,1)-line(2,2))/(line(1,1)-line(1,2)));
    newx = line(1, 1) - amount;
    newy = line(2,1) + amount*slope;
end

    
    