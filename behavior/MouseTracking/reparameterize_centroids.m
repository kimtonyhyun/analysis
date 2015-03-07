function [intermediate_centroids,param_centroids] = reparameterize_centroids(centroids)
    % Re-parametrize the spatial trace such that:
    %   south is (-*,0)
    %   north is (+*,0)
    %   east is (0,+*)
    %   west is (0,-*)
    %
    %   spatial_trace : [# of frames] x 2 array
    % 
    %   param_trace : parametrized spatial trace with the same size as
    %   spatial_trace
    %
    % Hakan Inan (Mar 15) 
    %
    T = size(centroids,1);
    if size(centroids,2)~=2,
        error('input must be an Nx2 array');
    end
    %Sometimes the initial location is not detected when obtaining the spatial
    %trace from the behavior camera, resulting in a bunch of (0,0)'s in the
    %beginning. We deal with it by replacing them with the first detected
    %location
    if sum(centroids(1,:) == [0,0])==2
        stop = 0;
        acc = 0;
        while stop==0
            acc = acc+1;
            if sum(centroids(acc,:) == [0,0])~=2 % index acc has a detected location
                stop = 1;
            end
        end
        %replace (0,0)'s with the first detected location
        centroids(1:acc-1,:) = centroids(acc,:);
    end

    % Use segmented least squares fitting with one knot. The knot is
    % found by trying every possible one and selecting the one with lowest
    % reconstruction error.
    coord_x = centroids(:,1);
    coord_y = -centroids(:,2);
    errorVec = zeros(T,1);
    for k = 1:T
        dum = zeros(T,1);
        dum(1:k) = 1;
        dum_prime = 1-dum;
        A = [dum.*coord_x+dum_prime.*coord_y,dum,dum_prime];
        b = -dum_prime.*coord_x+dum.*coord_y;
        u = A\b;
        errorVec(k) = norm(A*u-b);
    end
    [~,idx] = min(errorVec);
    
    % Do the fit again with the knot yielding the best error
    dum = zeros(T,1);
    dum(1:idx) = 1;
    dum_prime = 1-dum;
    A = [dum.*coord_x+dum_prime.*coord_y,dum,dum_prime];
    b = -dum_prime.*coord_x+dum.*coord_y;
    u = A\b;
    slope = u(1);
    intercept1 = u(2);
    intercept2 = u(3)/-slope;
    
    % Reconstruct y trace from x trace
    y_reconst = dum.*(slope*coord_x+intercept1)+...
        dum_prime.*(-1/slope *coord_x + intercept2);    
    intermediate_centroids = [coord_x,y_reconstruct];
    
    % Find x,y points where two lines intersect, that point is the knot
    knot = [-slope,1;1/slope,1] \ [intercept1;intercept2];
    knot_x = knot(1);
    knot_y = knot(2);
    
    % Translate
    y_shifted = y_reconst - knot_y;
    x_shifted = coord_x - knot_x;
    
    % Rotate
    angle = atan(slope);
    angToShift = angle;
    rotationMat = [cos(angToShift),-sin(angToShift);sin(angToShift),cos(angToShift)];
    param_centroids = [x_shifted,y_shifted]*rotationMat;

    %round anything in [-1,1] to zero
    param_centroids(param_centroids<1 & param_centroids>-1) = 0;

    % Median filtering to get rid of noise in the trajectory(optional, I guess)
    % param_centroids(:,1) = medfilt1(param_centroids(:,1),1);
    % param_centroids(:,2) = medfilt1(param_centroids(:,2),1);
