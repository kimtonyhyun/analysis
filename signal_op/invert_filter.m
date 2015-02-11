function g = invert_filter(h,g_len,impulse_loc)
%find an inverse of h for deconvolution:
%Idea is to work in the time domain and try to find a vector g that yields
%an 'almost' impulse when convolved with h. We fit Ag=y with least squares
%penalty.
h_len = length(h);
A = zeros(g_len+h_len,g_len);%a bunch of shifted h's
y = zeros(g_len+h_len,1); %output vector
%we want to force g to produce zero when convolved with noise.
%below noise matrix is to be appended to A
noise_length = h_len*4;%event_length*2 is arbitrary(no special reason)
noise = wgn(noise_length,g_len,0)*0.4; 
%append the noise terms and corresponding zeros in y
y = [y;zeros(noise_length,1)];
A = [A;noise];
dum = [zeros(1,g_len),zeros(1,h_len),fliplr(h)];
for k = 1:g_len+h_len
    dum = circshift(dum',1)'; %matlab only shifts column vectors..
    row_A = dum(1:g_len);
    A(k,:) = row_A;
    if k == floor(impulse_loc);%g_len/2+impulse_loc); %only case where we want y=1
        y(k) = 1;
    end
end
%do least squares fit to get g:(The form is so that regularization is possible)
g = (A'*A+eye(g_len)*0)\A'*y;
g=flipud(g');
%another approach : frequency domain filtering(Wiener-ish deconvolution)
% H = fft(h);
% G = 1./H .* (abs(H).^2 ./ (abs(H).^2+10));
% g = ifft(G);