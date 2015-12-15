function [ pseudoim ] = imcart2pseudopolar( im, rhoSampling, thetaSampling )
%FUNCTION [ pseudoim ] = imcart2pseudopolar( im, rhoSampling, thetaSampling )
%   Robert Cooper
%
% This function takes in an image and converts it to pseudopolar, where
% Rho (the radius) is as long as half of the shortest side of the image,
% and Theta is sampled every other degree (by default)
%
% Change rhoUpsampling and thetaSampling to increase/decrease the sampling of the image.
% If you wish to upsample, then lower the *Sampling.
% If you wish to downsample, then raise the *Sampling.

if ~exist('rhoSampling','var') || isempty(rhoSampling)
     rhoSampling =1;
end

if ~exist('thetaSampling','var') || isempty(thetaSampling)
     thetaSampling =1;
end

im = double(im);

[X, Y]= meshgrid( 1:size(im,2), 1:size(im,1) );

% X = X - (size(im,2)/2);
% Y = Y - (size(im,1)/2);

rho = 1:rhoSampling: floor(min(size(im))/2)-1;
theta_step = thetaSampling*pi/360;
theta = 0: theta_step: 2*pi-theta_step;

[R,T] = meshgrid(rho,theta);

[Rx, Ty] = pol2cart(T,R);

Rx = Rx + ceil(size(im,2)/2);
Ty = Ty + ceil(size(im,2)/2);

pseudoim = interp2(X,Y,im,Rx,Ty);

pseudoim(isnan(pseudoim)) = 0;
% imagesc( pseudoim ); colormap gray; axis image;
end

