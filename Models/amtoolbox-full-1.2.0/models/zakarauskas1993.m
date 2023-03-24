function [resp,hits,err,averr]=zakarauskas1993(ir,pol,stim,q10,do,bal,fstart,fend,fs,space)
% ZAKARAUSKAS1993 Localization model according to Zakarauskas et al. (1993)
%
%   Usage:
%     [resp,hits,err,averr]=zakarauskas1993(ir,pol)
%     [resp,hits,err,averr]=zakarauskas1993(ir,pol,stim)
%     [resp,hits,err,averr]=zakarauskas1993(ir,pol,stim,q10,do,bal,fstart,fend,fs,space)
%
%   Input parameters:
%     ir:       impulse responses of DFTs for all positions and
%               both ears (sorted)
%     pol:      polar angles of ir
%     stim:     magnitude spectrum of stimulus
%     q10:      relativ bandwidth of filter bands;  default: 10
%     do:       differential order;                 default: 2
%     bal:      balance of left to right channel;   default: 1
%     fstart:   start frequency; minimum: 0,5kHz;   default: 2kHz
%     fend:     end frequency;                      default: 16kHz
%     fs:       sampling frequency;                 default: 48kHz
%     space:    spacing factor of filter bands;     default: 1.03
%
%   Output parameters:
%     resp:     summed differences for response angles 
%               with respect to target positions
%     hits:     number of correct estimations
%     err:      error of estimation in degree
%     averr:    average absolute error of estimation in degree
%
%   localization model
%
%   Url: http://amtoolbox.org/amt-1.2.0/doc/models/zakarauskas1993.php

% Copyright (C) 2009-2022 Piotr Majdak, Clara Hollomey, and the AMT team.
% This file is part of Auditory Modeling Toolbox (AMT) version 1.2.0
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   #StatusDoc: Good
%   #StatusCode: Submitted
%   #Verification: Unknown
%   #Requirements: M-Signal
%   #Author: Robert Baumgartner, OEAW Acoustical Research Institute
% latest update: 2010-08-03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default settings
if ~exist('q10','var')
    q10=10;
end
if ~exist('do','var')
    do=2;
end
if ~exist('bal','var')
    bal=1;
end
if ~exist('fstart','var')
    fstart=2000;
end
if ~exist('fend','var')
    fend=16000;
end
if ~exist('fs','var')
    fs=48000;
end
if ~exist('space','var')
    space=1.03;
end

% model calculations 
if exist('stim','var')
    ir1=zeros(size(ir,1)+size(stim,1)-1,size(ir,2),size(ir,3));
    for ch=1:size(ir,3)
        for ind=1:size(ir,2)
            ir1(:,ind,ch) = conv(ir(:,ind,ch),stim);
        end
    end
end

[x,y]=zakarauskas1993_butterfb(ir1,ir,q10,fstart,fend,fs,space); % filter bank
resp=zeros(size(ir1,2)); 
err=zeros(size(ir1,2),1);
for ii=1:size(x,2) % analysis for each target position
    resp(:,ii)=zakarauskas1993_compare( x(:,ii,:),y,do,bal );
    [m,idx]=min(resp(:,ii));
    err(ii)=pol(idx)-pol(ii);
end
hits=sum(err==0);
averr=mean(abs(err));
end
