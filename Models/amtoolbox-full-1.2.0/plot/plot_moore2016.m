function out = plot_moore2016(ShortTermLoudness, LongTermLoudness,varargin)
%PLOT_MOORE2016 plots from Moore et al. 2016
%
%   Usage:
%     out = plot_moore2016(ShortTermLoudness, LongTermLoudness)
%
%   Input parameters:
%     ShortTermLoudness   : as calculated by Moore2016 [sone]
%     LongTermLoudness    : as calculated by Moore2016  [sone]
%
%   Output parameters:
%     out   : interpolated loudness [phon]
%
%   This function plots the long term- and short term loudness [sone] and
%   converts it to phon.
%
%   Url: http://amtoolbox.org/amt-1.2.0/doc/plot/plot_moore2016.php

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

%filenameSound = 'testSound';
%filenameFilter = 'default';
%dBMax = 100; 

figure;
plot( 0:(length(ShortTermLoudness)-1), ShortTermLoudness, 'b-' );
hold on;
plot( 0:(length(LongTermLoudness)-1), LongTermLoudness, 'r-' );
xlabel( 'time [ms]');
ylabel( 'Loudness [sone]');
legend( 'Short-term loudness','Long-term loudness' );

% strOutputFilename = [filenameSound ' ' num2str( dBMax ) ' dB calibration level TVL 2018.txt'];
% fid = fopen(strOutputFilename,'wt+');
% fprintf(fid,[filenameSound '\n\n']);
% fprintf(fid,['Calibration level:      ' num2str(dBMax) ' dB SPL (RMS level of a full-scale sinusoid)\n'] );
% fprintf(fid,['Filename of FIR filter: ' filenameFilter '\n\n']);
% fprintf(fid,['Maximum of long-term loudness:  ' num2str(max(LongTermLoudness),'%9.2f') ' sone\n'] );
% fprintf(fid,['                                ' num2str(max(Sone2Phon(LongTermLoudness)),'%9.2f') ' phon\n'] );
% fprintf(fid,['Maximum of short-term loudness: ' num2str(max(ShortTermLoudness),'%9.2f') ' sone\n'] );
% fprintf(fid,['                                ' num2str(max(Sone2Phon(ShortTermLoudness)),'%9.2f') ' phon\n\n'] );
% fprintf(fid,'Loudness over time\n');
% fprintf(fid,'1st column: time in milliseconds\n');
% fprintf(fid,'2nd column: short-term loudness in sone\n');
% fprintf(fid,'3rd column: short-term loudness level in phon\n');
% fprintf(fid,'4th column: long-term loudness in sone\n');
% fprintf(fid,'5th column: long-term loudness level in phon\n\n');
% fprintf(fid,'   time   short-t. loudness    long-t. loudness\n');
% fprintf(fid,'     ms      sone      phon      sone      phon\n');
% for i=1:length( LongTermLoudness )
%     fprintf(fid,'%7.0f %9.2f %9.1f %9.2f %9.1f\n',i-1, ShortTermLoudness(i), Sone2Phon(ShortTermLoudness(i)), LongTermLoudness(i), Sone2Phon(LongTermLoudness(i)) );
% end
% fprintf(fid,'max     %9.2f %9.1f %9.2f %9.1f\n', max(ShortTermLoudness), max(Sone2Phon(ShortTermLoudness)), max(LongTermLoudness), max(Sone2Phon(LongTermLoudness)) );
% fclose(fid);
% out = 1;

function out = Sone2Phon( in )

SonePhone = [ 0.000538 0
0.000566 0.1
0.000596 0.2
0.000627 0.3
0.000660 0.4
0.000694 0.5
0.000730 0.6
0.000768 0.7
0.000808 0.8
0.000849 0.9
0.000893 1
0.000939 1.1
0.000987 1.2
0.001037 1.3
0.001089 1.4
0.001144 1.5
0.001202 1.6
0.001262 1.7
0.001326 1.8
0.001392 1.9
0.001461 2
0.001533 2.1
0.001609 2.2
0.001688 2.3
0.001770 2.4
0.001857 2.5
0.001947 2.6
0.002043 2.7
0.002142 2.8
0.002245 2.9
0.002352 3
0.002464 3.1
0.002581 3.2
0.002703 3.3
0.002830 3.4
0.002963 3.5
0.003101 3.6
0.003245 3.7
0.003395 3.8
0.003551 3.9
0.003714 4
0.003883 4.1
0.004060 4.2
0.004243 4.3
0.004434 4.4
0.004632 4.5
0.004838 4.6
0.005052 4.7
0.005275 4.8
0.005506 4.9
0.005746 5
0.005995 5.1
0.006254 5.2
0.006522 5.3
0.006801 5.4
0.007089 5.5
0.007389 5.6
0.007699 5.7
0.008020 5.8
0.008353 5.9
0.008698 6
0.009040 6.1
0.009362 6.2
0.009693 6.3
0.010036 6.4
0.010390 6.5
0.010756 6.6
0.011090 6.7
0.011430 6.8
0.011780 6.9
0.012111 7
0.012434 7.1
0.012765 7.2
0.013105 7.3
0.013454 7.4
0.013813 7.5
0.014181 7.6
0.014559 7.7
0.014947 7.8
0.015346 7.9
0.015755 8
0.016176 8.1
0.016607 8.2
0.017049 8.3
0.017504 8.4
0.017970 8.5
0.018448 8.6
0.018905 8.7
0.019358 8.8
0.019822 8.9
0.020289 9
0.020727 9.1
0.021175 9.2
0.021631 9.3
0.022098 9.4
0.022574 9.5
0.023060 9.6
0.023556 9.7
0.024062 9.8
0.024579 9.9
0.025107 10
0.025645 10.1
0.026195 10.2
0.026756 10.3
0.027329 10.4
0.027913 10.5
0.028510 10.6
0.029118 10.7
0.029740 10.8
0.030373 10.9
0.031029 11
0.031689 11.1
0.032363 11.2
0.033050 11.3
0.033751 11.4
0.034466 11.5
0.035196 11.6
0.035941 11.7
0.036700 11.8
0.037464 11.9
0.038158 12
0.038862 12.1
0.039577 12.2
0.040303 12.3
0.041041 12.4
0.041790 12.5
0.042552 12.6
0.043326 12.7
0.044112 12.8
0.044910 12.9
0.045721 13
0.046545 13.1
0.047382 13.2
0.048233 13.3
0.049096 13.4
0.049973 13.5
0.050864 13.6
0.051769 13.7
0.052688 13.8
0.053621 13.9
0.054570 14
0.055533 14.1
0.056511 14.2
0.057504 14.3
0.058512 14.4
0.059536 14.5
0.060576 14.6
0.061632 14.7
0.062705 14.8
0.063793 14.9
0.064899 15
0.066022 15.1
0.067162 15.2
0.068320 15.3
0.069453 15.4
0.070597 15.5
0.071757 15.6
0.072934 15.7
0.074127 15.8
0.075308 15.9
0.076482 16
0.077670 16.1
0.078873 16.2
0.080090 16.3
0.081323 16.4
0.082570 16.5
0.083833 16.6
0.085116 16.7
0.086411 16.8
0.087721 16.9
0.089047 17
0.090389 17.1
0.091747 17.2
0.093122 17.3
0.094514 17.4
0.095923 17.5
0.097349 17.6
0.098792 17.7
0.100253 17.8
0.101731 17.9
0.103229 18
0.104743 18.1
0.106276 18.2
0.107828 18.3
0.109399 18.4
0.110988 18.5
0.112597 18.6
0.114225 18.7
0.115873 18.8
0.117541 18.9
0.119230 19
0.120938 19.1
0.122623 19.2
0.124325 19.3
0.126045 19.4
0.127785 19.5
0.129544 19.6
0.131323 19.7
0.133122 19.8
0.134912 19.9
0.136696 20
0.138498 20.1
0.140318 20.2
0.142156 20.3
0.144012 20.4
0.145887 20.5
0.147784 20.6
0.149696 20.7
0.151628 20.8
0.153578 20.9
0.155548 21
0.157538 21.1
0.159547 21.2
0.161576 21.3
0.163626 21.4
0.165696 21.5
0.167786 21.6
0.169903 21.7
0.172036 21.8
0.174190 21.9
0.176366 22
0.178563 22.1
0.180783 22.2
0.183024 22.3
0.185289 22.4
0.187576 22.5
0.189886 22.6
0.192219 22.7
0.194576 22.8
0.196957 22.9
0.199362 23
0.201791 23.1
0.204238 23.2
0.206669 23.3
0.209124 23.4
0.211602 23.5
0.214104 23.6
0.216609 23.7
0.219101 23.8
0.221615 23.9
0.224151 24
0.226708 24.1
0.229287 24.2
0.231888 24.3
0.234512 24.4
0.237158 24.5
0.239826 24.6
0.242518 24.7
0.245233 24.8
0.247972 24.9
0.250734 25
0.253520 25.1
0.256329 25.2
0.259169 25.3
0.262028 25.4
0.264912 25.5
0.267821 25.6
0.270756 25.7
0.273716 25.8
0.276702 25.9
0.279714 26
0.282753 26.1
0.285819 26.2
0.288911 26.3
0.292031 26.4
0.295179 26.5
0.298354 26.6
0.301558 26.7
0.304790 26.8
0.308052 26.9
0.311342 27
0.314663 27.1
0.318009 27.2
0.321331 27.3
0.324682 27.4
0.328056 27.5
0.331414 27.6
0.334799 27.7
0.338209 27.8
0.341646 27.9
0.345109 28
0.348600 28.1
0.352118 28.2
0.355663 28.3
0.359236 28.4
0.362838 28.5
0.366467 28.6
0.370125 28.7
0.373811 28.8
0.377527 28.9
0.381272 29
0.385049 29.1
0.388854 29.2
0.392689 29.3
0.396555 29.4
0.400452 29.5
0.404382 29.6
0.408341 29.7
0.412333 29.8
0.416357 29.9
0.420413 30
0.424503 30.1
0.428625 30.2
0.432782 30.3
0.436946 30.4
0.441110 30.5
0.445306 30.6
0.449533 30.7
0.453792 30.8
0.458083 30.9
0.462407 31
0.466763 31.1
0.471153 31.2
0.475576 31.3
0.480033 31.4
0.484524 31.5
0.489051 31.6
0.493612 31.7
0.498207 31.8
0.502838 31.9
0.507505 32
0.512167 32.1
0.516858 32.2
0.521583 32.3
0.526345 32.4
0.531141 32.5
0.535973 32.6
0.540841 32.7
0.545746 32.8
0.550687 32.9
0.555666 33
0.560682 33.1
0.565737 33.2
0.570793 33.3
0.575862 33.4
0.580967 33.5
0.586108 33.6
0.591285 33.7
0.596499 33.8
0.601750 33.9
0.607037 34
0.612362 34.1
0.617725 34.2
0.623126 34.3
0.628566 34.4
0.634045 34.5
0.639563 34.6
0.645120 34.7
0.650719 34.8
0.656357 34.9
0.662036 35
0.667756 35.1
0.673518 35.2
0.679322 35.3
0.685168 35.4
0.691057 35.5
0.696990 35.6
0.702966 35.7
0.708987 35.8
0.715052 35.9
0.721162 36
0.727299 36.1
0.733437 36.2
0.739619 36.3
0.745844 36.4
0.752112 36.5
0.758425 36.6
0.764782 36.7
0.771140 36.8
0.777537 36.9
0.783978 37
0.790463 37.1
0.796992 37.2
0.803566 37.3
0.810184 37.4
0.816848 37.5
0.823559 37.6
0.830315 37.7
0.837118 37.8
0.843968 37.9
0.850865 38
0.857811 38.1
0.864805 38.2
0.871847 38.3
0.878939 38.4
0.886081 38.5
0.893274 38.6
0.900517 38.7
0.907772 38.8
0.915052 38.9
0.922381 39
0.929759 39.1
0.937186 39.2
0.944662 39.3
0.952189 39.4
0.959766 39.5
0.967394 39.6
0.975074 39.7
0.982806 39.8
0.990590 39.9
0.998428 40
1.006318 40.1
1.014263 40.2
1.022262 40.3
1.030316 40.4
1.038425 40.5
1.046591 40.6
1.054813 40.7
1.063092 40.8
1.071429 40.9
1.079825 41
1.088279 41.1
1.096791 41.2
1.105299 41.3
1.113864 41.4
1.122486 41.5
1.131144 41.6
1.139833 41.7
1.148579 41.8
1.157382 41.9
1.166242 42
1.175160 42.1
1.184136 42.2
1.193172 42.3
1.202267 42.4
1.211422 42.5
1.220637 42.6
1.229914 42.7
1.239252 42.8
1.248652 42.9
1.258115 43
1.267642 43.1
1.277233 43.2
1.286888 43.3
1.296609 43.4
1.306395 43.5
1.316218 43.6
1.326069 43.7
1.335984 43.8
1.345962 43.9
1.356005 44
1.366114 44.1
1.376287 44.2
1.386528 44.3
1.396835 44.4
1.407210 44.5
1.417652 44.6
1.428164 44.7
1.438745 44.8
1.449396 44.9
1.460117 45
1.470910 45.1
1.481776 45.2
1.492713 45.3
1.503725 45.4
1.514810 45.5
1.525971 45.6
1.537207 45.7
1.548466 45.8
1.559783 45.9
1.571173 46
1.582638 46.1
1.594177 46.2
1.605791 46.3
1.617481 46.4
1.629248 46.5
1.641050 46.6
1.652924 46.7
1.664873 46.8
1.676900 46.9
1.689005 47
1.701189 47.1
1.713453 47.2
1.725796 47.3
1.738221 47.4
1.750727 47.5
1.763316 47.6
1.775988 47.7
1.788706 47.8
1.801474 47.9
1.814322 48
1.827252 48.1
1.840265 48.2
1.853360 48.3
1.866540 48.4
1.879804 48.5
1.893153 48.6
1.906588 48.7
1.920110 48.8
1.933720 48.9
1.947418 49
1.961206 49.1
1.975083 49.2
1.989052 49.3
2.003113 49.4
2.017266 49.5
2.031514 49.6
2.045832 49.7
2.060195 49.8
2.074650 49.9
2.089197 50
2.103836 50.1
2.118570 50.2
2.133399 50.3
2.148323 50.4
2.163344 50.5
2.178462 50.6
2.193679 50.7
2.208994 50.8
2.224411 50.9
2.239928 51
2.255548 51.1
2.271270 51.2
2.287098 51.3
2.303030 51.4
2.319006 51.5
2.335026 51.6
2.351147 51.7
2.367370 51.8
2.383696 51.9
2.400125 52
2.416658 52.1
2.433298 52.2
2.450043 52.3
2.466896 52.4
2.483858 52.5
2.500929 52.6
2.518111 52.7
2.535404 52.8
2.552810 52.9
2.570330 53
2.587965 53.1
2.605671 53.2
2.623459 53.3
2.641359 53.4
2.659373 53.5
2.677503 53.6
2.695747 53.7
2.714109 53.8
2.732589 53.9
2.751188 54
2.769908 54.1
2.788748 54.2
2.807712 54.3
2.826799 54.4
2.846012 54.5
2.865351 54.6
2.884817 54.7
2.904358 54.8
2.924000 54.9
2.943768 55
2.963663 55.1
2.983685 55.2
3.003837 55.3
3.024119 55.4
3.044533 55.5
3.065080 55.6
3.085761 55.7
3.106577 55.8
3.127531 55.9
3.148623 56
3.169855 56.1
3.191223 56.2
3.212657 56.3
3.234215 56.4
3.255879 56.5
3.277681 56.6
3.299622 56.7
3.321705 56.8
3.343930 56.9
3.366299 57
3.388813 57.1
3.411474 57.2
3.434283 57.3
3.457242 57.4
3.480352 57.5
3.503579 57.6
3.526910 57.7
3.550391 57.8
3.574022 57.9
3.597806 58
3.621744 58.1
3.645836 58.2
3.670086 58.3
3.694493 58.4
3.719061 58.5
3.743790 58.6
3.768682 58.7
3.793740 58.8
3.818907 58.9
3.844210 59
3.869677 59.1
3.895308 59.2
3.921106 59.3
3.947073 59.4
3.973209 59.5
3.999516 59.6
4.025997 59.7
4.052653 59.8
4.079486 59.9
4.106498 60
4.133633 60.1
4.160918 60.2
4.188380 60.3
4.216021 60.4
4.243843 60.5
4.271848 60.6
4.300038 60.7
4.328414 60.8
4.356978 60.9
4.385729 61
4.414629 61.1
4.443693 61.2
4.472890 61.3
4.502277 61.4
4.531854 61.5
4.561625 61.6
4.591591 61.7
4.621753 61.8
4.652115 61.9
4.682677 62
4.713443 62.1
4.744414 62.2
4.775538 62.3
4.806832 62.4
4.838332 62.5
4.870037 62.6
4.901952 62.7
4.934077 62.8
4.966414 62.9
4.998967 63
5.031737 63.1
5.064727 63.2
5.097897 63.3
5.131240 63.4
5.164802 63.5
5.198586 63.6
5.232593 63.7
5.266827 63.8
5.301288 63.9
5.335981 64
5.370907 64.1
5.406068 64.2
5.441391 64.3
5.476934 64.4
5.512714 64.5
5.548731 64.6
5.584989 64.7
5.621490 64.8
5.658236 64.9
5.695232 65
5.732478 65.1
5.769866 65.2
5.807475 65.3
5.845334 65.4
5.883445 65.5
5.921810 65.6
5.960432 65.7
5.999315 65.8
6.038460 65.9
6.077866 66
6.117452 66.1
6.157300 66.2
6.197415 66.3
6.237799 66.4
6.278455 66.5
6.319386 66.6
6.360594 66.7
6.402083 66.8
6.443782 66.9
6.485742 67
6.527983 67.1
6.570509 67.2
6.613323 67.3
6.656428 67.4
6.699828 67.5
6.743525 67.6
6.787435 67.7
6.831634 67.8
6.876133 67.9
6.920934 68
6.966040 68.1
7.011456 68.2
7.057185 68.3
7.103187 68.4
7.149452 68.5
7.196031 68.6
7.242928 68.7
7.290147 68.8
7.337691 68.9
7.385565 69
7.433741 69.1
7.482169 69.2
7.530895 69.3
7.579955 69.4
7.629352 69.5
7.679089 69.6
7.729171 69.7
7.779553 69.8
7.830233 69.9
7.881261 70
7.932641 70.1
7.984376 70.2
8.036471 70.3
8.088930 70.4
8.141658 70.5
8.194749 70.6
8.248207 70.7
8.302037 70.8
8.356244 70.9
8.410831 71
8.465724 71.1
8.520977 71.2
8.576616 71.3
8.632643 71.4
8.689065 71.5
8.745885 71.6
8.803023 71.7
8.860545 71.8
8.918471 71.9
8.976804 72
9.035552 72.1
9.094700 72.2
9.154185 72.3
9.214088 72.4
9.274415 72.5
9.335171 72.6
9.396361 72.7
9.457916 72.8
9.519880 72.9
9.582285 73
9.645135 73.1
9.708437 73.2
9.772135 73.3
9.836199 73.4
9.900720 73.5
9.965703 73.6
10.031154 73.7
10.097024 73.8
10.163318 73.9
10.230086 74
10.297336 74.1
10.365074 74.2
10.433229 74.3
10.501847 74.4
10.570961 74.5
10.640577 74.6
10.710689 74.7
10.781222 74.8
10.852266 74.9
10.923828 75
10.995916 75.1
11.068459 75.2
11.141503 75.3
11.215082 75.4
11.289204 75.5
11.363823 75.6
11.438938 75.7
11.514607 75.8
11.590837 75.9
11.667591 76
11.744852 76.1
11.822687 76.2
11.901103 76.3
11.980054 76.4
12.059540 76.5
12.139621 76.6
12.220305 76.7
12.301519 76.8
12.383314 76.9
12.465727 77
12.548752 77.1
12.632312 77.2
12.716504 77.3
12.801337 77.4
12.886747 77.5
12.972724 77.6
13.059357 77.7
13.146624 77.8
13.234481 77.9
13.323011 78
13.412218 78.1
13.502011 78.2
13.592495 78.3
13.683682 78.4
13.775474 78.5
13.867972 78.6
13.961192 78.7
14.055038 78.8
14.149611 78.9
14.244922 79
14.340880 79.1
14.437593 79.2
14.535043 79.3
14.633184 79.4
14.732105 79.5
14.831752 79.6
14.932147 79.7
15.033348 79.8
15.135251 79.9
15.237974 80
15.341474 80.1
15.445754 80.2
15.550875 80.3
15.656748 80.4
15.763487 80.5
15.871018 80.6
15.979404 80.7
16.088630 80.8
16.198701 80.9
16.309653 81
16.421449 81.1
16.534157 81.2
16.647719 81.3
16.762217 81.4
16.877587 81.5
16.993911 81.6
17.111131 81.7
17.229317 81.8
17.348401 81.9
17.468437 82
17.589443 82.1
17.711413 82.2
17.834402 82.3
17.958346 82.4
18.083329 82.5
18.209322 82.6
18.336338 82.7
18.464416 82.8
18.593519 82.9
18.723688 83
18.854944 83.1
18.987263 83.2
19.120666 83.3
19.255180 83.4
19.390808 83.5
19.527530 83.6
19.665376 83.7
19.804357 83.8
19.944477 83.9
20.085740 84
20.228148 84.1
20.371693 84.2
20.516391 84.3
20.662242 84.4
20.809247 84.5
20.957404 84.6
21.106710 84.7
21.257161 84.8
21.408742 84.9
21.561448 85
21.715274 85.1
21.870212 85.2
22.026348 85.3
22.183743 85.4
22.342405 85.5
22.502345 85.6
22.663571 85.7
22.826090 85.8
22.989914 85.9
23.155049 86
23.321504 86.1
23.489289 86.2
23.658411 86.3
23.828879 86.4
24.000697 86.5
24.173834 86.6
24.348341 86.7
24.524225 86.8
24.701494 86.9
24.880157 87
25.060222 87.1
25.241697 87.2
25.424591 87.3
25.608911 87.4
25.794667 87.5
25.981865 87.6
26.170516 87.7
26.360627 87.8
26.552207 87.9
26.745265 88
26.939809 88.1
27.135848 88.2
27.333392 88.3
27.532448 88.4
27.733027 88.5
27.935137 88.6
28.138788 88.7
28.343989 88.8
28.550751 88.9
28.759083 89
28.968994 89.1
29.180495 89.2
29.393596 89.3
29.608308 89.4
29.824641 89.5
30.042606 89.6
30.262213 89.7
30.483473 89.8
30.706399 89.9
30.931001 90
31.157290 90.1
31.385278 90.2
31.614978 90.3
31.846400 90.4
32.079558 90.5
32.314464 90.6
32.551129 90.7
32.789568 90.8
33.029792 90.9
33.271815 91
33.515649 91.1
33.761309 91.2
34.008808 91.3
34.258113 91.4
34.509281 91.5
34.762328 91.6
35.017266 91.7
35.274111 91.8
35.532878 91.9
35.793580 92
36.056233 92.1
36.320851 92.2
36.587450 92.3
36.856046 92.4
37.126653 92.5
37.399287 92.6
37.673965 92.7
37.950702 92.8
38.229514 92.9
38.510418 93
38.793431 93.1
39.078568 93.2
39.365848 93.3
39.655286 93.4
39.946901 93.5
40.240709 93.6
40.536729 93.7
40.834978 93.8
41.135473 93.9
41.438234 94
41.743278 94.1
42.050624 94.2
42.360291 94.3
42.672296 94.4
42.986661 94.5
43.303402 94.6
43.622541 94.7
43.944096 94.8
44.268088 94.9
44.594535 95
44.923459 95.1
45.254880 95.2
45.588818 95.3
45.925294 95.4
46.264329 95.5
46.605944 95.6
46.950159 95.7
47.296998 95.8
47.646481 95.9
47.998630 96
48.353468 96.1
48.711016 96.2
49.071287 96.3
49.434277 96.4
49.800043 96.5
50.168610 96.6
50.539999 96.7
50.914234 96.8
51.291340 96.9
51.671340 97
52.054258 97.1
52.440119 97.2
52.828947 97.3
53.220767 97.4
53.615605 97.5
54.013485 97.6
54.414432 97.7
54.818474 97.8
55.225634 97.9
55.635941 98
56.049420 98.1
56.466098 98.2
56.886001 98.3
57.309158 98.4
57.735595 98.5
58.165340 98.6
58.598421 98.7
59.034866 98.8
59.474704 98.9
59.917963 99
60.364672 99.1
60.814860 99.2
61.268557 99.3
61.725792 99.4
62.186596 99.5
62.650999 99.6
63.119031 99.7
63.590723 99.8
64.066106 99.9
64.545212 100
65.028072 100.1
65.514718 100.2
66.005182 100.3
66.499497 100.4
66.997696 100.5
67.499811 100.6
68.005876 100.7
68.515926 100.8
69.029993 100.9
69.548112 101
70.070318 101.1
70.596645 101.2
71.127130 101.3
71.661806 101.4
72.200673 101.5
72.743839 101.6
73.291107 101.7
73.842711 101.8
74.398729 101.9
74.958982 102
75.523511 102.1
76.092524 102.2
76.666060 102.3
77.244157 102.4
77.826854 102.5
78.414201 102.6
79.006061 102.7
79.602640 102.8
80.203979 102.9
80.810119 103
81.421102 103.1
82.037093 103.2
82.657737 103.3
83.283406 103.4
83.913866 103.5
84.549380 103.6
85.189990 103.7
85.835741 103.8
86.486677 103.9
87.142844 104
87.804285 104.1
88.471047 104.2
89.143188 104.3
89.820562 104.4
90.503395 104.5
91.191735 104.6
91.885629 104.7
92.585126 104.8
93.290273 104.9
94.001121 105
94.717850 105.1
95.440072 105.2
96.168144 105.3
96.902187 105.4
97.641959 105.5
98.387734 105.6
99.139562 105.7
99.897496 105.8
100.661590 105.9
101.431936 106
102.208330 106.1
102.991046 106.2
103.780137 106.3
104.575659 106.4
105.377668 106.5
106.186221 106.6
107.001374 106.7
107.823207 106.8
108.651502 106.9
109.486567 107
110.328463 107.1
111.177250 107.2
112.032986 107.3
112.895733 107.4
113.765552 107.5
114.642505 107.6
115.526463 107.7
116.417678 107.8
117.316214 107.9
118.222136 108
119.135506 108.1
120.056391 108.2
120.984994 108.3
121.920762 108.4
122.864240 108.5
123.815497 108.6
124.774601 108.7
125.741620 108.8
126.716623 108.9
127.699714 109
128.690696 109.1
129.689872 109.2
130.697315 109.3
131.713096 109.4
132.737289 109.5
133.770101 109.6
134.811133 109.7
135.860799 109.8
136.919173 109.9
137.986332 110
139.062353 110.1
140.147474 110.2
141.241242 110.3
142.344106 110.4
143.456145 110.5
144.577440 110.6
145.708071 110.7
146.848230 110.8
147.997565 110.9
149.156483 111
150.325067 111.1
151.503403 111.2
152.691575 111.3
153.889655 111.4
155.097539 111.5
156.315519 111.6
157.543815 111.7
158.782105 111.8
160.030767 111.9
161.289648 112
162.559073 112.1
163.839133 112.2
165.129880 112.3
166.431412 112.4
167.743661 112.5
169.066920 112.6
170.401284 112.7
171.746948 112.8
173.103583 112.9
174.471616 113
175.851147 113.1
177.242277 113.2
178.645038 113.3
180.059430 113.4
181.485726 113.5
182.924029 113.6
184.374391 113.7
185.836782 113.8
187.311496 113.9
188.798795 114
190.298233 114.1
191.810318 114.2
193.335160 114.3
194.872939 114.4
196.423380 114.5
197.986913 114.6
199.563651 114.7
201.153644 114.8
202.756885 114.9
204.373679 115
206.004147 115.1
207.648141 115.2
209.306042 115.3
210.978006 115.4
212.663948 115.5
214.363888 115.6
216.078247 115.7
217.806833 115.8
219.550062 115.9
221.308039 116
223.080666 116.1
224.868481 116.2
226.671016 116.3
228.488835 116.4
230.322114 116.5
232.170623 116.6
234.034812 116.7
235.914708 116.8
237.810388 116.9
239.722137 117
241.649850 117.1
243.593972 117.2
245.554271 117.3
247.531295 117.4
249.524740 117.5
251.535206 117.6
253.562366 117.7
255.606821 117.8
257.668260 117.9
259.747214 118
261.843499 118.1
263.957563 118.2
266.089321 118.3
268.239072 118.4
270.406979 118.5
272.592931 118.6
274.797399 118.7
277.020447 118.8
279.261979 118.9
281.522484 119
283.802007 119.1
286.100689 119.2
288.418433 119.3
290.755797 119.4
293.112774 119.5
295.489605 119.6
297.886109 119.7
300.302678 119.8
302.739465 119.9
305.196455 120 ] ;

SonePhone = SonePhone';

out = interp1(SonePhone(1,:), SonePhone(2,:), in, 'cubic');
