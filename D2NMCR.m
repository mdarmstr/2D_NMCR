function [confidencelimit] = D2NMCR(clust_1, clust_2)

%Please read LICENSE.md before using, and attribute Michael D. Sorochan Armstrong and Dr. James J. Harynuk
%
%Current Version 1.0 2020-10-02
%
%(c) M. Sorochan Armstrong 2020; Inspired by previous code: LA_gaschrom_HoppingCRII.m
%An efficient numerical method for estimation of cluster resolution for a
%series of points in two dimensions, usually PCA scores.
%
%Usage: [confidencelimit] = 2D_NMCR(clust_1,clust_2)
%
%inputs:
%clust_1: N1x2 matrix of N1 points attributed to cluster 1 with 2 double-precision values describing their positions in 2D space.
%clust_2: N2x2 matrix of N2 points attributed to cluster 2 with 2 double-precision values describing their positions in 2D space.
%
%Requires the Statistics and Machine Learning Toolbox for chi2cdf.m and the Optimization toolbox for fminsearch.m
%

%%Version History
%1.0 - First version uploaded on GitHub - 2020-10-02

%find the center points of the data
cntr1 = mean(clust_1); cntr2 = mean(clust_2);

%SVD for ranked eigenvalues from covariance matrices
covar1 = cov(clust_1); covar2 = cov(clust_2);
[eivec1, eival1] = svd(covar1); [eivec2, eival2] = svd(covar2);

%angle of rotation
phi1 = atan2(eivec1(1,2),eivec1(1,1)); phi2 = atan2(eivec2(1,2),eivec2(1,1));

%Conversion from 4quad to standard form for angles
if phi1 < 0
    phi1 = phi1 + 2*pi;
end

if phi2 < 0
    phi2 = phi2 + 2*pi;
end

%Reorganizing vectors into individual values for anonymous function
eival11 = eival1(1,1); eival12 = eival1(2,2);
eival21 = eival2(1,1); eival22 = eival2(2,2);

d21x = cntr2(1) - cntr1(1);
d21y = cntr2(2) - cntr1(2);

%Cost function
costfn = @(theta) sqrt((-(sqrt(eival21)*cos(theta(2))*sin(phi2) + sqrt(eival22)*sin(theta(2))*cos(phi2))*d21x/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) - sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) + sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) - sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) + sqrt(eival12)*sin(theta(1))*cos(phi1)))) + (sqrt(eival21)*cos(theta(2))*cos(phi2) - sqrt(eival22)*sin(theta(2))*sin(phi2))*d21y/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) - sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) + sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) - sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) + sqrt(eival12)*sin(theta(1))*cos(phi1)))))^2 + (-(sqrt(eival11)*cos(theta(1))*sin(phi1) + sqrt(eival12)*sin(theta(1))*cos(phi1))*d21x/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) - sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) + sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) - sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) + sqrt(eival12)*sin(theta(1))*cos(phi1)))) + (sqrt(eival11)*cos(theta(1))*cos(phi1) - sqrt(eival12)*sin(theta(1))*sin(phi1))*d21y/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) - sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) + sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) - sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) + sqrt(eival12)*sin(theta(1))*cos(phi1)))))^2);

%fminsearch options
options = optimset('TolX',1e-8,'MaxFunEvals',1e3);

%fminsearch
exit = 5; iter = 1;

while exit ~= 1 && iter < 100 %catch in case the solution does not converge. Will keep iterating until convergence is reached.
    
    theta0 = rand(1,1).*[2*pi, 2*pi+pi/2];  %Diffence of pi/2 for between theta1, theta2. Calculated from line between critical points.
    
    [~,fval,exit] = fminsearch(costfn,theta0,options);%options);
    
    chisq = (fval^2)/2;
    confidencelimit = chi2cdf(chisq,2);
    
    iter = iter + 1;
    
    if iter == 100
        disp('NM_CR : NM_CR failed to reach convergence.') %error message
        
        disp('NM_CR failed to converge, setting CR = 0')
        
        confidencelimit = 0;
        
    end
end

end
