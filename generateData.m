function [data] = generateData(clusters)

%Please read LICENSE.md before using, and attribute Michael D. Sorochan Armstrong and Dr. James J. Harynuk
%
%Current Version 1.0 2020-10-02
%
%(c) M. Sorochan Armstrong 2020
%Generates random data for testing the cluster resolution determination in 2D_NMCR.m
%
%
%Usage: [data] = generateData(clusters)
%
%inputs:
%clusters: integer value for the number of clusters to be generated.
%
%output:
%data: a 200*Nx3 matrix containing 200 randomly generated data points with 2 dimensions per class. Column 3 contains class information.
%
%
%%Version History
%1.0 - First version uploaded on GitHub - 2020-10-02

for qq = 1:clusters
    
  Pi = rand(1)*2*pi;

  max_drift = 60;

  X_0 = [rand(1)*max_drift,rand(1)*max_drift];

  max_eig = 20;

  eig_maj = normrnd(0.5,0.1)*max_eig;
  eig_min = abs(normrnd(0.5,0.1)*eig_maj);

  ranx = normrnd(0,1,[2,200]).*[eig_maj; eig_min];

  data(:,:,qq) = [(X_0 + ranx'*[cos(Pi),-sin(Pi); sin(Pi), cos(Pi)]), repmat(qq,size(ranx,2),1)]; %#ok %generating a linear relationship here

end

data = permute(data, [1 3 2]);
data = reshape(data,[max(size(ranx))*qq,3]);

data(:,1:2) = (data(:,1:2) - mean(data(:,1:2)))./std(data(:,1:2));

[u,s,~] =svd(data(:,1:2));

data(:,1:2) = u*s;

end
