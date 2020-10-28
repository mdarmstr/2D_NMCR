# 2D_NMCR
2-Dimensional Nedler-Mead simplex Cluster Resolution algorithm. Returns a value between 0 and 1, where a value of 1 describes two clusters which are completely resolved. A metric of model quality that can be used in feature selection routines for variable subsets projected into a linear subspace.

Created by Michael D. Sorochan Armstrong and Dr. James J. Harynuk at the University of Alberta

To test the D2NMCR.m use the following lines of code in a directory containing generateData.m and D2NMCR.m:

data = generateData(2);

D2NMCR(data(data(:,3)==1,1:2),data(data(:,3)==2,1:2));

A value for cluster resolution between 0 and 1 should display.
