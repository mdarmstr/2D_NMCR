# Importing the statistics module
import math
import scipy.optimize
import numpy as np
from scipy.stats import chi2
# here is the transalter version of clust resolution from matlab to python
# here the input type for  data clust1 and clust2 is matrix
def clustResolution(clust1, clust2):
    # make the array list data clust1 and clust2 matrix
    c1 = clust1
    c2 = clust2


    # find the center points of the data
    c1_mean = c1.mean(axis=0)
    c2_mean = c2.mean(axis=0)

    #SVD for ranked eigenvalues from covariance matrices
    c1_cov = np.cov(c1.T)
    c2_cov = np.cov(c2.T)
    eivec1, eival1, Vh1 = np.linalg.svd(c1_cov)
    eivec2, eival2, Vh2 = np.linalg.svd(c2_cov)

    eival1 = np.diag(eival1)
    eival2 = np.diag(eival2)

    # print("###############")

    # print(eivec1)
    #angle of rotation
    phi1 = math.atan2(eivec1[0][1], eivec1[0][0])
    phi2 = math.atan2(eivec2[0][1], eivec2[0][0])

    #Conversion from 4quad to standard form for angles
    if phi1 < 0:
        phi1 = phi1 + 2*math.pi
    if phi2 < 0:
        phi2 = phi2 + 2*math.pi

    #Reorganizing vectors into individual values for anonymous function



    eival1_1 = eival1[0][0]
    eival1_2 = eival1[1][1]
    eival2_1 = eival2[0][0]
    eival2_2 = eival2[1][1]

    d21x = c2_mean[0] - c1_mean[0]
    d21y = c2_mean[1] - c1_mean[1]
    #Cost function
    # sqrt((-(sqrt(eival21)*cos(theta(2))*sin(phi2) +
    # sqrt(eival22)*sin(theta(2))*cos(phi2))*d21x/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) -
    # sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) +
    # sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) -

    # sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) +

    # sqrt(eival12)*sin(theta(1))*cos(phi1)))) + (sqrt(eival21)*cos(theta(2))*cos(phi2) -

    # sqrt(eival22)*sin(theta(2))*sin(phi2))*d21y/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) -

    # sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) +

    # sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) -

    # sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) +

    # sqrt(eival12)*sin(theta(1))*cos(phi1)))))^2 + (-(sqrt(eival11)*cos(theta(1))*sin(phi1) +

    # sqrt(eival12)*sin(theta(1))*cos(phi1))*d21x/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) -

    # sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) +

    # sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) -

    # sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) +

    # sqrt(eival12)*sin(theta(1))*cos(phi1)))) + (sqrt(eival11)*cos(theta(1))*cos(phi1) -

    # sqrt(eival12)*sin(theta(1))*sin(phi1))*d21y/((-(sqrt(eival11)*cos(theta(1))*cos(phi1) -

    # sqrt(eival12)*sin(theta(1))*sin(phi1))*(sqrt(eival21)*cos(theta(2))*sin(phi2) +

    # sqrt(eival22)*sin(theta(2))*cos(phi2))) + ((sqrt(eival21)*cos(theta(2))*cos(phi2) -

    # sqrt(eival22)*sin(theta(2))*sin(phi2))*(sqrt(eival11)*cos(theta(1))*sin(phi1) +
    # sqrt(eival12)*sin(theta(1))*cos(phi1)))))^2);

    costfn = lambda theta: math.sqrt((-(math.sqrt(eival2_1)*math.cos(theta[1])*math.sin(phi2)+
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.cos(phi2))*d21x/((-(math.sqrt(eival1_1)*math.cos(theta[0])*math.cos(phi1)-
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.sin(phi1))*(math.sqrt(eival2_1)*math.cos(theta[1])*math.sin(phi2)+
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.cos(phi2)))+((math.sqrt(eival2_1)*math.cos(theta[1])*math.cos(phi2)-
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.sin(phi2))*(math.sqrt(eival1_1)*math.cos(theta[0])*math.sin(phi1)+
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.cos(phi1)))) + (math.sqrt(eival2_1)*math.cos(theta[1])*math.cos(phi2)-
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.sin(phi2))*d21y/((-(math.sqrt(eival1_1)*math.cos(theta[0])*math.cos(phi1)-
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.sin(phi1))*(math.sqrt(eival2_1)*math.cos(theta[1])*math.sin(phi2)+
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.cos(phi2)))+((math.sqrt(eival2_1)*math.cos(theta[1])*math.cos(phi2)-
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.sin(phi2))*(math.sqrt(eival1_1)*math.cos(theta[0])*math.sin(phi1)+
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.cos(phi1)))))**2 +(-(math.sqrt(eival1_1)*math.cos(theta[0])*math.sin(phi1)+
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.cos(phi1))*d21x/((-(math.sqrt(eival1_1)*math.cos(theta[0])*math.cos(phi1)-
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.sin(phi1))*(math.sqrt(eival2_1)*math.cos(theta[1])*math.sin(phi2)+
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.cos(phi2)))+((math.sqrt(eival2_1)*math.cos(theta[1])*math.cos(phi2) -
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.sin(phi2))*(math.sqrt(eival1_1)*math.cos(theta[0])*math.sin(phi1)+
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.cos(phi1))))+(math.sqrt(eival1_1)*math.cos(theta[0])*math.cos(phi1) -
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.sin(phi1))*d21y/((-(math.sqrt(eival1_1)*math.cos(theta[0])*math.cos(phi1)-
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.sin(phi1))*(math.sqrt(eival2_1)*math.cos(theta[1])*math.sin(phi2)+
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.cos(phi2)))+((math.sqrt(eival2_1)*math.cos(theta[1])*math.cos(phi2)-
                          math.sqrt(eival2_2)*math.sin(theta[1])*math.sin(phi2))*(math.sqrt(eival1_1)*math.cos(theta[0])*math.sin(phi1)+
                          math.sqrt(eival1_2)*math.sin(theta[0])*math.cos(phi1)))))**2)

    #fminsearch
    exit = 5
    iter = 1
    while exit != 1 and iter < 100:
        pi_matrix = np.array([2*math.pi, ((2*math.pi + math.pi)/2)])
        theta0 = np.array(np.random.rand(1)) * pi_matrix.T
        # fminsearch options

        minimum = scipy.optimize.fmin(costfn, theta0, xtol=10**-5, maxfun=10**3, disp=False)
        min = costfn(minimum)
        # print("minimum")
        chisq = (min**2)/2
        # changed the cdf -> sf
        confidencelimit = chi2.cdf(chisq, 2)

        iter = iter + 1
        if iter == 100:
            confidencelimit = 0
        # print(confidencelimit)
        exit = 1
    return confidencelimit