#!/usr/bin/python

import numpy as np
import math

class Params:
    def base_stack(self, atom1, atom2, atom3, siteId):
        if siteId[atom2] == 'A':
            if siteId[atom3] == 'A':
                return [13.810,3.58,100.13]
            elif siteId[atom3] == 'T':
                return [15.050,3.56,90.48]
            elif siteId[atom3] == 'G':
                return [13.320,3.85,104.39]
            elif siteId[atom3] == 'C':
                return [15.820,3.45,93.23]
        elif siteId[atom2] == 'T':
            if siteId[atom3] == 'A':
                return [9.250,4.15,102.59]
            elif siteId[atom3] == 'T':
                return [12.44,3.93,93.32]
            elif siteId[atom3] == 'G':
                return [9.58,4.32,103.70]
            elif siteId[atom3] == 'C':
                return [13.11,3.93,94.55]
        elif siteId[atom2] == 'G':
            if siteId[atom3] == 'A':
                return [13.76,3.51,95.45]
            elif siteId[atom3] == 'T':
                return [14.59,3.47,87.63]
            elif siteId[atom3] == 'G':
                return [14.77,3.67,106.36]
            elif siteId[atom3] == 'C':
                return [15.17,3.42,83.12]
        elif siteId[atom2] == 'C':
            if siteId[atom3] == 'A':
                return [9.25,4.15,102.69]
            elif siteId[atom3] == 'T':
                return [12.42,3.99,96.05]
            elif siteId[atom3] == 'G':
                return [8.83,4.34,100.46]
            elif siteId[atom3] == 'C':
                return [14.01,3.84,100.68]
        else:
            print("ERROR: No base found!\n")

    def cross_stack_1(self, atom2, atom5, siteId):
        if siteId[atom2] == 'A':
            theta3 = 110.92
            if siteId[atom5] == 'A':
                theta1 = 154.04
                sigma = 6.420
                eps = 2.139*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'T':
                theta1 = 158.77
                sigma = 6.770
                eps = 2.714*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'G':
                theta1 = 153.88
                sigma = 6.270
                eps = 2.772*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'C':
                theta1 = 157.69
                sigma = 6.840
                eps = 1.909*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom2] == 'T':
            theta3 = 110.92
            if siteId[atom5] == 'A':
                theta1 = 148.62
                sigma = 6.770
                eps = 2.714*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'T':
                theta1 = 155.05
                sigma = 7.210
                eps = 2.139*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'G':
                theta1 = 147.54
                sigma = 6.530
                eps = 2.485*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'C':
                theta1 = 153.61
                sigma = 7.080
                eps = 2.916*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom2] == 'G':
            theta3 = 120.45
            if siteId[atom5] == 'A':
                theta1 = 153.91
                sigma = 6.270
                eps = 2.772*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'T':
                theta1 = 155.72
                sigma = 6.530
                eps = 2.485*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'G':
                theta1 = 151.84
                sigma = 5.740
                eps = 3.693*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'C':
                theta1 = 157.80
                sigma = 6.860
                eps = 1.104*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
        elif siteId[atom2] == 'C':
            theta3 = 120.45
            if siteId[atom5] == 'A':
                theta1 = 152.04
                sigma = 6.840
                eps = 1.909*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'T':
                theta1 = 157.72
                sigma = 7.080
                eps = 2.916*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'G':
                theta1 = 151.65
                sigma = 6.860
                eps = 1.104*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]
            elif siteId[atom5] == 'C':
                theta1 = 154.49
                sigma = 6.790
                eps = 4.699*0.88
                return [theta3 * math.pi / 180.0, theta1 * math.pi /180.0, sigma, eps]


#Second cross stacking interaction parameters
    def cross_stack_2(self, atom4, atom6, siteId):
        if siteId[atom4] == 'A':
            if siteId[atom6] == 'A':
                theta2 = 116.34
                sigma = 5.580
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'T':
                theta2 = 119.61
                sigma = 6.140
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'G':
                theta2 = 115.19
                sigma = 5.630
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'C':
                theta2 = 120.92
                sigma = 6.180
                return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom4] == 'T':
            if siteId[atom6] == 'A':
                theta2 = 107.40
                sigma = 6.140
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'T':
                theta2 = 110.76
                sigma = 6.800
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'G':
                theta2 = 106.33
                sigma = 6.070
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'C':
                theta2 = 111.57
                sigma = 6.640
                return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom4] == 'G':
            if siteId[atom6] == 'A':
                theta2 = 121.61
                sigma = 5.630
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'T':
                theta2 = 124.92
                sigma = 6.070
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'G':
                theta2 = 120.52
                sigma = 5.870
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'C':
                theta2 = 124.88
                sigma = 5.660
                return [theta2 * math.pi /180.0, sigma]
        elif siteId[atom4] == 'C':
            if siteId[atom6] == 'A':
                theta2 = 112.45
                sigma = 6.180
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'T':
                theta2 = 115.43
                sigma = 6.640
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'G':
                theta2 = 110.51
                sigma = 5.660
                return [theta2 * math.pi /180.0, sigma]
            elif siteId[atom6] == 'C':
                theta2 = 115.80
                sigma = 6.800
                return [theta2 * math.pi /180.0, sigma]

