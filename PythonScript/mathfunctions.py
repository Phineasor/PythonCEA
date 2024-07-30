# Supporting math functions for things
# Might even eliminate the need for external math libraries
"""
@author: phineas
"""

import math


def sind(Theta):
    return math.sin((math.pi / 180) * Theta)


def cosd(Theta):
    return math.cos((math.pi / 180) * Theta)


def tand(Theta):
    return math.tan((math.pi / 180) * Theta)


def linspace(start, end, number):
    i = 0
    output = []
    dx = (end - start) / number

    while i < number:
        output.append(start + (i * dx))
        i += 1
    return output
