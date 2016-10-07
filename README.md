# PHYS220 CW 6

**Author(s):** _\<your name(s)\>_

[![Build Status](https://travis-ci.org/chapman-phys220-2016f/cw-06-saktill.svg?branch=master)](https://travis-ci.org/chapman-phys220-2016f/cw-06-saktill)

**Due date:** 2016/10/11

## Specification

**Reminder: We have switched to Python3 officially.**

1. In a python module ```calculus.py``` create a vectorized version of a central finite difference. That is, given a discretized function $f = (f_1, \ldots, f_N)$ defined as $f_i = f(x_i)$ on a grid of points $x_i$ from $i=0$ to $i=N$ with fixed spacing $x_{i+1} - x_i = h$, the central difference at the point $f_i$ is: $$Df_i = (f_{i+1} - f_{i-1}) / (2h)$$. Implement this central difference operation as a matrix product of a derivative matrix $D$ with $f$ such that the central difference formula is satisfied for all non-boundary points of $f$. What dimensions must this matrix have to produce the correct answer? How should you handle the boundary points? Test your implementation to make sure it works as you expect. 
1. Create a vectorized version of a second-order central finite difference (i.e., repeat the first-order finite difference formula from above twice - what second-order formula do you get?). Implement this second-order derivative as a matrix product of a second-order derivative matrix $D2$ with $f$. Is $D2$ the same thing as squaring the matrix $D$ from the previous exercise?
1. Create a vectorized version of the trapezoidal integration rule. That is, given a discretized function $f$ defined as in the previous two exercises, define the integral as a matrix product of an integral matrix $I$ with $f$ such that $$If = \int f(x) dx$$, using the trapezoid rule for the discrete approximation to the integral. How should you handle the boundary points? What dimensions must this matrix have to produce the correct answer? 
1. Read the following abstract from the journal article [Diabetes Care 17, 152-154 (1994)](http://care.diabetesjournals.org/content/17/2/152.abstract).  Based on the abstract (and [paper](TaisMethod.pdf)), explain in your own words what the results of the author are, and how they relate to the exercises in this assignment. Note that according to the [ISI Web of Science](http://apps.webofknowledge.com/), this article has been cited 213 times since 1994, according to the chart below:
![Citation chart](citations.jpg)

To cleanly present your work, create a Jupyter notebook ```cw6-discretecalc.ipynb``` that imports ```calculus.py``` and carefully discusses your answers to each question. The notebook should be formatted professionally using Markdown headings, including your name and a brief description of each homework section. 

Pro-tip: using git to manage conflicts on Jupyter notebooks is a pain. I recommend delegating one person from your group to edit the notebook, to avoid merge conflicts.

## Assessment

Analyze in this section what you found useful about this assignment in your own words. Include any lingering questions or comments that you may have.

**CHANGEME**

## Honor Pledge

I pledge that all the work in this repository is my own with only the following exceptions:

* Content of starter files supplied by the instructor;
* Code borrowed from another source, documented with correct attribution in the code and summarized here.

Signed,

Sakthi and Will
