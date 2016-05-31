---
title: "DDP Project"
author: "Arinze Akutekwe"
date: "30 May 2016"
output: html_document
---

# Comprehensive Learning Particle Swarm Optimization

## Synopsis:
The Comprehensive Learning Particle Swarm Optimization was introduced by J.J Liang et al (2006) and is a variant of the popular Particle Swarm Optimization by Kennedy and Eberhart (1995); a population-based
stochastic optimization technique for global optimization of multimodal functions. I implemented the Algorithm in R for a different optimization scenario in one of my papers during my PhD. Citation [below](http://tinyurl.com/z79yctn):

A. Akutekwe and H. Seker, "Inferring the dynamics of gene regulatory networks via optimized recurrent neural network and dynamic Bayesian network," Computational Intelligence in Bioinformatics and Computational Biology (CIBCB), 2015 IEEE Conference on, Niagara Falls, ON, 2015, pp. 1-8.

In simple terms the algorithm performs a minimization towards the global best solution in a D-dimensional search space. In this project, the user can input the desired number of iterations, number of populations and number of dimensions. Under the hood, a Toroidal saturation function keeps particles constrained within a pre-defined search space.

The algorithm outputs initial fitness, initial solutions, global best fitness and global best solutions and a plot of the best function values with the number of iterations.

To use the algorithm, select the desired number of population, number of dimensions and number of iterations. The select how many iterations you want and press the run button. The results are display with the resulting minimization plot result.
