---
title       : Developing Data Products
subtitle    : Comprehensive Learning Particle Swarm Optimization
author      : Arinze Akutekwe
job         : Data Scientist
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : []            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---
## What is CLPSO?

 - A variant of Particle SWarm Optimization.
 
 - Used for global optimization of multimodal functions.
 
 - Was used for different optimization scenarios during PhD.
 
 - Performs much better when compared to others (in my opinion).
 
---

---
## About the package
 
 - Implements eight popular test problems. Allows user to enter desired number of iterations
 
 - The higher the number, the bigger the problem and hence more time!
 
 - User can also select number of populations (particles) and number of number of dimensions.
 
 ### Result Outputs
 
 - Initial Best Fitness
 
 - Initil Best Solution
 
 - Final Best Fitness
 
 - Final Best Solution
 
 - A plot of best values versus number of iterations.

---
## Discussion
 
 - Various parameters affect the output of the algorithm
 
- Many of the parameters have been fixed for simplicity, see [J.J Liang et al (2006)](http://tinyurl.com/hkvhtuj) for more details.
 
- Was implemented to fine-tune the parameters of the Elman Recurrent Neural Network with Dynamic Bayesian Network
 
- More details in the paper [A. Akutekwe and H. Seker (2015)](http://tinyurl.com/z79yctn)
 
 - The Schwefel function is complex, with many local minima and for a much larger dimensional space (-500,500)
 than what was fixed here in this work (-5.12,5.12) for simplicity. Hence it looks as if nothing happened.
 
---

---
## Conclusion
 
 - No free lunch in the world of optimization.
 
 - No "best" algorithm. Differential Evolution algorithm could perform better in some cases.
 
 - Simple test problems such as DeJobg and Parabola are minimized the most.
 
 - The code can be made better, faster and parallelized. Please contribute to it :)
 
---




