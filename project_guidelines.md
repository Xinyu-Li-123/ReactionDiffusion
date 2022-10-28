# Modeling and Simulation in Science, Engineering, and Economics, Fall 2022

## PROJECT WRITEUP GUIDELINES

The writeup should be one self-contained document, not a collection of files.
Any still figures should be part of the document, and any animations should
be included by providing links to the animation that the reader can use
to view the animation.   The reader should be able to view and/or
download the animation without going through any login procedure.

Title Page --state the **project title**, with **full names** and **email** 
             addresses of the authors.  Be sure to give your project 
             a descriptive title, not something like "Project 1".

Introduction -- explain **the problem** and **how you are going to address it**.

Equations -- state **the equations that govern the system** you are simulating.
             You do not need to derive the equations, but **say what the
             variables and parameters mean**, and **describe in words what
             each equation says**.  Besides the fundamental governing
             equations, there may be **some mathematical or geometric
             problems that you needed to solve to set up you model**,
             and this would be a good place to describe them.
             It is OK to write the equations by hand instead of
             typesetting them, but please do the writing yourself,
             don't scan my notes.  Try to get in the habit of thinking
             about what the equations mean as you write them.

Numerical Method -- this will look a lot like the equations, but involving
             **finite differences instead derivatives**.  In some projects,
             e.g., stochastic simulation, this may look more like an
             algorithm.  In any case, this section should **use mathematical
             language, not computer code, to describe what the computer
             will be doing as it simulates your system**.

Validation -- in this section, provide **evidence that your code is
             correct**.  If you can identify **quantities that should be
             conserved**, check that they are.  If there is a special
             case of your problem for which **an exact solution** is known,
             run your code for that case, and compare your results
             to the known exact solution.  **If your code has a timestep
             parameter, make sure it is small enough that the results
             are essentially independent of that parameter**.  If you are
             doing stochastic simulation, there is no timestep parameter,
             but you can compare long-time averages to the predictions
             obtained by considering microscopic equilibrium, or you
             can scale up and compare results to macroscopic theory.

Results and Discussion -- report **what computer experiments you did**,
             and **what you found**.  Use **graphs** where possible, generated
             by your computer program, in preference to tables of numbers.
             If your project produced **animations**, include links to the
             animations themselves, but also **try to capture the
             essence of what happened in some still figures**, e.g., by
             plotting the trajectories over time of some selected points.
             The results should include **at least one study in which you
             vary some parameter systematically to see  what happens**.
             To the extent that you can, **explain your findings**.

Summary and Conclusions -- write **a short summary of what you did and
             what you found**.  **Avoid discussions of future work** -- if there
             is something you would really like to do in the future, get
             it done now!

References -- be careful to **give full credit for anything you use from
             any source, including online sources**.  **This applies especially
             to computer code**.  If your project involves code written by
             anyone else, you should not only give credit but also make
             it clear in your writeup what you took from the source and
             what you added to it or how you modified it as you did your
             project.  Put a statement about this, if applicable, in the
             Introduction, with corresponding references in the reference
             section.

Computer code as an Appendix -- include some version of your **code, with
             comments** to make it clear that you understand what the
             code is doing.  If there are multiple versions, choose one
             that is representative (unless there are important differences
             between versions that call for explanation.)

**As you work on the writeup, keep in mind a reader who is *not* in our
class.  Your mission is to make your project understandable to such
a reader.**

Downloading of figures and diagrams from the internet is discouraged!
The figures in your writeup should be the results of your simulations,
or explanatory diagrams that you made yourself, one way or another.
Drawing a diagram by hand and scanning it is fine.  Likewise for
equations, but in both cases do the drawing or writing yourself.

Although the topic of the project is up to you, if you choose
something that we did not discuss in class, make sure you
understand the topic thoroughly and can explain what you are
doing (both to me and to the class) in a coherent way.

GRADING:
The Grader will grade each writeup on a scale of 0 - 25 points,
and I will grade each project on a scale of 0 - 25 points,
so that the maximum score for each project will be 50 points,
and since there are two projects to be done during the semester,
your maximum overall score will be 100 points.

Note that the Grader will only grade the writeup, but my grade will
be more of an overall impression of the project.  I will use the
writeups, however, at least to remind me what the project was about.
My grading will probably be done at the end of the semester, but
you should get more timely feedback from the Grader.

You can submit a revised writeup to the Grader, and the revision will
be graded as well.  What will count, then, as the grade assigned by
the Grader will be the average of the original grade and the grade
assigned to the revision.  Only one revision of any one writeup will
be allowed.  Revisions will only be accepted during the semester,
i.e., up to the last day of exam period, so it may not be practical
submit a revision on the writeup of project 2.  The Grader is under no
obligation to give you time to revise writeup 2, but if there happens
to be enough time, you may do so.

Late submission will not affect the grade assigned to a writeup,
but the date of all submissions will separately be noted and may affect
your final grade for the course.  Also, no work will be accepted
after the last day of exam period.  (If there is a medical reason
why you cannot complete the word during the semester, then you
need to *request* an Incomplete.)

RULES ABOUT PROJECT TEAMS:
Maximum team size = 2 People.
A team submits one writeup, by both team members, and makes
one joint presentation, in which both members participate.
Both members of a team receive the same grade for the project.
No two people are allowed to work together on both projects.

--C. Peskin
