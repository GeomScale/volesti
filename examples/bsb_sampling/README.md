## Explanation
This is a variant of Running Shake and Bake algorithm implemented in 'shake_and_bake_walk.hpp' that belongs to the class of boundary sampling algorithms. It follows the steps as described in [1], but after step 1 (direction sampling) it does number of reflections defined by upper bound of reflection (nr). For each step, the number of reflections is sampled by using inverse exponential distribution, i.e. the number of reflections is (1-z)*nr. 

[1] C. G. E. Boender, R. J. Caron, J. F. McDonald, A. H. G. Rinnooy Kan,H. E. Romeijn, R. L. Smith, J. Telgen i A. C. F. Vorst,  *Shake-And-Bake Algorithms for Generating Uniform Points on the Boundary of Bounded Polyhedra*, 1991.  
Available at: https://doi.org/10.1016/0166-218X(91)90006-7

