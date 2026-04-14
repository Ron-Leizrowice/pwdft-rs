# Proposal: Replace hand-rolled special functions with validated libraries

**Status: COMPLETED (Option D — extended in-house, no new dependencies)**

## Result

Extended all three special functions to production quality without adding dependencies:

1. **Spherical Bessel j_l(x)**: Upward recurrence from j_0, j_1 for arbitrary l (was l <= 3, panicked otherwise). Tested against scipy for l=0-5.

2. **Legendre P_l(x)**: Bonnet's recurrence for arbitrary l (was l <= 3). Tested against exact polynomials for l=4-6, orthogonality validated to l=6, P_l(1)=1 and P_l(-1)=(-1)^l for l=0-10.

3. **erfc(x)**: Replaced Abramowitz & Stegun 5-term approximation (~1e-6 accuracy) with Chebyshev expansion (~1e-14). Taylor series for |x| < 0.5, Numerical Recipes Chebyshev for |x| >= 0.5.

No new dependencies added. All functions remain in their original files.

## Key decision: recurrence over library

The recurrence relations are 10-15 lines each, well-understood numerically, and called O(n_pw²) times in the hot path. A library call would add function pointer overhead without improving accuracy. The erfc Chebyshev expansion is 28 coefficients — more code than a library wrapper, but avoids a dependency for a single function.
