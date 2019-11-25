# A nonconvex optimization approach to IMRT planning with dose-volume constraints

Fluence map optimization for intensity-modulated radiation therapy planning can be formulated as a large-scale inverse problem with multi-objectives on the tumors and organs-at-risk. Unfortunately, clinically relevant dose-volume constraints are nonconvex, so convex formulations and algorithms cannot be directly applied to the problem. We propose a novel approach to handle dose-volume constraints while preserving their nonconvexity, as opposed to previous efforts which focused on iterative convexification. The proposed method is amenable to efficient algorithms based on partial minimization and naturally adapts to handle maximum and mean dose constraints, which are prevalent in current practice, and cases of infeasibility. We demonstrate our approach using the CORT dataset, and show that it is easily adaptable to radiation treatment planning with dose-volume constraints for multiple tumors and organs-at-risk.

## Documents
* [Preprint](https://arxiv.org/abs/1907.10712) (July 2019)
