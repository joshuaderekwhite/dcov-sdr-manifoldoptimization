# Sequential SDR via DCOV Optimized on a Stiefel Manifold
## Inputs
- X: An nxp matrix containing numeric values (no limitations on n or p)
- Y: An nXd matrix containing numberic values (no limitations on n or d; may fail if d is too large)
- labels: p labels representing the column names of the X matrix

## Outputs
- Syx: The reduced subspace for the X data matrix
- beta: The prediction coefficients for the regression
- labels: The labels for the reduced X data matrix

## How to Run
- Requires manopt package for MATLAB/Octave to be loaded to local machine
- Save all .m files to local repository
- Run ssdr(X,Y,labels) to receive the reduced matrix
