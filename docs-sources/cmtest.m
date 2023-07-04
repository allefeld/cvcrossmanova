% ---
% title: Test of Cross-MANOVA using Examples
% author: Carsten Allefeld
% date: 2023-3-23
% jupyter:
%   jupytext:
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.14.5
%   kernelspec:
%     display_name: Matlab
%     language: matlab
%     name: matlab
% ---

%| echo: false
clear all %#ok<CLALL>
format short

% # Perfect pattern stability
%
% ## Define true situation and simulate data
%
% Create patterns across two variables for the difference between conditions 1 & 2 and 3 & 4. They are identical, i.e. there is perfect pattern stability.

b21 = [0 1];
b43 = [0 1];

% Define error covariance matrix as identity matrix.

Sigma = eye(2);

% Define contrasts for the difference between conditions 1 & 2 and 3 & 4.

C21 = [-1 1 0 0 0]';
C43 = [0 0 -1 1 0]';

% We repeat the basic design matrix in each session many times, to well estimate pattern distinctness and pattern stability.

nRep = 1000;

% -   Create an indicator-style design matrix for a design with four conditions and a constant regressor.
% -   Create parameters for the five regressors across two variables.
% -   Create data with four sessions, with the same design matrix and parameters.

% +
X = repmat([eye(4) ones(4, 1)], nRep, 1);
[n, q] = size(X);

b1 = [1, 2];
b2 = b1 + b21;
b3 = [3, 4];
b4 = b3 + b43;
b5 = [5, 6];
B = [b1; b2; b3; b4; b5];
[q, p] = size(B);

m = 4;
Ys = cell(m, 1);
Xs = cell(m, 1);
for k = 1 : m
  Ys{k} = X * B + mvnrnd(zeros(1, p), Sigma, n);
  Xs{k} = X;
end
% -

% ## Define analyses
%
% 1)  A leave-one-session-out cross-validated MANOVA with contrast `C21`.
% 2)  A leave-one-session-out cross-validated MANOVA with contrast `C43`.
% 3)  A leave-one-session-out cross-validated cross-MANOVA from contrast `C21` to contrast `C43`.
% 4)  A leave-one-session-out cross-validated cross-MANOVA from contrast `C43` to contrast `C21`.
% 5)  A custom analysis with two folds, where in fold 1 contrast `C21` is applied to the first two sessions and contrast `C43` to the last two sessions, and vice versa in fold 2.

analysis1 = {C21};
analysis2 = {C43};
analysis3 = {C21, C43};
analysis4 = {C43, C21};
analysis5 = Analysis(C21, C43, [1 1 0 0 ; 0 0 1 1], [0 0 1 1 ; 1 1 0 0]);
analyses = {analysis1, analysis2, analysis3, analysis4, analysis5};

% The custom analysis is specified as an object of class `Analysis`:

analysis5

% An `Analysis` object has a method `show` which graphically displays the analysis.

analysis5.show;

% ## Prepare analyses
%
% Create a `CvCrossManova` object, passing data matrices, design matrices, and analysis specifications.

cm = CvCrossManova(Ys, Xs, analyses)

% Internally, it converts all analysis specifications into `Analysis` objects.

cm.analyses

% Analysis details can be listed textually or graphically.

cm.dispAnalyses

cm.showAnalyses

% ## Run analyses

cm.runAnalyses()

% The true values for the five analyses are:
%
% 1) pattern distinctness 0.125
% 2) pattern distinctness 0.125
% 3) pattern stability 0.125
% 4) pattern stability 0.125
% 5) pattern stability 0.125
%
% # Perfect pattern instability
%
% ## Define modified true situation and simulate data
%
% Create patterns across two variables for the difference between conditions 1 & 2 and 3 & 4. They are orthogonal, the first only involving variable 1 and the second only involving variable 2.

b21 = [1 0];
b43 = [0 1];

% +
X = repmat([eye(4) ones(4, 1)], nRep, 1);
[n, q] = size(X);

b1 = [1, 2];
b2 = b1 + b21;
b3 = [3, 4];
b4 = b3 + b43;
b5 = [5, 6];
B = [b1; b2; b3; b4; b5];
[q, p] = size(B);

m = 4;
Ys = cell(m, 1);
Xs = cell(m, 1);
for k = 1 : m
  Ys{k} = X * B + mvnrnd(zeros(1, p), Sigma, n);
  Xs{k} = X;
end
% -

% ## Prepare and run analyses

cm = CvCrossManova(Ys, Xs, analyses)
cm.runAnalyses()

% The true values for the five analyses are:
%
% 1) pattern distinctness 0.125
% 2) pattern distinctness 0.125
% 3) pattern stability 0
% 4) pattern stability 0
% 5) pattern stability 0
%
%
% # Orthogonal patterns with a non-identity error covariance matrix
%
%
% ## Define modified true situation and simulate data
%
% Modify the error covariance matrix to include a positive correlation between the two variables.

Sigma = [1 0.5; 0.5 1]

% +
X = repmat([eye(4) ones(4, 1)], nRep, 1);
[n, q] = size(X);

b1 = [1, 2];
b2 = b1 + b21;
b3 = [3, 4];
b4 = b3 + b43;
b5 = [5, 6];
B = [b1; b2; b3; b4; b5];
[q, p] = size(B);

m = 4;
Ys = cell(m, 1);
Xs = cell(m, 1);
for k = 1 : m
  Ys{k} = X * B + mvnrnd(zeros(1, p), Sigma, n);
  Xs{k} = X;
end
% -

% ## Prepare and run analyses

cm = CvCrossManova(Ys, Xs, analyses)
cm.runAnalyses()

% The true values for the five analyses are:
%
% 1) pattern distinctness 0.1667
% 2) pattern distinctness 0.1667
% 3) pattern stability -0.0833
% 4) pattern stability -0.0833
% 5) pattern stability -0.0833
