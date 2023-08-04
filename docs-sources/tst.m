clear
clc

m = 3;
n = 10;
p = 5;

C = [0 1]';
CB = [1 0]';
q = size(C, 1);

analysis1 = Analysis.leaveOneSessionOut(m, C);
analysis1.addPermutations();
analysis2 = Analysis.leaveOneSessionOut(m, C, CB);
analyses = {analysis1, analysis2};

Ys = cell(1, m);
Xs = cell(1, m);
for i = 1 : m
    Ys{i} = randn(n, p);
    Xs{i} = randi(2, n, q) - 1;
end

ccm = CvCrossManova(Ys, Xs, analyses)

% ccm.dispAnalyses()
% % ccm.showAnalyses()
% 
% ccm.runAnalyses()


% sessions k = 1 : m
% folds l = 1 : L
% 'training' sessions a in A_l
% 'validation' sessions b in B_l
% analyses i = 1 : nAnalyses
% permutations j = 1 : nPerms