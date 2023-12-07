# TODO

## important for first publication

-   fix bug in quartoToScript which de-indents code

-   `reference.qmd` / `helpToMarkdown`: references from method / function documentation to General Notes. Implemented as a number of search strings.

-   Consistent capitalization of 'Cross-validated MANOVA'
    as well as 'Cross-validated Cross-MANOVA'
    consistent numerus for 'data'

-   integrate documentation from <https://github.com/allefeld/cvmanova>
    consider moving explanation of negative D-values from `cvmanova.qmd` to `reference.qmd`

-   check whether regressors have the same meaning across `sessionsA` and `sessionsB`, and print information labeling analyses as cross or not based on that.

-   `preparation.qmd`: do regions masks need to be aligned to the anatomical or the mean BOLD image?

-   Check Thomas' similarity definitions
    ```
    D = [AA AB BA BA];
    test_signal_estimate(jj,1) = sqrt(abs(D(4,jj) / D(1,jj)));
    pattern_similarity_estimate(jj,1) = ((D(2,jj) + D(3,jj)) / 2) / D(1,jj) /test_signal_estimate(jj,1);
    ```
    Email 2023–11–3

-   release script: zip and attach documentation HTML files

## nice to have

-   give ids to errors (asserts) and warnings

-   use strings vs character arrays consistently

-   use dedicated tests instead of `asserts` within the tutorial
    see MATLAB Testing Framework
    https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html

-   do this again:
    The implementation contains a hard-coded limit on the number of voxels within a searchlight or ROI regardless of regularization, of 90% of the available error degrees of freedom. That is already a rather large threshold, which one should normally not get close to.
    or rather, check numerical / statistical stability

-   Implementation notes?
    classes are derived from `handle & matlab.mixin.Scalar`
    `handle`: prevents unnecessarily copying data and enables methods to change the object
    `matlab.mixin.Scalar`: prevents forming object arrays, for wich the classes are not designed

-   hide quasi-private methods, indicated by "_"
    Seems not to be necessary, except maybe `checkEstimability`

-   nvargs for undocumented functions: `searchlight`, `loadDataSPM`, `spmReadVolMatched`, `spmWriteImage`
    generally, revise undocumented functions

-   Actually, Cross-MANOVA should be called "canonical correlation".
