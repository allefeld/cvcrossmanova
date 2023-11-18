# TODO


-   new function which creates the searchlight
    as a logical 3d sprite with an anchor point:
    `searchlight(radius, mat = eye(3))`.
    `mat` is extension of functionality to create a sphere in physical space.

    for arbitrary mat, the sphere is an ellipsoid in voxel index space
    how to compute its bounding box? are the 6 cardinal points enough?
    https://math.stackexchange.com/questions/3926884/smallest-axis-aligned-bounding-box-of-hyper-ellipsoid

    this is called by `ccmSearchlight`
    the nonzero elements are counted for display
    it is passed to runsearchlight to use as a running window

    the `searchlightSizes` function uses `searchlight`
    which therefore has to return the admissible radii (distances)

    If we do this, is the table under 'Searchlight radius and size' still useful? Maybe return the table from `ccmSearchlight` if the radius is not specified?

-   helpToMarkdown: automatic references?

-   Consistent capitalization of 'Cross-validated MANOVA'
    as well as 'Cross-validated Cross-MANOVA'
    consistent numerus for 'data'

-   Actually, Cross-MANOVA should be called "canonical correlation".

-   give ids to errors (asserts) and warnings

-   use strings vs character arrays consistently

-   use dedicated tests instead of `asserts` within the tutorial
    see MATLAB Testing Framework
    https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html

-   integrate documentation from cvmanova repository,
    especially explanation negative D-values
    (move from `cvmanova.m`)

-   do this again:
    The implementation contains a hard-coded limit on the number of voxels within a searchlight or ROI regardless of regularization, of 90% of the available error degrees of freedom. That is already a rather large threshold, which one should normally not get close to.

-   integrate new repository with old one?
    making cvcrossmanova an evolution of cvmanova

-   Implementation notes?
    classes are derived from `handle & matlab.mixin.Scalar`
    `handle`: prevents unnecessarily copying data and enables methods to change the object
    `matlab.mixin.Scalar`: prevents forming object arrays, for wich the classes are not designed

-   in documented parameter list. specification of default value with = suggests that it is a keyword parameter. fix. by making them actually keyword parameters?

-   see also / cross references to and within `reference.qmd`

-   hide quasi-private methods, indicated by "_"

-   check whether regressors have the same meaning across `sessionsA` and `sessionsB`

-   "of logical 3D volumes or filenames specifying region masks"
    → "of logical region masks specified as three-dimensional arrays or filenames"

-   "to apply whitening and high-pass filtering (set up in SPM)"
    → "to apply the whitening and high-pass filtering set up in SPM"

-   `Ds` add ", either a scalar value or an array of permutation values."

-   `preparation.qmd`: do regions masks need to be aligned to the anatomical or the mean BOLD image?

-   Check Thomas' similarity definitions
    ```
    D = [AA AB BA BA];
    test_signal_estimate(jj,1) = sqrt(abs(D(4,jj) / D(1,jj)));
    pattern_similarity_estimate(jj,1) = ((D(2,jj) + D(3,jj)) / 2) / D(1,jj) /test_signal_estimate(jj,1);
    ```
    Email 2023–11–3
    