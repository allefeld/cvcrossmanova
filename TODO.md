# TODO

-   helpToMarkdown: automatic references?

-   Consistent capitalization of 'Cross-validated MANOVA'
    as well as 'Cross-validated Cross-MANOVA'
    consistent numerus for 'data'

-   Actually, Cross-MANOVA should be called canonical correlation.

-   give ids to errors (asserts) and warnings

-   use strings vs character arrays consistently

-   Rename class `CvCrossManova` / clarify class `Analysis`
    use namespace / package?

-   use dedicated tests instead of `asserts` within the tutorial
    see MATLAB Testing Framework
    https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html

-   Two tutorials:
    -   High-level using SPM and searchlight and region functions,
        applied to Haxby data
    -   Low-level using class directly, using ERP CORE.
        
-   Do we need an option to omit high-pass filtering & whitening?

-   pre-render script transforms `qmd` to `m` with
    `jupytext file.qmd --to m`
    can we attach these to the pages?

-   integrate documentation of cvmanova,
    especially explanation negative D-values

-   integrate new repository with old one?
    making cvcrossmanova an evolution of cvmanova

-   use `obj` instead of `self`

-   Implementation notes?
    classes are derived from `handle & matlab.mixin.Scalar`
    `handle`: prevents unnecessarily copying data and enables methods to change the object
    `matlab.mixin.Scalar`: prevents forming object arrays, for wich the classes are not designed

-   in documented parameter list. specification of default value with = suggests that it is a keyword parameter. fix. by making them actually keyword parameters?
