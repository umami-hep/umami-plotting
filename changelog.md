# Changelog

### [Latest]

### [v0.2.5] (2023/05/12)

- Update requirements to include atlas-ftag-tools v0.1.3 [!180](https://github.com/umami-hep/puma/pull/180)
- Use `VarVsVar` as a base class for `VarVsEff` [!179](https://github.com/umami-hep/puma/pull/179)
- Adding fraction scans to high level API [!178](https://github.com/umami-hep/puma/pull/178)
- Update pre-commit (using ruff) and store tagger scores as structured array [!177](https://github.com/umami-hep/puma/pull/177)
- Remove dev image [!176](https://github.com/umami-hep/puma/pull/176)
- Fix bug in ratio axis limits [!175](https://github.com/umami-hep/puma/pull/175)
- Add `VarVsVar` plot [!172](https://github.com/umami-hep/puma/pull/172)

### [v0.2.4] (2023/04/06)

- Replace `dijets` category with `QCD` category [!170](https://github.com/umami-hep/puma/pull/170)

### [v0.2.3] (2023/03/28)

- Integrate [atlas-ftag-tools](https://github.com/umami-hep/atlas-ftag-tools/) package [!168](https://github.com/umami-hep/puma/pull/168)
- HLAPI and CI Updates [!165](https://github.com/umami-hep/puma/pull/165)
- Extend format saving options [!160](https://github.com/umami-hep/puma/pull/160)

### [v0.2.2] (2023/02/28)

- Cast scores to full precision [!159](https://github.com/umami-hep/puma/pull/159)
- Add Xbb Support [!157](https://github.com/umami-hep/puma/pull/157)
- Improvements to the high level API [!155](https://github.com/umami-hep/puma/pull/155)
- Fixate the python container version [!153](https://github.com/umami-hep/puma/pull/153)
- Improve ROC format [#146](https://github.com/umami-hep/puma/pull/149)
- Fix for CI [!152](https://github.com/umami-hep/puma/pull/152)

### [v0.2.1] (2022/12/15)

- Change legend label of `dijets` [#146](https://github.com/umami-hep/puma/pull/146)

### [v0.2.0] (2022/12/09)

- Adding new high level API [#128](https://github.com/umami-hep/puma/pull/128)

### [v0.1.9] (2022/11/30)

- Adding boosted categories for Xbb to utils [!138](https://github.com/umami-hep/puma/pull/138)
- Running pylint also for tests [#133](https://github.com/umami-hep/puma/pull/133)
- Fix handling of nan values in histograms [#125](https://github.com/umami-hep/puma/pull/125)
- Adding support for under- and overflow bins in `puma.HistogramPlot` [#124](https://github.com/umami-hep/puma/pull/124)
- (Documentation) Adding copy-button to code cells in documentation [#131](https://github.com/umami-hep/puma/pull/131)

### [v0.1.8] (2022/08/30)

- Fix `set_ylim` in `puma.PlotBase` such that y-limits are correctly modified in cases with a y-offset [#119](https://github.com/umami-hep/puma/pull/119)
- Adding example for `puma.Line2DPlot` to the docs [#117](https://github.com/umami-hep/puma/pull/117)
- Adding support for ROC plots without ratio panels (was not possible until now) [#114](https://github.com/umami-hep/puma/pull/114)
- Lines with `label=None` (which is the default) will not appear in the legend anymore [#113](https://github.com/umami-hep/puma/pull/113)
- Adding new function `puma.utils.get_good_linestyles()` for easier linestyle management [#116](https://github.com/umami-hep/puma/pull/116)
- Adding the method `make_linestyle_legend()` which allows to specify an additional legend for linestyles [#113](https://github.com/umami-hep/puma/pull/113)

### [v0.1.7] (2022/08/10)

- Adding new option to place rejection label legend in ROC plots [#109](https://github.com/umami-hep/puma/pull/109)

### [v0.1.6] (2022/07/26)

- Adding support for weighted histograms (`puma.Histogram` now has an optional argument `weights`) [#86](https://github.com/umami-hep/puma/pull/86)
- Fixing bug where code crashed when histograms with discrete values + ratio panel were drawn [#99](https://github.com/umami-hep/puma/pull/99)
- Adding `h5py` to the Docker images [#97](https://github.com/umami-hep/puma/pull/97)
- Adding `transparent` attribute to `PlotObject` class. This allows to specify transparent background when initialising the plot [#96](https://github.com/umami-hep/puma/pull/96)

### [v0.1.5] (2022/07/05)

- Add `linewidth` and `alpha` to legend handles + set `alpha=1` by default (in `puma.Histogram`) [#92](https://github.com/umami-hep/puma/pull/92)
- Decreased default `figsize` for plots with zero or one ratio panels [#90](https://github.com/umami-hep/puma/pull/90)

### [v0.1.4] (2022/06/30)

- Renamed the `puma.FractionScan` and `puma.FractionScanPlot` classes to more general `puma.Line2DPlot` and `pumal.Line2D` [#84](https://github.com/umami-hep/puma/pull/84)
- Splitting `force` argument of `set_log()` method into `force_x` and `force_y` [#83](https://github.com/umami-hep/puma/pull/83)
- Adding `puma.PiePlot` class. Pie chart plots with `puma.HistogramPlot` are no longer possible [#70](https://github.com/umami-hep/puma/pull/70)
- Change default labels of `singlebjets` and `singlecjets` [#82](https://github.com/umami-hep/puma/pull/82)
- Support linestyles for variable vs. efficiency plots [#78](https://github.com/umami-hep/puma/pull/78)

### [v0.1.3] (2022/06/23)

- Adding more flavours to the global config  [#73](https://github.com/umami-hep/puma/pull/73)
- `ratio_group` in `puma.Histogram` objects can no longer be set via `flavour` argument [#74](https://github.com/umami-hep/puma/pull/74)
- Adding example for `plt.show` replacement + adding theme switcher button to docs [#72](https://github.com/umami-hep/puma/pull/72)
- Adding `atlas_tag_outside` and change default for `atlas_second_tag` [#71](https://github.com/umami-hep/puma/pull/71)
- Change default *bb*-jets colour to dark red and vlines to black [#69](https://github.com/umami-hep/puma/pull/69)
- Adding more general `ratio_group` argument to `puma.Histogram` [#67](https://github.com/umami-hep/puma/pull/67)
- Adding `calc_separation()` to `puma.metrics`, which allows to calculate the separation between two distributions [#27](https://github.com/umami-hep/puma/pull/27)
- Adding Zenodo link

### [v0.1.2] (2022/06/02/02)

- Adding automated coverage comment for pull request [#58](https://github.com/umami-hep/puma/pull/58)
- Fix that colour and legend label can be individually modified in case of flavoured histogram [#57](https://github.com/umami-hep/puma/pull/57)

### [v0.1.1] (2022/05/30)

- Adding documentation for updating the version switcher in the docs [#49](https://github.com/umami-hep/puma/pull/49)
- Adding version support in docs [#42](https://github.com/umami-hep/puma/pull/42)[#45](https://github.com/umami-hep/puma/pull/45)
- Adding development guidelines to the docs [#41](https://github.com/umami-hep/puma/pull/41)
- Adding `logx` [#40](https://github.com/umami-hep/puma/pull/40)
- Adding example page for the fraction scan plots [#38](https://github.com/umami-hep/puma/pull/38)
- Add warning for plotting normalised histograms without saying to in the ylabel [#34](https://github.com/umami-hep/puma/pull/34)
- Adding Docker images [#33](https://github.com/umami-hep/puma/pull/33)
- Remove warning when adding more than one reference histogram is case is allowed [#32](https://github.com/umami-hep/puma/pull/32)
- Update documentation [#24](https://github.com/umami-hep/puma/pull/24)

### [v0.1.0] (2022/05/16)

- Set default number of ratio panels to 0 [#17](https://github.com/umami-hep/puma/pull/17)
- Adding uncertainties to Roc main plots and improving dummy data generator [#14](https://github.com/umami-hep/puma/pull/14)
- Fixing all PEP8 related variable namings [#13](https://github.com/umami-hep/puma/pull/13)
- Adding FractionScan to puma [#8](https://github.com/umami-hep/puma/pull/8)
