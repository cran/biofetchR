# biofetchR 0.1.2

## Bug fixes

* Standardised failure handling for explicitly requested terrestrial and
  freshwater vector overlays. TEOW, FEOW, HydroLAKES, HydroRIVERS,
  HydroBASINS, Natural Earth urban and administrative layers, RESOLVE
  ecoregions, WDPA and Ramsar resources now use the pipeline's common
  auditable overlay-loading pathway.

* When `strict_overlay_loading = TRUE`, failures while loading, validating or
  preparing these requested overlays are recorded in `gbif_summary.csv` and
  stop the pipeline before avoidable GBIF download submission.

* Existing fail-soft behaviour is retained when
  `strict_overlay_loading = FALSE`: overlay failures are recorded explicitly
  while downstream processing may continue where possible.

* Added regression tests covering strict and fail-soft overlay-loading
  behaviour, including verification that strict failures occur before the GBIF
  downloader is reached.

* Updated the package-managed Standardising and Integrating Alien Species
  (SInAS) default from dataset 3.1.1 to dataset 3.2 (Zenodo record 21933976;
  DOI 10.5281/zenodo.21933976), while retaining SInAS 3.1.1 as an explicitly
  selectable legacy release.

* Updated SInAS resource resolution and native-origin workflows for the 3.2
  file layout, including `SInAS_3.2.csv`, `All_Config_Files_v3.2.zip` and
  `All_Output_Files_v3.2.zip`.

* Fixed command-line curl downloads on Windows when cache or destination paths
  contain spaces or standalone hyphens, and corrected success detection when
  `system2()` captures curl output in quiet mode.
