# popsicleR 0.2.1

## Seurat 5 Compatibility

* Added backward-compatible support for Seurat 5.x
* Package now automatically detects Seurat version and uses appropriate API:
  - Seurat < 5: Uses `slot` parameter in `GetAssayData()` and `DoHeatmap()`
  - Seurat >= 5: Uses `layer` parameter in `GetAssayData()` and `DoHeatmap()`
* All existing functionality remains unchanged for users with Seurat < 5
* No changes required to user code - compatibility is handled internally

## Internal Changes

* Added internal compatibility wrapper functions:
  - `GetAssayData_compat()`: Version-aware wrapper for GetAssayData
  - `DoHeatmap_compat()`: Version-aware wrapper for DoHeatmap
* Updated all internal calls to use compatibility wrappers
