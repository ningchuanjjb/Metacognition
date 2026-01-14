code changes

* ...\default_ops.py
> line 81: add `'PCAdenoiseThres': 0.5,  # PCA denoise threshold`
> line 85: add `'smooth_masks': False,  # whether or not to run smooth_masks (only for sourcery mode)`

* ...\io\binary.py
> line 498: change `500` to `100*bin_size`

* ...\detection\detect.py
> line 35: change `500` to `100*bin_size`
> line 127: change `0.5` to `ops['PCAdenoiseThres']`

* ...\detection\denoise.py
> line 25: add `n_oversamples=10` (`10` could be changed)

* ...\detection\metrics.py
> line 87: change `0.5` to `ops['PCAdenoiseThres']`

* ...\detection\utils.py
> line 203: change `ix+batch_size` to `min(ix+batch_size+1, nbins)`
