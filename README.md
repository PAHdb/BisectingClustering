# BisectingClustering

This is an implementation of a Divisive Bisecting Clustering algorithm.

This is the IDL implementation by Dr. Christiaan Boersma (Christiaan.Boersma@nasa.gov) and used in the publications:

--- Boersma, C., Bregman, J., Allamandola, L.J., "The Charge State of Polycyclic
Aromatic Hydrocarbons Across Reflection Nebulae: PAH Charge Balance and
Calibration", 2016, ApJ, 832, 51 [https://doi.org/10.3847/0004-637X/832/1/51](https://doi.org/10.3847/0004-637X/832/1/51) ---

--- Boersma, C., Bregman, J., Allamandola, L.J., "Properties of Polycyclic
Aromatic Hydrocarbons in the Northwest Photon Dominated Region of NGC 7023. II.
Traditional PAH Analysis Using k-means as a Visualization Tool", 2014, ApJ, 795,
110 [https://doi.org/10.1088/0004-637X/795/2/110](https://doi.org/10.1088/0004-637X/795/2/110) ---

When using this code, please refer to:

--- Boersma, C., Bregman, J., Allamandola, L.J., "The Charge State of Polycyclic
Aromatic Hydrocarbons Across Reflection Nebulae: PAH Charge Balance and
Calibration", 2016, ApJ, 832, 51 [https://doi.org/10.3847/0004-637X/832/1/51](https://doi.org/10.3847/0004-637X/832/1/51) ---

Note that the sklearn package offers a similar, though not identical, Python implementation of a Bisecting K-Means algorithm in its cluster module ([sklearn.cluster.BisectingKMeans](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.BisectingKMeans.html#sklearn.cluster.BisectingKMeans)).

## Inputs

1. data: the data to cluster (samples x features).
2. img2: bisection_max: maximum number of bisections. The expected number of bins is 2 x bisection_max.
3. min_samples: minimum number of samples needed to be considered a bin.

## Outputs

1. bin_labels: the label for each sample it belongs to (samples).

## Basic Usage

Given 2 test data formatted as samples x features

```(IDL)
 bin_labels = BISECTINGCLUSTERING(data, 2, 10)
```

## Advanced Usage

A spectral cube (n, m, w). For example

```(IDL)
bin_labels = REFORM( $
               BISECTINGCLUSTERING( $
                 REFORM(cube, n * m, w), 2, 10), n, m)
```

## Authors

* **Christiaan Boersma** - *Initial work* - [PAHdb](https://github.com/PAHdb)

## License

This project is licensed under the BSD 3-Clause License - see the
[LICENSE](LICENSE) file for details

## Acknowledgments

* The NASA Ames PAH IR Spectroscopic Database Team
* The Astrophysics & Astrochemistry Laboratory at NASA Ames Research
  Center - [www.astrochemistry.org](http://www.astrochemistry.org)
