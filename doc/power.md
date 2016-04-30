Power
=====

<span class="kcl">class</span> **PowerSpectrum**

```cpp
PowerSpectrum {
  int n;
  double *k, *P
}
```




*ps* = *PowerSpectrum*(<span class="kt">char</span> *filename[]*)

> Reads tabulated power spectrum file k, P(k). Throws PowerFileError
> if unable to open the file

<span class="kt">double</span> *sigma_R* = ps&rarr;compute_sigma(<span class="kt">double</span> *R*)

> Computes rms fluctuation smoothed on length scale *R*

