# Toy Retrieval

Simple Jupyter Notebook toy retrieval where we generate synthetic radiances and then try to retrieve them.

## Dependencies

* absco_lookup.py: code to read the absorption coefficient file (absco.h5)
* find_nearest.py: simple function to find the index of the nearest value in an array
* mie.py: Bohren and Huffman Mie scattering theory
* retrieval.py: contains the majority of the forward model and retrieval code

## Setup

* You'll need to make sure that `absco_table` in `site_settings.ini` points to the file `absco.h5`, wherever that resides on your file system
* You can also change the band ranges, spectral resolutions, geometry, SNR, etc. within settings.py

## Executing program

* In a shell:

```
jupyter notebook toy_retrieval.ipynb
```

* If you change settings.py, you may need to restart your ipynb kernel for it to see those changes
* You can also run the Jupyter Notebook in a GUI such as Visual Studio

