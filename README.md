# WaveOpticsBrdf

## Description

This code implements the key ideas of the SIGGRAPH 2018 paper:

**Rendering Specular Microgeometry with Wave Optics**, by

Ling-Qi Yan, Miloš Hašan, Bruce Walter, Steve Marschner and Ravi Ramamoorthi.

It generates a BRDF image, given:

(1) A heightfield, 

(2) A query (center, size) as pixel footprint or coherence area,

(3) An incident direction, and 

(4) A diffraction BRDF model (for wave optics only).

It supports both geometric optics and wave optics. For geometric optics, it uses 
the binning method. For wave optics, you can switch between single and multiple 
wavelengths.

## Usage:

In order to build, you need the eigen library and the OpenEXR library. Both can 
be installed using a package manager such as apt in Ubuntu.

To compile for Linux, use 
```makefile
make -f makefile.linux
```

For macOS, use 
```makefile
make -f makefile.osx
```
Note that, if you are using macOS, by default the clang compiler does not 
have OpenMP support. To enable it, you need to install the libomp library ("brew 
install libomp" if you use Homebrew).

To generate a BRDF image, use the genBrdf command. You'll find the examples 
below useful. Also, be sure to read the command line options in genBrdf.cpp.

We provide three different kinds of heightfields (download separately): 
[Isotropic](https://www.dropbox.com/s/siepjp35pfw218i/isotropic.exr?dl=1)(1.6 MB), 
[Scratched](https://www.dropbox.com/s/p3mm6ws2o18kh3u/scratched.exr?dl=1)(4.3 MB), and  
[Brushed](https://www.dropbox.com/s/vykp2bravkp3tzv/brushed_8K.exr?dl=1)(99.3 MB). 
Note that, to run the examples below, you need to download these heightfields and put them in a *heightfields* folder.

To view the generated BRDF image in EXR format, we recommend the [tev EXR 
viewer](https://github.com/Tom94/tev).

## Examples:

(1) A colored wave optics BRDF generaged from an isotropic heightfield, queried 
at (512, 512) with a size of 10 as one standard deviation. The incident light 
comes from (-1, -1, 1). The output has a resolution of 128 by 128.
```
./genBrdf -m Wave -d OHS -i heightfields/isotropic.exr -o isotropicBrdf.exr -x 512 -y 512 -p 10.0 -w 1.0 -s -1.0 -t -1.0 -r 128
```
(2) A grayscale wave optics BRDF generaged from a scratched isotropic 
heightfield, where each texel has a length of 0.5 microns. The incident light 
comes from top to down by default. The wave length is 0.5 microns (500 nm). The diffraction model is ROHS.
```
./genBrdf -m Wave -d ROHS -i heightfields/scratched.exr -o scratchedBrdf.exr -w 0.5 -x 8512 -y 422 -p 10.0 -l 0.5 -r 128
```
(3) A geometric optics BRDF generated from a large brushed heightfield, with 
five million samples.
```
./genBrdf -m Geom -i heightfields/brushed_8K.exr -o brushed_geom_brdf.exr -r 256 -x 20000 -y 10000 -p 20.0 -n 5000000
```
