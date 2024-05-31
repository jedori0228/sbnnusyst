# sbnnusyst

This package enables user to update caf::SRTrueInteraction::wgt to include new cross-section reweights evaluated by [nusystematics](https://github.com/NuSystematics/nusystematics).

# dependencies

## GENIE

```
setup genie v3_04_00 -qe20:prof
setup genie_xsec v3_04_00 -qAR2320i00000:e1000:k250
```

## sbnanaobj

### 1) Using tagged version

```
setup sbnanaobj <version> -q<qaul>
```

### 2) Using local product

As of 31May2024, we don't have sbnanaobj release that includes relevant PRs, so we need a local build.
Below is an example that uses sbnanaobj of v09_20_06_03

#### i) Initial setup

```
mkdir -p sbnanaobj; cd sbnanaobj;
dirName=v09_20_06_03
mkdir ${dirName}; cd ${dirName}
mrb newDev -v v09_72_00 -q e20:prof # sbnanaobj of v09_20_06_03 matches with LarSoft of v09_72_00
source localProducts_larsoft_v09_72_00_e20_prof/setup
cd srcs
mrb g -t v09_20_06_03 sbnanaobj
mrbsetenv
# make local updates, and compile
mrb i -j4
```

#### ii) New shell after initialization

You only need to source the setup file:

```
cd ${mywd}/sbnanaobj/v09_20_06_03/
source localProducts_larsoft_v09_72_00_e20_prof/setup
mrbsetenv
```

## nusystematics

### Installation

```
# #{mywd} is your working area
cd ${mywd} # go to your working directory
mkdir nusystematics; cd nusystematics
git clone git@github.com:NuSystematics/nusystematics.git nusystematics-src
mkdir build; cd build
cmake ../nusystematics-src/
make install
```

### setup script

Whenever you open a new shell, run
```
source ${mywd}/nusystematics/build/Linux/bin/setup.fhicl_cpp_standalone.sh
source ${mywd}/nusystematics/build/Linux/bin/setup.systematicstools.sh
source ${mywd}/nusystematics/build/Linux/bin/setup.nusystematics.sh
```

# build

```
# #{mywd} is your working area
cd ${mywd} # go to your working directory
mkdir sbnnusyst; cd sbnnusyst;
git clone git@github.com:jedori0228/sbnnusyst.git sbnnusyst-src
mkdir build; cd build
cmake ../sbnnusyst-src/
make install
```

# setup script

Whenever you open a new shell, run
```
source ${mywd}/sbnnusyst-src/build/Linux/bin/setup.sbnnusyst.sh
```

# Running UpdateReweight

You need a nusyst configuration fhicl file. An example is located in `example/zexpansion_weighter.ParameterHeader.fcl`.

Then write a txt file that contains the list of input CAF files:
```
$ cat input_cafs.txt
path/to/your/input/caf.root
```

Then run `UpdateReweight` using following command:
```
UpdateReweight -c example/zexpansion_weighter.ParameterHeader.fcl -i input_cafs.txt -o output_flat.caf.root
```



