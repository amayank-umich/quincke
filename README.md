# Quincke Plugin for HOOMD-blue

This plugin is built using the example-plugin directory in HOOMD-blue. Check out (https://github.com/glotzerlab/hoomd-blue/tree/master/example_plugin). The plugin in itself is a mini makeshift HOOMD-blue program.

[CMakeLists.txt](CMakeLists.txt) and [FindHOOMD.cmake](FindHOOMD.cmake) connects and configures this plugin for hoomd installed in your system.
(quincke) is the main module. It is analogous to, for instance, [md](md) module of hoomd (https://github.com/glotzerlab/hoomd-blue/tree/master/hoomd/md)

In [quincke](quincke), [compute.py](compute.py) is the python interface of the plugin that interacts with the C++ files and cuda files. 

### Note
Currently, this plugin is not developed for GPU. Files .cu and .cuh are more like placeholders.

## Installation
Clone the plugin in some folder
`git clone https://github.com/amayank-umich/quincke.git` Then,
```bash
cd quincke
mkdir build
cd build
cmake ../
make -j4
make install
```
Now you can `import hoomd.quincke` in your simulation script to use the plugin

## Example script
[run.py](run.py) is the example script that uses

[init.gsd](init.gsd) to initialize the system producing

[dump.gsd](dump.gsd) as the output