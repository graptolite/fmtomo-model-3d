Generate a 3D Wavefront/.obj model of a velocity anomaly isosurface as produced by an FMTOMO run (i.e. with grid files in FMTOMO format). Also options to automatically render the result looking from the SE, SW, NE and NW.

A wireframe map will be placed at the top and bottom of the model volume.

# Dependencies

Linux system with `gmt`, `inskcape` and `blender` (Blender 4.3.2) installed. `fmtomo` working directory using default velocity grid filenames (`vgrids.in` for final, `vgridsref.in` for reference, `vgridstrue.in` for true if it's a recovery test).

Python packages: `json`, `numpy`, `os`, `pandas`, `re`, `scipy`, `shutil`, `skimage`, `subprocess`, `sys`, `tkinter`

Python packages within Blender (e.g., via https://blender.stackexchange.com/a/122337 (where the python executable might instead be python3 or python3.11 etc.) or https://stackoverflow.com/a/68964354) only if `./scripting/script.py` is to be run: `pyyaml`.

# Basic Usage
After running `gui.py`, a Tkinter GUI window will pop up with various options.

If this is the first time running the GUI, or the FMTOMO result of interest has changed, then the desired FMTOMO files (vgrids.in, vgridsref.in, propgrid.in and vgridstrue.in if present) must be copied over. This is done by inputting the filepath to the folder that contains these FMTOMO files (i.e. the FMTOMO working directory) into the "Input Filepath" field, and the pressing the "Load files" button. This moves the relevant files into a temporary (`./tmp`) folder. Unless another filepath is loaded, the files in there will remain the same and so 3D loading runs on the same FMTOMO model when the GUI is closed and reopened.

Once FMTOMO files are loaded, the options for the 3D modelling can be set:
- "dv isosurface (km/s)": the relative velocity isosurface level (in km/s) to be plotted in 3D. Can be positive or negative.
- "material": corresponding material to assign to the relative velocity isosurface specified above.
- "Blender downscale (km to 1 m)": the scale-down from the FMTOMO model domain to the Blender model domain, which can affect how easy it to move around the Blender model domain (the Blender model domain shouldn't span many 10s of m or more).
- "Render": whether to automatically render the result or not. When checked, this will render the result looking from the SW, SE, NE, and NW, saving render results to "./tmp/".
- "Blender GUI": whether to open the Blender GUI or not. Can be unchecked to avoid opening the Blender GUI, which is most useful when render is checked.

Once all the options are set as desired, the button "Load 3D" can be pressed to create and load the 3D model via Blender.

# More Involved Usage
Different parameters can be changed for different types of plot.

## `blender-plot.sh`
As long as the options `-JM${width}c -R$bounds -Bf --MAP_FRAME_TYPE=plain --MAP_TICK_LENGTH=0` and the map file output name are retained, different `gmt coast` calls can be used e.g. rivers and/or lakes can be added in the call, and additional plotting functions can be called.

## `materials.py`
The wireframe map colours can be set in `materials.py` by changing `upper_map_rgb` and `lower_map_rgb` tuples to modify the respective maps. By default, the upper map is green and the lower map is red. The RGB tuple is fractional (i.e. (1,1,1) = 100% R, 100% G, 100% B).

The render background can be set using the `background` RGBA fractional tuple.

Some materials for isosurfaces are defined already in the `materials` variable. More can be defined in the format `{<name> : dict(rgb=<RGB>,alpha=<alpha [optional]>)}` where the RGB tuple is not fractional (i.e. (255,255,255) = 100% R, 100% G, 100% B).

# Script Usage
The script (`python script.py` in the `./scripting` folder) can be used to investigate tomography models in the context of subduction zones as it also handles the loading of Slab depth models from Slab2.0 (Hayes 2018; https://www.usgs.gov/data/slab2-a-comprehensive-subduction-zone-geometry-model). Depths can also be plotted. The controls for this script are found in `./scripting/config.yml`. Copy over the velocity grid files into this directory to use.
