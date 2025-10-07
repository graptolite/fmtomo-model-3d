#!/usr/bin/env python3

'''
FMTOMO Result 3D Modelling | Generate a 3D Wavefront/.obj model of a velocity anomaly isosurface as produced by an FMTOMO run (i.e. with grid files in FMTOMO format).
    Copyright (C) 2025 Yingbo Li

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

# Example script that also includes the handling of slab plotting.

import numpy as np
import re
import matplotlib as mpl
from scipy.interpolate import LinearNDInterpolator
import shutil
from itertools import combinations
import yaml

# Load configuration file.
with open("config.yml") as infile:
    config = yaml.safe_load(infile)["plot_3d"]

import sys
sys.path.insert(0,"../")
from processing_functions import *

# Load in configuration options.
blender_downscale = config["blender_downscale"]
render = config["render"]
isosurface_specs = {c["dv"]:Material(k,c["RGB"]) for k,c in config["isosurfaces"].items()}
slab2_txt_dir = config["slab2_txt_dir"]

def load_slab_data(domain,slab_f="slabs.json",slab2_txt_dir=slab2_txt_dir):
    ''' Identify slabs from Slab2.0 (with input files in txt format as downloaded and unzipped from Slab2Distribute_Mar2018.tar.gz on https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467) that lie wholly or partially within the Blender model domain and collate them into one data structure (dict) to save to disk.

    domain        | <BlenderGeoDomain> | object containing information on the domain of interest.
    slab_f        | <str>              | json file to save the collated slab data to.
    slab2_txt_dir | <str>              | input folder containing Slab2 data in txt format (with depths in km).

    Returns: <dict> {<str>:<dict>} | collated data in dictionary format {<slab name>:<slab geometry (x,y,z) data in dict format>}.
    '''
    # Identify domain bounds.
    min_lon,max_lon,min_lat,max_lat = domain.get_map_bounds()
    domain_str = "/".join([str(round(x,3)) for x in [min_lon,max_lon,min_lat,max_lat]])
    # Check whether an existing slab datafile exists.
    if os.path.exists(slab_f):
        # Load the collated slabs from disk.
        with open(slab_f) as infile:
            accepted = json.load(infile)
        # Check whether the domain of the existing slab datafile is suitable, and if so return the data and end the function early.
        if "domain" in accepted and accepted["domain"] == domain_str:
            return accepted
    # Initialize dictionary to store accepted slabs within the domain bounds.
    accepted = {"slabs":dict(),
                "domain":domain_str,}
    # Get list of Slab2 depth files.
    fs = [f for f in os.listdir(slab2_txt_dir) if "_dep_" in f]
    # Iterate through the individual Slab2 depth files.
    for f in fs:
        print("Processing slab file",f)
        # Load the slab depth data.
        df = pd.read_csv(os.path.join(slab2_txt_dir,f),names=["lon","lat","dep"])
        # Remove empty entries from the slab depth data.
        df = df[df.dep.notna()]
        # Check whether there are any slab depth entries within the model domain.
        bound_check = (df.lat>min_lat) & (df.lat<max_lat) & (df.lon>min_lon) & (df.lon<max_lon)
        if any(bound_check):
            # Store the slab as accepted if there are.
            accepted["slabs"][f.split("_")[0]] = df.to_dict()
    # Save the collated slabs to disk.
    with open(slab_f,"w") as outfile:
        json.dump(accepted,outfile,indent=2)
    return accepted

def construct_slabs(domain,interp_interval=0.25,obj_file="./tmp/slabs.obj",critical_delta_depth=50,cmap_name="viridis",cmap_norm_add=0,cmap_idx_offset=0):
    ''' Generate obj files for the slabs that have been isolated within the requested domain.

    domain               | <BlenderGeoDomain> | object containing information on the domain of interest.
    interp_interval      | <float>            | interval in degrees to perform linear interpolation in a grid on across the lat-lon domain.
    obj_file             | <str>              | path to save the obj data to.
    critical_delta_depth | <float>            | maximum permissible difference in depth (in km) between neighbouring points (to avoid sharp differences caused arising as interpolation artefacts). Note that this is not well tested at the moment.
    cmap_name            | <str>              | name of the matplotlib colormap to color the slabs using. The assigment of color works as follows (ignoring the scaling to 255 and index offsetting): cmap[(slab_index+cmap_idx_offset)/(N_slabs_total+cmap_norm_add)].
    cmap_norm_add        | <int>              | number to control the cmap application (see above).
    cmap_idx_offset      | <int>              | number to control the cmap application (see above).

    Returns: None
    '''
    # Load the within-domain slab data.
    accepted = load_slab_data(domain)["slabs"]
    # Identify the domain region.
    min_lon,max_lon,min_lat,max_lat = domain.get_map_bounds()
    # Compute the scaled z interval (i.e. distance in the Blender domain per unit in the slab depth domain).
    z_interval = 1/domain.blender_downscale
    # Initialize obj constructor for the slab with a corresponding material file.
    slab_obj = WavefrontObj(domain.get_dimensions(),mtl_f=os.path.join(os.path.dirname(obj_file),"slabs.mtl"))
    # Identify the cmap to apply to the slabs.
    cmap = mpl.colormaps[cmap_name]
    # Compute the norm from the zero-index corrected length and additional component supplied.
    cmap_norm = (len(accepted) - 1) + cmap_norm_add
    # Assign each slab a color for its material.
    for i,slab in enumerate(accepted.keys()):
        mtl = Material(slab,255*np.array(cmap(int(255 * (i+cmap_idx_offset)/cmap_norm))[:3]))
        # Parse the slab depth data into a pandas DataFrame.
        slab_depths = pd.DataFrame(accepted[slab])
        # Perform crop (with 2*interp_interval margin) to avoid interpolation on points close to the margin.
        crop_n = 2
        # Remove marginal depth points.
        slab_depths = slab_depths[(slab_depths.lat>(min_lat-crop_n*interp_interval)) &
                                  (slab_depths.lat<(max_lat+crop_n*interp_interval)) &
                                  (slab_depths.lon>(min_lon-crop_n*interp_interval)) &
                                  (slab_depths.lon<(max_lon+crop_n*interp_interval))]
        # Check if parts of the active slab remain after margin cropping.
        if len(slab_depths):
            # Initialize interpolator with the slab depths and the positions at which those were computed.
            interp = LinearNDInterpolator(list(zip(slab_depths.lon,slab_depths.lat)),slab_depths.dep,fill_value=np.nan)
            # Construct regularly-spaced grid of points at which to perform interpolation.
            ilons,ilats = np.meshgrid(np.arange(min_lon,max_lon,interp_interval),np.arange(min_lat,max_lat,interp_interval))
            # Interpolate slab depths on the interpolation grid.
            z_interp = interp(ilons,ilats).T
            # Subsetting indices for the grid of depths for comparison of neighbouring depth interpolations.
            idx_mapping = {1:[-1,-1],
                           2:[-1,1],
                           3:[1,-1],
                           4:[1,1]}
            # Construct iterable of unique comparisons.
            comparisons = combinations(idx_mapping.keys(),2)
            # Iterate through the unique comparisons.
            for idx_a,idx_b in comparisons:
                # Determine the relevant matrix slicing indices in a not so efficient manner.
                mat_a_slicing = idx_mapping[idx_a]
                mat_b_slicing = idx_mapping[idx_b]
                idxs = []
                for slicing in [mat_a_slicing,mat_b_slicing]:
                    subidxs = []
                    for idx in slicing:
                        # Positive index means remove from the start ([idx:None]). Negative index means remove from the end ([None:idx]). This orders the slice indices the right way round.
                        slice_idx = [idx,None][::idx]
                        subidxs.append(slice_idx)
                    idxs.append(subidxs)
                # Isolate the submatrices for comparison.
                mat_a = z_interp[idxs[0][0][0]:idxs[0][0][1],idxs[0][1][0]:idxs[0][1][1]]
                mat_b = z_interp[idxs[1][0][0]:idxs[1][0][1],idxs[1][1][0]:idxs[1][1][1]]
                # Compute difference between the two comparison submatrices.
                deltas = np.array(mat_a-mat_b)
                # Pad the matrix with zeros to match the size of the full depth interpolation matrix.
                ## Use the slicing to get mat_a as the reference. If the slicing along the axis is 1 (i.e. cutting off the start), then pad the zeros at the start (i.e. do not change the order of the stack); if -1 (cutting off the end), then pad zeros at the end (reverse order of the stack).
                deltas = np.hstack([np.zeros([len(deltas[:,0]),1]),deltas][::mat_a_slicing[0]])
                deltas = np.vstack([np.zeros([1,len(deltas[0,:])]),deltas][::mat_a_slicing[1]])
                # Remove sudden differences in point height greater than the threshold.
                z_interp[abs(deltas)>critical_delta_depth] = np.nan
            # Add material for the active slab.
            slab_obj.add_material(mtl)
            # Set the material for the active slab as the one that was just added.
            slab_obj.set_material(mtl)
            # Add the slab height map at the correct transformed location.
            slab_obj.add_height_map(z_interp,z_interval,object_name=slab,
                                    translate=[0,0,domain.blender_z_range])
    # Save the slabs object to disk.
    slab_obj.write_obj(obj_file)
    return

if os.path.exists("tmp"):
    # Clean files from the tmp dir for fresh reprocessing.
    for f in os.listdir("tmp"):
        os.remove(os.path.join("tmp",f))

# Construct a Blender domain from the vgrids.in file in the current directory.
blender_domain = construct_blender_domain(blender_downscale)
# Attempt to load the required grid files.
for f in ["vgrids.in","vgridsref.in"]:
    shutil.copy(f,"tmp")
# Load the optional grid file if present.
v_true = "vgridstrue.in"
if os.path.exists(v_true):
    shutil.copy(v_true,"tmp")

# Generate the slabs object.
construct_slabs(blender_domain,interp_interval=0.1,cmap_name="magma_r",cmap_norm_add=2,cmap_idx_offset=1)

# If an additional plotting file exists, execute that file with a large degree of map downscaling and copy over the output svg it generates to the tmp dir.
map_downscale = 100
if os.path.exists("blender-plot-frame.sh"):
    subprocess.call(["bash","blender-plot-frame.sh","%.2f/%.2f/%.2f/%.2f" % blender_domain.get_map_bounds(),"%.2f" % (blender_domain.ew_range/map_downscale),"%.2f" % (blender_domain.ns_range/map_downscale)])
    fix_svg_scale("frame.svg",blender_domain,map_downscale)
    os.rename("frame.svg","tmp/frame.svg")

# Run the 3D modelling.
exec_3d(blender_downscale,isosurface_specs,render=render,open_gui=config["open_gui"],blender_script="blender_load_tomo_depth_plots.py",force_rerun=True,map_downscale=map_downscale)
