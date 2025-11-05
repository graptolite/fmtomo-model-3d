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

# Utility functions and classes.

import numpy as np
import skimage
import pandas as pd
from scipy.ndimage import zoom
import os
import json
import re
import subprocess
import shutil

r_earth = 6371

class Material:
    ''' Data class to handle the conversion of material with specified id, RGB color and transparency.

    name  | <str>            | id for the material, which must be unique.
    RGB   | <list> [<int>*3] | color of the material specified in using RGB where each channel has a value between 0 and 255.
    alpha | <float>          | transparency of the material specified as a value between 0 and 1.
    '''
    def __init__(self,name,rgb=[255,255,255],alpha=1):
        self.name = name
        self.rgb = np.array(rgb).astype(int)
        self.alpha = alpha
    def hex(self):
        ''' Convert material's RGB (with which it was initialized) to a hex code.

        Returns: <str> | hex code of material's RGB color.
        '''
        return "#%02x%02x%02x" % tuple(self.rgb)
    def mtl(self):
        ''' Convert material class specification into .mtl format.

        Returns: <list> [<str>] | list containing lines of material's .mtl formatted specification.
        '''
        color_string = " ".join((self.rgb/255).astype(str))
        spec = ["Ka %s" % color_string,"Kd %s" % color_string]
        if self.alpha != 1:
            spec.extend(["Tr %.2f" % (1-self.alpha),"d %.2f" % self.alpha])
        return spec

class BlenderGeoDomain():
    ''' Wrapper class for the storage and easy retrieval of parameters describing the spatial relation between the tomography model's velocity grid (in earth lat-lon coodinates) and the blender model (downscaled, free-floating cartesian coordinates.

    lat0              | <float> | minimum (southernmost) latitude in the tomography model.
    lat1              | <float> | maximum (northernmost) latitude in the tomography model.
    lon0              | <float> | minimum (westernmost) longitude in the tomography model.
    lon1              | <float> | maximum (easternmost) longitude in the tomography model.
    z0                | <float> | minimum depth in the tomography model (i.e. the top of the model).
    z1                | <float> | maximum depth in the tomography model (i.e. the bottom of the model).
    n_lon             | <int>   | number of longitude nodes in the tomography model
    n_lat             | <int>   | number of latitude nodes in the tomography model
    n_z               | <int>   | number of depth nodes in the tomography model
    blender_downscale | <float> | downscale to apply to the tomography model domain (after conversion into cartesian with km units) to get into a cartesian blender domain with m units. This can just be set to 100 in most cases.
    '''
    def __init__(self,lat0,lat1,lon0,lon1,z0,z1,n_lon,n_lat,n_z,blender_downscale):
        # Declare the downscaling applied for the active model.
        print("Every 1 m in blender is %.2f km in real life" % blender_downscale)
        self.lat0 = lat0
        self.lat1 = lat1
        self.lon0 = lon0
        self.lon1 = lon1
        self.lat_range = lat1 - lat0
        self.lon_range = lon1 - lon0
        self.z0 = z0
        self.z1 = z1
        self.z_range = z1 - z0 # Positive is downwards
        self.blender_downscale = blender_downscale
        # Compute latitude and longitude cartesian spacing for the model (assuming the model is small enough that it's uniform throughout the model domain)
        self.km_per_lat_deg,self.km_per_lon_deg = latlon_grid_to_km((self.lat0+self.lat1)/2)
        # Compute the EW (lon) and NS (lat) extent of the model in km.
        self.ew_range = self.lon_range * self.km_per_lon_deg
        self.ns_range = self.lat_range * self.km_per_lat_deg
        # Compute the Blender model extent after applying downscaling.
        self.blender_ew_range = self.ew_range/self.blender_downscale
        self.blender_ns_range = self.ns_range/self.blender_downscale
        self.blender_z_range = self.z_range/self.blender_downscale
        # Declare the depth range in Blender units (m)
        print("Blender depth range:",self.blender_z_range,"Model depth min max:",self.z0,self.z1)
        # Compute node intervals in Blender units.
        self.scale_ew = self.blender_ew_range/n_lon
        self.scale_ns = self.blender_ns_range/n_lat
        self.scale_z = self.blender_z_range/n_z
        return
    def get_scale(self):
        ''' Retrieve velocity grid node intervals in Blender units.

        Returns: <tuple> (<float>) | EW interval, NS interval, Z interval
        '''
        return (self.scale_ew,self.scale_ns,self.scale_z)
    def get_dimensions(self):
        ''' Retrieve dimensions of the Blender model.

        Returns: <tuple> (<float>) | EW extent, NS extend, Z extent
        '''
        return (self.blender_ew_range,self.blender_ns_range,self.blender_z_range)
    def get_map_bounds(self):
        ''' Retrieve extreme latitude and longitude positions of the velocity grid.

        Returns: <tuple> (<float>) | minimum lon, maximum lon, minimum lat, maximum lat.
        '''
        return (self.lon0,self.lon1,self.lat0,self.lat1)

class WavefrontObj():
    ''' Handle Wavefront-format (.obj) file generation procedures for a 3D model.

    dimensions | <list> [<float>] | extent of the 3D model in x,y,z directions.
    mtl_f      | <str>            | .mtl filename to store material specifications in if materials are required.
    '''
    def __init__(self,dimensions=[1,1,1],mtl_f=None):
        # Initialize obj (3D model) command list.
        self.cmd_list = []
        # Initialize mtl (materials) command list.
        self.mtl_list = []
        # Initialize vertices counter.
        self.vert_counter = 1
        self.dimensions = np.array(dimensions)
        self.mtl_f = mtl_f
        # Link mtl file to the obj commands list if required.
        if mtl_f:
            self.cmd_list.append("mtllib %s" % mtl_f)
        return
    def add_height_map(self,height_matrix,z_interval,translate=[0,0,0],object_name=""):
        ''' Add a height map surface spanning the entire 3D model's map-view extent to the obj.

        height_matrix | <np.array>       | 2D array of heights spanning the 3D model's map extent (and with the same orientation).
        z_interval    | <float>          | vertical units of the heightmap relative to blender's 1 unit = 1 m conversion.
        translate     | <list> [<float>] | x,y,z offset of the heightmap.
        object_name   | <str>            | name of the heightmap display in the Blender layers panel.

        Returns: None
        '''
        # Avoid adding blank height matrices.
        if all(np.isnan(height_matrix).flatten()):
            return
        # Declare the object's name in the command list.
        self.cmd_list.append("o %s" % object_name)
        # Pre-check that the height matrix is rectangular.
        n_rows = len(height_matrix)
        const_cols = len(set([len(row) for row in height_matrix])) == 1
        if const_cols:
            n_cols = len(height_matrix[0])
        else:
            print("Error: non-constant number of columns for the rows of the matrix")
        # Compute the x and y intervals in the obj model's domain.
        x_interval = self.dimensions[0]/n_rows
        y_interval = self.dimensions[1]/n_cols
        # Add vertex commands.
        for i,row in enumerate(height_matrix):
            for j,z in enumerate(row):
                coords = [i*x_interval+translate[0],z*z_interval+translate[2],-j*y_interval+translate[1]]
                self.cmd_list.append("v %.2f %.2f %.2f" % tuple(coords))
        # Add face/surface commands.
        for i in range(n_rows-1):
            for j in range(n_cols-1):
                near_nans = [np.isnan(height_matrix[i,j]),
                             np.isnan(height_matrix[i+1,j]),
                             np.isnan(height_matrix[i+1,j+1]),
                             np.isnan(height_matrix[i+1,j-1]),
                             np.isnan(height_matrix[i-1,j]),
                             np.isnan(height_matrix[i-1,j+1]),
                             np.isnan(height_matrix[i-1,j-1]),
                             np.isnan(height_matrix[i,j+1]),
                             np.isnan(height_matrix[i,j-1]),
                             ]
                # Avoid adding a face if self or neighbouring points (in square ring) have nan values).
                if not any(near_nans):
                    c = i*n_cols+self.vert_counter + j
                    face = [c,c+1,c+n_cols,c+1+n_cols]
                    self.cmd_list.append("f %u %u %u" % tuple(face[:-1]))
                    self.cmd_list.append("f %u %u %u" % tuple(face[1:]))
        self.vert_counter += height_matrix.shape[0]*height_matrix.shape[1]
        return
    def add_isosurface(self,grid,isosurface_level,translate=[0,0,0],name_append=""):
        ''' Add an isovalue surface from a 3D grid of values spanning the 3D model's extent to the obj.

        grid             | <np.array>       | 3D array of values in the order [z,y,x]
        isosurface_level | <float>          | value to interpolate an isosurface of in the 3D grid of values.
        translate        | <list> [<float>] | x,y,z offset of the heightmap.
        name_append      | <str>            | string to append to the isosurface object, which will be named by the isosurface level (e.g. units for the value).
        '''
        # Compute the vertices and faces that describe the isosurface.
        verts,faces,_,_ = skimage.measure.marching_cubes(grid,isosurface_level,step_size=1)
        # Declare the object's name in the command list.
        self.cmd_list.append("o %.2f %s" % (isosurface_level,name_append))
        # Convert grid dimensions in order of Z, NS, EW -> EW, NS, Z
        grid_dims = np.array(grid.shape)[::-1]
        # Get the scaling factors (from isosurface model domain to the 3D model domain) in order of EW, NS, Z (same order as self.dimensions)
        scale = self.dimensions/grid_dims
        # Since the marching cubes algo returns coordinates in a zero indexed node system whereas the model dimension counting system is 1 indexed, need to account for the difference in scaling of the maximum node location.
        index_mismatch_scaling = grid_dims/(grid_dims-1)
        scale *= index_mismatch_scaling
        scale = scale[::-1]
        # Add vertex commands.
        for v in verts:
            # Scale in order of depth, lon, lat
            vz,vy,vx = tuple(v*scale)
            # Fix the order of axes for Blender coordinate system.
            self.cmd_list.append("v %.2f %.2f %.2f" % (vx+translate[0],vz+translate[2],-vy+translate[1]))
        # Add face commands.
        for f in faces:
            self.cmd_list.append("f %u %u %u" % tuple(f+self.vert_counter))
        self.vert_counter += len(verts)
        return
    def add_material(self,mtl):
        ''' Add a material specification to the mtl commands list.

        mtl | <Material> | material object containing material specification.
        '''
        self.mtl_list.append("newmtl %s" % mtl.name)
        self.mtl_list += mtl.mtl()
        return
    def set_material(self,mtl):
        ''' Assign a material to be used for the following obj commands.

        mtl | <Material> | material object containing material specification.
        '''
        self.cmd_list += ["usemtl %s" % mtl.name]
        return
    def write_obj(self,obj_file="obj.obj"):
        ''' Write obj (and mtl if present) commands list(s) to disk.

        obj_file | <str> | file to save the obj commands list to.
        '''
        # Save obj commands list to obj file.
        with open(obj_file,"w") as outfile:
            outfile.write("\n".join(self.cmd_list))
        if self.mtl_f:
            # Save mtl commands list to mtl file if requested when constructing this 3D obj model class.
            with open(self.mtl_f,"w") as outfile:
                outfile.write("\n".join(self.mtl_list))
        return

def haversine_dist(r,theta):
    ''' Compute distance between two points separated by an angular distance (`theta`) on a sphere of radius (`r`) using the haversine method.

    r     | <float> | radius of the sphere on which the points lie.
    theta | <float> | angular distance between the two points.

    Returns: <float> | distance between the two points in the same units as `r`.
    '''
    # Compute haversine of the angular separation.
    hav = (1-np.cos(theta))/2
    # Compute and return the distance on a sphere given hav(angular separation).
    return 2 * r * np.arcsin(np.sqrt(hav))

def latlon_grid_to_km(at_lat):
    ''' Compute the cartesian grid spacing (in km per degree) for a given latitude.

    at_lat | <float> | Latitude (in degrees) to compute the spacing at.

    Returns: <float>, <float> | N-S, E-W grid spacing in km per degree.
    '''
    # Convert one degree to radians.
    rad_d = np.radians(1)
    # Compute the N-S grid spacing (latitudinal grid spacing).
    km_per_lat_deg = haversine_dist(r_earth,rad_d)
    # Compute the E-W grid spacing at the latitude of interest by setting `r` to the radius of the small circle at that line of latitude.
    r_small = r_earth * np.sin(np.radians(90-at_lat))
    km_per_lon_deg = haversine_dist(r_small,rad_d)
    return km_per_lat_deg,km_per_lon_deg

def relative_vgrid(vgrids_f,vgridsref_f):
    ''' Compute relative velocity differences between a "observed" velocity grid and reference velocity grid. These velocity grids will have cushion nodes that get removed throughout this workflow.

    vgrids_f    | <str> | filepath to observed velocity grid in fmtomo format (i.e. structured 1D list of velocities).
    vgridsref_f | <str> | filepath to reference velocity grid in fmtomo format.

    Returns: <np.array> | 3D (grid) array with relative velocity values (observed - reference).
    '''
    # Load velocity grids.
    vs_inv = pd.read_csv(vgrids_f,sep=r"\s+",names=["v","u_v"],skiprows=4)
    vs_ref = pd.read_csv(vgridsref_f,sep=r"\s+",names=["v","u_v"],skiprows=4)
    # Identify the (common) velocity grid shape.
    with open(vgrids_f) as infile:
        shape = infile.read().split("\n")[1].split(" ")
    shape = [int(s) for s in shape if s.strip()]
    # Compute relative velocity at each grid point.
    vs = pd.DataFrame({"v":vs_inv["v"]-vs_ref["v"]})
    # Shape the structured 1D velocity grid list into a 3D grid array.
    vs_3d = vs.to_numpy().reshape(shape)
    # Remove the cushion nodes.
    vs_3d = vs_3d[1:-1,1:-1,1:-1]
    return vs_3d

def is_recovery_test():
    ''' Check whether the fmtomo working directory's velocity grids (with default filenames) reflect a recovery test run or not.

    Returns: <bool> | if the velocity grids represent a recovery test.
    '''
    try:
        with open("vgridstrue.in") as infile:
            vtrue = infile.read()
        with open("vgridsref.in") as infile:
            vref = infile.read()
        # If both a reference and true velocity grid file are present, return the boolean for if they are different (if different, then this is a recovery test).
        return vref != vtrue
    except FileNotFoundError:
        # If no vgridstrue.in (vgridsref.in must be present), then can't be a recovery test.
        return False

def split_vgrids(vgrid):
    ''' Split a fmtomo velocity grid that may contain multiple layers into multiple single-layer velocity grid files.

    vgrid | <str> | path to the (potentially multi-layer) velocity grid.

    Returns: <list> [<str>] | list of files containing exclusively single-layer fmtomo velocity grid files.
    '''
    # Load contents of velocity grid file and identify the header (with the first number in this header being the number of layers).
    with open(vgrid) as infile:
        vgrid_data = infile.read().split("\n")
        header,data = vgrid_data[0],vgrid_data[1:]
    # Check if the number of layers in the velocity grid file is not 1.
    if re.split(r"\s+",header.strip())[0].strip() != "1":
        # Identify the velocity grid as being multi-layered.
        print("Multilayer",vgrid)
        # Identify velocity grid header line locations (which specify the number of grid nodes i.e. contain only integer numbers).
        vgrid_splits = [i for i,l in enumerate(data) if "." not in l and l.strip()]
        vgrid_splits.append(len(data))
        # Initialise list to store filenames for the split single-layer velocity grids.
        files = []
        # Iterate through the header line locations.
        for i,split in enumerate(vgrid_splits[:-1]):
            # Construct single-layer velocity grid specifications.
            out = "      1      1\n"
            out += "\n".join(data[split:vgrid_splits[i+1]])
            # Write single-layer velocity grid specifications to disk with a unique filename.
            outpath = "%s-%u.in" % (vgrid,i)
            with open(outpath,"w") as outfile:
                outfile.write(out.strip())
            # Store the unique filename for the single-layer velocity grid.
            files.append(outpath)
    else:
        # Identify the velocity grid as being single-layered and perform no splitting operations.
        print("Monolayer",vgrid)
        files = [vgrid]
    return files

def process_dv(vgrid_f,domain,obj_f,isosurface_specs,name_append="",vgrid_ref="vgridsref.in",upsample_factor=1):
    ''' Process a velocity grid file in order to plot an isosurface of relative velocity (into an obj 3D model) from it with requested domain mappings.

    vgrid_f          | <str>                       | velocity grid file for input.
    domain           | <BlenderGeoDomain>          | domain properties for the conversion of the velocity grid domain to 3D model domain.
    obj_f            | <str>                       | obj filename for writing the 3D model output.
    isosurface_specs | <dict> {<float>:<Material>} | isosurface plotting specifications in dictionary format: {isosurface level : material to apply to isosurface}.
    name_append      | <str>                       | string to append to the name of each isosurface after the default units (km s-1).
    vgrid_ref        | <str>                       | reference velocity grid file for the computation of relative velocities.
    upsample_factor  | <int>                       | factor by which to upsample the 3D velocity grid.

    Returns: <np.array> | 3D relative velocity grid.
    '''
    # Initialize wavefront object with a material filename supplied based on the name of the object file.
    obj = WavefrontObj(domain.get_dimensions(),mtl_f=obj_f + "-mtl.mtl")
    # Obtain the 3D relative velocity grid.
    dv = relative_vgrid(vgrid_f,vgrid_ref)
    # Upsample this relative velocity grid using 2nd order spline interpolation if requested.
    if upsample_factor > 1:
        dv = zoom(dv,upsample_factor,order=2)
    # Iterate through the requested isosurface levels.
    for level,mtl in isosurface_specs.items():
        # Add and set for use the isosurface level's material.
        obj.add_material(mtl)
        obj.set_material(mtl)
        # Extract and add the isosurface from the relative velocity 3D grid.
        try:
            obj.add_isosurface(dv,level,name_append="km s-1 %s" % name_append)
        except ValueError:
            pass
    # Write the obj (and mtl) files.
    obj.write_obj(obj_f)
    return dv

def parse_in_line(l,t=float):
    ''' Convert a line of space-separated strings to a list of specified types.

    l | <str>  | line of space-separated strings that can be cast to floats.
    t | <type> | type to cast the list items to

    Returns: <list> [<type>] | list of specified types.
    '''
    return [t(x) for x in l.split(" ") if x.strip()]

def construct_blender_domain(blender_downscale,vgrid_f="vgrids.in"):
    ''' Construct BlenderGeoDomain object from the grid file of an fmtomo run.

    blender_downscale | <float> | downscale to apply to the tomography model domain (after conversion into cartesian with km units) to get into a cartesian blender domain with m units. This can just be set to 100 in most cases.
    '''
    # Identify the number of nodes in each direction in the velocity grid.
    with open("vgridsref.in") as infile:
        vdata = infile.read().split("\n")
    # y (radians) -> lat (degrees), x (radians) -> lon (degrees)
    n_z,n_y,n_x = parse_in_line(vdata[1],int)
    dz,dy,dx = parse_in_line(vdata[2],float)
    z0,y0,x0 = parse_in_line(vdata[3],float)
    dy,dx,y0,x0 = (np.degrees(p) for p in [dy,dx,y0,x0])
    get_bounds = lambda n,d,p0 : p0 + d * np.array([1,n-2])
    z0,z1 = get_bounds(n_z,dz,z0)
    z0 -= r_earth
    z1 -= r_earth
    lat0,lat1 = get_bounds(n_y,dy,y0)
    lon0,lon1 = get_bounds(n_x,dx,x0)
    BGD = BlenderGeoDomain(lat0=lat0,lat1=lat1,
                           lon0=lon0,lon1=lon1,
                           z0=z0,z1=z1,
                           n_lon=n_x,n_lat=n_y,n_z=n_z,
                           blender_downscale=blender_downscale)
    return BGD

def make_maps(blender_domain,map_downscale=100):
    ''' Generate vector (svg) map covering the horizontal extent of the 3d velocity model.

    blender_domain | <BlenderGeoDomain> | BlenderGeoDomain object for the velocity model and 3D model domains.
    map_downscale  | <int>              | a strong GMT map (document) downscale factor just to avoid too large a document dimension to be created by GMT (this downscaling will be reversed before plotting in Blender).

    Returns: None
    '''
    # Write the GMT map plotting script if not already present.
    if not os.path.exists("blender-plot.sh"):
        shutil.copy("../blender-plot.sh","blender-plot.sh")
    # Execute the map plotting script to cover the relevant domain.
    # N-S and E-W range are redundant for the Blender workflow (as long as they are strongly downscaled) but make sure that the map can be viewed at roughly the correct aspect ratio outside of this workflow. map_downscale is an arbitrary, large number here and is not used to set dimensions here (only the SVG internal coordinate system).
    subprocess.call(["bash","blender-plot.sh","%.2f/%.2f/%.2f/%.2f" % blender_domain.get_map_bounds(),"%.2f" % (blender_domain.ew_range/map_downscale),"%.2f" % (blender_domain.ns_range/map_downscale)])
    # This filename shouldn't need changing otherwise this workflow will break at the Blender loading stage also.
    map_svg = "map.svg"
    if not os.path.exists(map_svg):
        raise FileNotFoundError("Map could not be generated. Make sure Inkscape is installed to the system path.")
    # Fix the scale of the SVG output to be directly loadable into Blender.
    fix_svg_scale(map_svg,blender_domain,map_downscale)
    return

def fix_svg_scale(svg_f,blender_domain,map_downscale):
    ''' Modify the SVG header of a file to match its size to that specified by the BlenderGeoDomain.

    svg            | <str>              | name of the svg input file.
    blender_domain | <BlenderGeoDomain> | BlenderGeoDomain object for the velocity model and 3D model domains.
    map_downscale  | <int>              | a strong GMT map (document) downscale factor just to avoid too large a document dimension to be created by GMT (this downscaling will be reversed before plotting in Blender).

    Returns: None
    '''
    # Load the contents of the svg file.
    with open(svg_f) as infile:
        svg = infile.read()
    # Read the svg tag properties (referred to as the header here).
    svg_header = svg_header_original = re.search(r"<svg[\S\s]+?>",svg).group(0)
    # Update svg dimensions with the blender model domain dimensions (in mm) after reversible map downscaling.
    m2mm = 1e3
    effective_scale = m2mm/map_downscale
    width = blender_domain.blender_ew_range*effective_scale
    height = blender_domain.blender_ns_range*effective_scale
    svg_header = re.sub("width=\"(.*?)\"","width=\"%.2fmm\"" % width,svg_header)
    svg_header = re.sub("height=\"(.*?)\"","height=\"%.2fmm\"" % height,svg_header)
    # Update SVG file with updated dimensions.
    with open(svg_f,"w") as outfile:
        outfile.write(svg.replace(svg_header_original,svg_header))
    return

def exec_3d(blender_downscale,isosurface_specs,render,open_gui,blender_script="blender_load_tomo.py",force_rerun=False,map_downscale=None):
    ''' Execute preprocessing and Blender script execution in a temp dir containing all the necessary grid files (at least vgrids.in and vgridsref.in plus vgridstrue.in if the fmtomo working directory is for a recovery test) copied over from an fmtomo working dir.

    blender_downscale | <float>                     | downscale to apply to the tomography model domain (after conversion into cartesian with km units) to get into a cartesian blender domain with m units. This can just be set to 100 in most cases.
    isosurface_specs  | <dict> {<float>:<Material>} | isosurface plotting specifications in dictionary format: {isosurface level : material to apply to isosurface}.
    render            | <bool>                      | whether to autorender (looking from SE, SW, NW, NE directions).
    open_gui          | <bool>                      | whether to launch the Blender GUI.
    blender_script    | <str>                       | path to script to load into Blender.
    map_downscale     | <int>                       | a strong GMT map (document) downscale factor just to avoid too large a document dimension to be created by GMT (this downscaling will be reversed before plotting in Blender).
    '''
    # Avoid chdir into tmp if already there (e.g. due to a previous failure to load).
    if os.path.basename(os.path.dirname(os.getcwd())) != "tmp":
        # Move to the temporary directory
        if not os.path.exists("tmp"):
            os.mkdir("tmp")
        os.chdir("tmp")
    # Get a list of obj files in the tmp dir.
    fs_obj = [f for f in os.listdir() if f.endswith(".obj")]
    # Create an fmtomo-Blender domain mapping using the grid files in the temp dir.
    blender_domain = construct_blender_domain(blender_downscale)
    if not map_downscale:
        # If no map downscaling provided, set one. To avoid rounding errors, the map_downscale factor should not be set to too high a value, e.g., target a map width of around 10 cm.
        target_map_width = 10 # cm
        map_downscale = blender_domain.ew_range/target_map_width
    # Only remake obj files if there are no obj files in there already.
    if len(fs_obj) == 0 or force_rerun:
        # Split the velocity grids as necessary
        vgrids_fs = split_vgrids("vgrids.in")
        vgridsref_fs = split_vgrids("vgridsref.in")
        print("Recovery",is_recovery_test())
        if is_recovery_test():
            vgridstrue_fs = split_vgrids("vgridstrue.in")
        # Get isosurfaces for each of the grids.
        for grid_f,ref_f in zip(vgrids_fs,vgridsref_fs):
            dv_rec = process_dv(grid_f,blender_domain,"%s.obj" % grid_f,isosurface_specs,name_append="recovered",vgrid_ref=ref_f,upsample_factor=5)
        if is_recovery_test():
            for true_f,ref_f in zip(vgridstrue_fs,vgridsref_fs):
                dv_true = process_dv(true_f,blender_domain,"%s.obj" % true_f,isosurface_specs,name_append="true",vgrid_ref=ref_f,upsample_factor=5)
        # Create vector GMT map for the fmtomo model.
        make_maps(blender_domain,map_downscale)
    # Dump a bunch of parameters for passing on to the model loading code to be run via Blender.
    with open("tmp.txt","w") as outfile:
        outfile.write(" ".join([str(map_downscale),str(blender_domain.blender_z_range),str(blender_domain.z1/blender_domain.blender_downscale),
                                str(blender_domain.blender_ew_range),str(blender_domain.blender_ns_range),
                                str(int(render)),str(int(is_recovery_test()))]))
    os.chdir("../")
    # Run blender on the processing script.
    exec_blender_with_script(blender_script,open_gui)
    return

def exec_blender_with_script(script,open_gui=True):
    ''' Execute Blender with a Blender Python script, and choose whether to open the Blender GUI or not. This also handles Blenders that are on the console bin path but not installed onto the default bin path.

    script   | <str>  | path to the Blender Python script.
    open_gui | <bool> | whether to open the Blender GUI or not.

    Returns: None
    '''
    # Default list of commands used to run the fmtomo model loading code through Blender.
    cmd_list = ["blender","-P",script]
    # Avoid opening the Blender GUI if requested.
    if not open_gui:
        cmd_list.append("--background")
    try:
        # Try loading under the assumption of a system install of blender.
        subprocess.call(cmd_list)
    except:
        # Otherwise load blender under the assumption of a bash alias.
        subprocess.call(["/bin/bash","-i","-c"] + [" ".join(cmd_list)])
    return
