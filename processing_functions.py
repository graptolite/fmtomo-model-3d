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

import numpy as np
import skimage
import pandas as pd
from scipy.ndimage import zoom
import os
import json
import re
import subprocess

r_earth = 6371

class Material:
    def __init__(self,name,rgb=[255,255,255],alpha=1):
        self.name = name
        self.rgb = np.array(rgb).astype(int)
        self.alpha = alpha
    def hex(self):
        return "#%02x%02x%02x" % tuple(self.rgb)
    def mtl(self):
        color_string = " ".join((self.rgb/255).astype(str))
        spec = ["Ka %s" % color_string,"Kd %s" % color_string]
        if self.alpha != 1:
            spec.extend(["Tr %.2f" % (1-self.alpha),"d %.2f" % self.alpha])
        return spec

class BlenderGeoDomain():
    def __init__(self,lat0,lat_range,lon0,lon_range,z0,z_range,n_lon,n_lat,n_z,blender_downscale):
        print("Every 1 m in blender is %.2f km in real life" % blender_downscale)
        self.lat0 = lat0
        self.lat_range = lat_range
        self.lat1 = lat0 + lat_range
        self.lon0 = lon0
        self.lon_range = lon_range
        self.lon1 = lon0 + lon_range
        self.z0 = z0
        print(self.z0)
        self.z_range = z_range
        self.z1 = z0 + z_range # Positive is downwards
        print(self.z1)
        self.blender_downscale = blender_downscale
        self.km_per_lat_deg,self.km_per_lon_deg = latlon_grid_to_km(lat0)
        self.ew_range = lon_range * self.km_per_lon_deg
        self.ns_range = lat_range * self.km_per_lat_deg
        self.blender_ew_range = self.ew_range/blender_downscale
        self.blender_ns_range = self.ns_range/blender_downscale
        self.blender_z_range = self.z_range/blender_downscale
        print(self.blender_z_range)
        self.scale_ew = self.blender_ew_range/n_lon
        self.scale_ns = self.blender_ns_range/n_lat
        self.scale_z = self.blender_z_range/n_z
        return
    def get_scale(self):
        return (self.scale_ew,self.scale_ns,self.scale_z)
    def get_dimensions(self):
        return (self.blender_ew_range,self.blender_ns_range,self.blender_z_range)
    def get_map_bounds(self):
        return (self.lon0,self.lon1,self.lat0,self.lat1)

class WavefrontObj():
    def __init__(self,dimensions=[1,1,1],mtl_f=None):
        self.cmd_list = []
        self.mtl_list = []
        self.vert_counter = 1
        self.dimensions = np.array(dimensions)
        self.mtl_f = mtl_f
        if mtl_f:
            self.cmd_list.append("mtllib %s" % mtl_f)
        return
    def add_height_map(self,height_matrix,z_interval,translate=[0,0,0],object_name=""):
        # Avoid adding blank height matrices.
        if all(np.isnan(height_matrix).flatten()):
            return
        self.cmd_list.append("o %s" % object_name)
        n_rows = len(height_matrix)
        const_cols = len(set([len(row) for row in height_matrix])) == 1
        if const_cols:
            n_cols = len(height_matrix[0])
        else:
            print("Error: non-constant number of columns for the rows of the matrix")
        x_interval = self.dimensions[0]/n_rows
        y_interval = self.dimensions[1]/n_cols
        for i,row in enumerate(height_matrix):
            for j,z in enumerate(row):
                coords = [i*x_interval+translate[0],z*z_interval+translate[2],-j*y_interval+translate[1]]
                self.cmd_list.append("v %.2f %.2f %.2f" % tuple(coords))
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
        verts,faces,_,_ = skimage.measure.marching_cubes(grid,isosurface_level,step_size=1)
        self.cmd_list.append("o %.2f %s" % (isosurface_level,name_append))
        # Grid dimensions in order of Z, NS, EW -> EW, NS, Z
        nx,ny,nz = grid.shape[::-1]
        # In order of EW, NS, Z
        scale = self.dimensions/np.array([nx,ny,nz])
        # In order of Z, EW, NS
        scale = np.array([scale[2],scale[0],scale[1]])
        for v in verts:
            # Scale in order of depth, lon, lat
            vx,vy,vz = tuple(v*scale)
            # Fix the order of axes.
            self.cmd_list.append("v %.2f %.2f %.2f" % (vz+translate[0],vx+translate[2],-vy+translate[1]))
        for f in faces:
            self.cmd_list.append("f %u %u %u" % tuple(f+self.vert_counter))
        self.vert_counter += len(verts)
        return
    def add_material(self,mtl):
        self.mtl_list.append("newmtl %s" % mtl.name)
        self.mtl_list += mtl.mtl()
        return
    def set_material(self,mtl):
        self.cmd_list += ["usemtl %s" % mtl.name]
        return
    def write_obj(self,obj_file="obj.obj"):
        with open(obj_file,"w") as outfile:
            outfile.write("\n".join(self.cmd_list))
        if self.mtl_f:
            with open(self.mtl_f,"w") as outfile:
                outfile.write("\n".join(self.mtl_list))
        return

def haversine_dist(r,theta):
    ''' Compute distance between two points separated by an angular distance (`theta`) on a sphere of radius (`r`) using the haversine method.

    r     | <Numerical> | radius of the sphere on which the points lie.
    theta | <Numerical> | angular distance between the two points.

    Returns: <Numerical> | distance between the two points in the same units as `r`.
    '''
    # Compute haversine of the angular separation.
    hav = (1-np.cos(theta))/2
    # Compute and return the distance on a sphere given hav(angular separation).
    return 2 * r * np.arcsin(np.sqrt(hav))

def latlon_grid_to_km(at_lat):
    ''' Compute the cartesian grid spacing (in km per degree) for a given latitude.

    at_lat | <Numerical> | Latitude to compute the spacing at.

    Returns: <Numerical>, <Numerical> | N-S, E-W grid spacing in km per degree.
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
    vs_inv = pd.read_csv(vgrids_f,sep=r"\s+",names=["v","u_v"],skiprows=4)
    vs_ref = pd.read_csv(vgridsref_f,sep=r"\s+",names=["v","u_v"],skiprows=4)
    with open(vgrids_f) as infile:
        shape = infile.read().split("\n")[1].split(" ")
    shape = [int(s) for s in shape if s.strip()]
    vs = pd.DataFrame({"v":vs_inv["v"]-vs_ref["v"]})
    vels = vs.to_numpy().flatten()
    vs_3d = vs.to_numpy().reshape(shape)
    return vs_3d

def is_checker_test():
    try:
        with open("vgridstrue.in") as infile:
            vtrue = infile.read()
        with open("vgridsref.in") as infile:
            vref = infile.read()
        return vref != vtrue
    except FileNotFoundError:
        return False

def split_vgrids(vgrid):
    with open(vgrid) as infile:
        vgrid_data = infile.read().split("\n")
        header,data = vgrid_data[0],vgrid_data[1:]
    if re.split(r"\s+",header.strip())[0].strip() != "1":
        print("Multilayer",vgrid)
        vgrid_splits = [i for i,l in enumerate(data) if "." not in l and l.strip()]
        vgrid_splits.append(len(data))
        files = []
        for i,split in enumerate(vgrid_splits[:-1]):
            out = "      1      1\n"
            out += "\n".join(data[split:vgrid_splits[i+1]])
            outpath = "%s-%u.in" % (vgrid,i)
            with open(outpath,"w") as outfile:
                outfile.write(out.strip())
            files.append(outpath)
    else:
        print("Monolayer",vgrid)
        files = [vgrid]
    return files

def process_dv(vgrids_f,domain,obj_f,isosurface_specs,name_append="",vgrid_ref="vgridsref.in",upsample_factor=1):
    obj = WavefrontObj(domain.get_dimensions(),mtl_f=obj_f + "-mtl.mtl")
    dv = relative_vgrid(vgrids_f,vgrid_ref)
    if upsample_factor > 1:
        dv = zoom(dv,upsample_factor,order=2)
    for level,mtl in isosurface_specs.items():
        obj.add_material(mtl)
        obj.set_material(mtl)
        try:
            obj.add_isosurface(dv,level,translate=[0,0,-(domain.z0/domain.blender_downscale)],name_append="km s-1 %s" % name_append)
        except ValueError:
            pass
    obj.write_obj(obj_f)
    return dv

def parse_in_line(l):
    return [float(x) for x in l.split(" ") if x.strip()]

def construct_blender_domain(blender_downscale,vertical_exaggeration=1):
    with open("vgridsref.in") as infile:
        vdata = infile.read().split("\n")
    n_vz,n_vlat,n_vlon = parse_in_line(vdata[1])
    with open("propgrid.in") as infile:
        pdata = infile.read().split("\n")
    n_pz,n_plat,n_plon = parse_in_line(pdata[0])
    d_pz,d_plat,d_plon = parse_in_line(pdata[1])
    pz0,plat0,plon0 = parse_in_line(pdata[2])
    plat_range = n_plat * d_plat
    plon_range = n_plon * d_plon
    pz_range = n_pz * d_pz
    return BlenderGeoDomain(plat0,plat_range,plon0,plon_range,-pz0,pz_range,n_plat,n_plon,n_pz,blender_downscale)

def make_maps(blender_domain):
    # This is just to avoid too large a document dimension to be created by GMT.
    map_downscale = 100
    if not os.path.exists("blender-plot.sh"):
        with open("blender-plot.sh","w") as outfile:
            outfile.write("""bounds=$1
    width=$2

    gmt begin map pdf
    gmt coast -W1p,black -A100 -JM${width}c -R$bounds -Bf --MAP_FRAME_TYPE=plain --MAP_TICK_LENGTH=0
    gmt end

    inkscape --without-gui --file=map.pdf --export-plain-svg=map.svg""")
        subprocess.call(["bash","blender-plot.sh","%.2f/%.2f/%.2f/%.2f" % blender_domain.get_map_bounds(),"%.2f" % (blender_domain.ns_range/map_downscale)])
    return map_downscale

def exec_3d(blender_downscale,isosurface_specs,render,open_gui):
    os.chdir("tmp")
    fs_obj = [f for f in os.listdir() if f.endswith(".obj")]
    if len(fs_obj) == 0:
        blender_domain = construct_blender_domain(blender_downscale)
        vgrids_fs = split_vgrids("vgrids.in")
        vgridsref_fs = split_vgrids("vgridsref.in")
        print("Checker",is_checker_test())
        if is_checker_test():
            vgridstrue_fs = split_vgrids("vgridstrue.in")
        for grid_f,ref_f in zip(vgrids_fs,vgridsref_fs):
            dv_rec = process_dv(grid_f,blender_domain,"%s.obj" % grid_f,isosurface_specs,name_append="recovered",vgrid_ref=ref_f,upsample_factor=5)
        if is_checker_test():
            for true_f,ref_f in zip(vgridstrue_fs,vgridsref_fs):
                dv_true = process_dv(true_f,blender_domain,"%s.obj" % true_f,isosurface_specs,name_append="true",vgrid_ref=ref_f,upsample_factor=5)
        map_downscale = make_maps(blender_domain)
        with open("tmp.txt","w") as outfile:
            outfile.write(" ".join([str(map_downscale),str(blender_domain.blender_z_range),
                                    str(blender_domain.blender_ew_range),str(blender_domain.blender_ns_range),
                                    str(int(render)),str(int(is_checker_test()))]))
    os.chdir("../")
    cmd_list = ["blender","-P","blender_load_tomo.py"]
    if not open_gui:
        cmd_list.append("--background")
    subprocess.call(cmd_list)
