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

import pandas as pd
import numpy as np
import os
from io import StringIO

xyz_f = "vel_xyz.csv"

out_folder = "generated"
# REQUIRED: even gridding
# REQUIRED UNITS: lat lon in degrees, dep in km with positive being down, dv in kms-1; change column order for parsing the input velocity model as necessary
# LEAVE EMPTY IF COLUMN HEADERS ARE SET IN THE CSV
column_order = ["lat","lon","dep","dv"]
sep = ","

###########

r_earth = 6371

if len(column_order):
    df = pd.read_csv(xyz_f,names=column_order,sep=sep)
else:
    df = pd.read_csv(xyz_f,sep=sep)
n_x = nlon = len(set(df.lon))
n_y = nlat = len(set(df.lat))
n_z = ndep = len(set(df.dep))
x0 = min(df.lon)
y0 = min(df.lat)
z0 = r_earth - max(df.dep)
lats = sorted(set(df.lat))
lons = sorted(set(df.lon))
deps = sorted(set(df.dep))
deltas = lambda l : [b-a for a,b in zip(l[:-1],l[1:])]
dlats = deltas(lats)
dlons = deltas(lons)
ddeps = deltas(deps)

if len(set(dlats)) > 1 or len(set(dlons)) > 1 or len(set(ddeps)) > 1:
    print("WARNING: model not uniformly gridded; output may look deformed")

dz = np.mean(ddeps)
dy = np.mean(dlats)
dx = np.mean(dlons)

vgrids_header = []
gap = "     "
vgrids_header.append(gap.join(["1","1"]))
vgrids_header.append(gap.join(list(map(str,[n_z,n_y,n_x]))))
vgrids_header.append(gap.join(list(map(str,[dz,np.radians(dy),np.radians(dx)]))))
vgrids_header.append(gap.join(list(map(str,[z0,np.radians(y0),np.radians(x0)]))))

vel_str_io = StringIO()
df = df.sort_values(by=["dep","lat","lon"],ascending=[False,True,True])
df["dv"].to_csv(vel_str_io,index=None,header=None)
vel_str_io.seek(0)
vel_str = vel_str_io.read()
vel_str_io.close()

velref_str = "0\n"*len(df)

if not os.path.exists(out_folder):
    os.mkdir(out_folder)
with open(os.path.join(out_folder,"vgridsref.in"),"w") as outfile:
    outfile.write("\n".join(vgrids_header + [velref_str]))
with open(os.path.join(out_folder,"vgrids.in"),"w") as outfile:
    outfile.write("\n".join(vgrids_header + [vel_str]))
