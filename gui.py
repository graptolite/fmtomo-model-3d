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

from tkinter import *
import os
import shutil
from materials import materials
from processing_functions import *

if not os.path.exists("tmp"):
    os.mkdir("tmp")

class GUI(Tk):
    def __init__(self):
        super().__init__()
        w,h = 400, 600
        self.geometry("%ux%u" % (w,h))
        self.title("FMTOMO 3D")
        self.inputs_frame = Frame(self,bg="lightgrey",borderwidth=5)
        self.inputs_frame.columnconfigure(0,weight=3)
        inputs_w = 0.6
        self.inputs_frame.place(relheight=1,relwidth=inputs_w,relx=0,y=0)
        self.l_input = Label(self.inputs_frame,text="Input Filepath:",font=("bold"))
        self.inp = Entry(self.inputs_frame)
        self.btn = Button(self.inputs_frame,text="Load files",width=10,height=1,command=self.load_files)
        self.l_isosurf = Label(self.inputs_frame,text="dv isosurface (km/s)")
        self.isosurface_val = DoubleVar(value=0.2)
        self.isosurface = Entry(self.inputs_frame,textvariable=self.isosurface_val)
        self.l_mat = Label(self.inputs_frame,text="material")
        self.material_val = StringVar(value=list(materials.keys())[0])
        self.material = OptionMenu(self.inputs_frame,self.material_val,*materials.keys())
        self.l_downscale = Label(self.inputs_frame,text="Blender downscale (km to 1 m):")
        self.blender_downscale_val = DoubleVar(value=100)
        self.blender_downscale = Entry(self.inputs_frame,textvariable=self.blender_downscale_val)
        self.do_render = IntVar(value=0)
        self.render = Checkbutton(self.inputs_frame,text="Render",variable=self.do_render,onvalue=1,offvalue=0)
        self.open_gui = IntVar(value=1)
        self.blender_gui = Checkbutton(self.inputs_frame,text="Blender GUI",variable=self.open_gui,onvalue=1,offvalue=0)
        self.btn_3d = Button(self.inputs_frame,text="Load 3D",width=10,height=1,command=self.load_3d)
        input_widget_list = [self.l_isosurf,self.isosurface,
                             self.l_input,self.inp,self.btn,
                             self.l_mat,self.material,
                             self.l_downscale,self.blender_downscale,
                             self.render,self.blender_gui,
                             self.btn_3d]
        self.stack_widgets(input_widget_list)
        self.outputs_frame = Frame(self,bg="lightblue",borderwidth=5)
        self.outputs_frame.columnconfigure(0,weight=3)
        self.outputs_frame.place(relheight=1,relwidth=1-inputs_w,relx=inputs_w,y=0)
        self.l_msg = Label(self.outputs_frame,text="Messages",bg="lightblue")
        self.msg = Canvas(self.outputs_frame,bg="white")
        output_widget_list = [self.l_msg,self.msg]
        self.stack_widgets(output_widget_list)
        self.update_idletasks()
        wm,hm = self.msg.winfo_width(),self.msg.winfo_height()
        self.text_placeholder = self.msg.create_text(wm/2,hm/2,anchor=CENTER,width=(1-inputs_w)*w)
        self.protocol("WM_DELETE_WINDOW",self.destroy)
        self.fp = ""
        return
    def load_files(self):
        ''' Copy the required grid files from the fmtomo working directory to the temp dir.

        Returns: <None>
        '''
        # Retrieve the path to the fmtomo working dir.
        fp = self.inp.get()
        self.fp = fp
        # Validate that the fmtomo working dir exists.
        if fp and os.path.exists(fp):
            try:
                # Clear any old files that may be present in the temp dir.
                for f in os.listdir("tmp"):
                    os.remove(os.path.join("tmp",f))
                # Attempt to load the required grid files.
                for f in ["vgrids.in","vgridsref.in","propgrid.in"]:
                    shutil.copy(os.path.join(fp,f),"tmp")
                # Load the optional grid file if present.
                v_true = os.path.join(fp,"vgridstrue.in")
                if os.path.exists(v_true):
                    shutil.copy(v_true,"tmp")
                self.update_msg("Folder successfully loaded")
            except:
                self.update_msg("Folder does not contain one or more of vgrids.in, vgridsref.in or propgrid.in")
        else:
            self.update_msg("Folder does not exist.")
        return
    def load_3d(self):
        ''' Handle the 3D modelling from the grid files present in the temp dir.

        Returns: <None>
        '''
        self.update_msg("Loading 3D...")
        # Get all of the input and control parameters.
        isosurface = float(self.isosurface.get())
        material = materials[self.material_val.get()]
        blender_downscale = self.blender_downscale_val.get()
        render = bool(self.do_render.get())
        open_gui = bool(self.open_gui.get())
        isosurface_specs = {isosurface:material}
        print(isosurface,material,blender_downscale,render,open_gui)
        # Try running the 3D modelling.
        try:
            exec_3d(blender_downscale,isosurface_specs,render,open_gui)
            self.update_msg("Finished blender rendering")
        except FileNotFoundError:
            self.update_msg("No FMTOMO folder loaded yet, or input params are missing")
        return
    def stack_widgets(self,widget_list):
        ''' Stack tkinter widgets above each other.

        widget_list | <list> | list of tkinter widgets.

        Returns: <None>
        '''
        for i,w in enumerate(widget_list):
            w.grid(column=0,row=0+i)
        return
    def update_msg(self,text):
        ''' Replace the content of the message in the GUI with new text.

        text | <str> | text to replace the previous message in the GUI with.

        Returns: <None>
        '''
        self.msg.itemconfig(self.text_placeholder,text=text)


if __name__=="__main__":
    root = GUI()
    root.mainloop()
