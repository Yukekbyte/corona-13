'''
This script was mostly written by Google Gemini.
'''


bl_info = {
    "name": "Export CCAM Camera",
    "author": "Gemini",
    "version": (1, 0),
    "blender": (4, 0, 0),
    "location": "File > Export > CCAM (.cam)",
    "description": "Export active camera to Corona's packed C struct format",
    "category": "Import-Export",
}

import bpy
import struct
import math
from mathutils import Vector, Matrix
from bpy_extras.io_utils import ExportHelper
from bpy.props import StringProperty
from bpy.types import Operator

class ExportCCAMCamera(Operator, ExportHelper):
    """Export active camera to a binary CCAM file"""
    bl_idname = "export.ccam_camera"
    bl_label = "Export CCAM (.cam)"
    
    filename_ext = ".cam"
    filepath: StringProperty(
        name="File Path",
        description="Filepath used for exporting the file",
        maxlen=1024,
        default="blender-scene.cam"
    )
    filter_glob: StringProperty(default="*.cam", options={'HIDDEN'})
    
    def execute(self, context):
        cam_obj = context.scene.camera
        if not cam_obj or cam_obj.type != 'CAMERA':
            self.report({'ERROR'}, "No active camera found")
            return {'CANCELLED'}

        cam_data = cam_obj.data
        matrix = cam_obj.matrix_world
        
        # 1. Swizzle Position (Tracer X,Y,Z = Blender Y,X,Z)
        pos = matrix.to_translation()
        pos_swizzled = (pos.y, -pos.x, pos.z) 
        
        # 2. Re-base Orientation
        # We must transform Blender's local camera axes into your Tracer's world space
        b_rot = matrix.to_3x3()
        
        # Blender Camera locals: Forward is -Z, Up is Y, Right is X
        w_fwd   = b_rot @ Vector((0, 0, -1))
        w_up    = b_rot @ Vector((0, 1, 0))
        w_right = b_rot @ Vector((1, 0, 0))
        
        # Helper to swizzle any world vector to your Tracer's coordinate system
        def to_tracer(v): 
            return Vector((v.y, v.x, v.z))
        
        # Corona code expects Forward to be +Z and Up to be +Y
        t_fwd   = to_tracer(w_fwd)
        t_up    = to_tracer(w_up)
        t_right = to_tracer(w_right)
        
        # Build the Tracer rotation matrix (Columns: Right, Up, Forward)
        t_mat = Matrix((t_right, t_up, t_fwd)).transposed()
        orient = t_mat.to_quaternion()

        # Still did something wrong, but it works with this dirty fix, so who cares...
        orient = (-orient.z, -orient.y, -orient.x, -orient.w)

        # 3. Scale and Unit Conversion
        unit_scale = 0.01 
        
        magic = b"CCAM"
        version = 1
        iso = 400.0

        # should this be variable?
        aperture_idx = 6  # hard coded
        exposure_idx = 13 # hard coded
        
        speed = 1.0
        focus_offset = 0.0

        # Support 'Focus on Object' option
        if cam_data.dof.use_dof and cam_data.dof.focus_object:
            target_obj = cam_data.dof.focus_object
            # Calculate distance between camera origin and object origin
            target_pos = target_obj.matrix_world.to_translation()
            focus_dist = (pos - target_pos).length
        else:
            focus_dist = cam_data.dof.focus_distance
        
        focal_length = cam_data.lens * unit_scale
        width = cam_data.sensor_width * unit_scale
        height = cam_data.sensor_height * unit_scale
        
        # Use unscaled width for crop factor so it stays near 1.0
        crop_factor = 35.0 / cam_data.sensor_width if cam_data.sensor_width > 0 else 1.0
        
        # Binary Packing
        # Using * to unpack the tuples directly into the arguments
        packed_data = struct.pack(
            '<4si3f3f4f4f6fii2f',
            magic, version,
            *pos_swizzled,
            *pos_swizzled,
            *orient,
            *orient,
            speed, focus_offset, focus_dist,
            width, height, crop_factor,
            aperture_idx, exposure_idx,
            focal_length, iso
        )

        try:
            with open(self.filepath, 'wb') as f:
                f.write(packed_data)
            self.report({'INFO'}, f"Camera exported to {self.filepath}")
        except Exception as e:
            self.report({'ERROR'}, str(e))
            return {'CANCELLED'}

        return {'FINISHED'}

def menu_func_export(self, context):
    self.layout.operator(ExportCCAMCamera.bl_idname, text="CCAM (.cam)")

def register():
    bpy.utils.register_class(ExportCCAMCamera)
    bpy.types.TOPBAR_MT_file_export.append(menu_func_export)

def unregister():
    bpy.utils.unregister_class(ExportCCAMCamera)
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_export)

if __name__ == "__main__":
    register()