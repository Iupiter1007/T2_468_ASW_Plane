# make_737_vsp.py
import openvsp as vsp
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def safe_set_parm(geom_id, parm_name, group_name, value, verbose=True):
    """Try SetParmVal and print friendly error + available parms if not found."""
    try:
        vsp.SetParmVal(geom_id, parm_name, group_name, value)
        if verbose:
            print(f"Set {geom_id} / {group_name} / {parm_name} = {value}")
    except Exception as e:
        print(f"WARNING: Could not set parm '{parm_name}' in group '{group_name}' for geom '{geom_id}': {e}")
        # Try to show available parms for that geom to help you find the right name
        try:
            print("Available parm names (group / name) for this geom (best-effort):")
            # Many builds offer FindParm/GetParm helpers; we'll iterate common groups
            groups = ["Design", "XForm", "WingGeom", "XSec_1", "XSec_2", "XSec_3", "Airfoil", "Geometry"]
            shown = False
            for g in groups:
                try:
                    pname = vsp.FindParm(geom_id, parm_name, g)
                    # If FindParm found something then the parm exists; skip (we failed earlier)
                except Exception:
                    pass
            # If more detailed listing functions exist you can add them here; otherwise ask the API for doc
            print("  - run the OpenVSP GUI, open this geom and expand the Parameters tree to find exact names")
        except Exception:
            pass

# Start
print("OpenVSP version:", vsp.GetVSPVersion())
vsp.ClearVSPModel()

# --- Fuselage ---
fuse = vsp.AddGeom("FUSELAGE")   # returns geom id string
# 737-800-ish: length ~39.5 m, max radius ~2.6 m
safe_set_parm(fuse, "Length", "Design", 39.5)
safe_set_parm(fuse, "Max_Radius", "Design", 2.6)   # guess — if missing adjust below
# If the exact "Max_Radius" parm isn't in your build, edit the fuselage XSec widths in GUI or inspect parms.

# --- Main Wing ---
wing = vsp.AddGeom("WING")
# 737-800-ish geometry
wing_span = 35.8
wing_root_chord = 6.5
wing_tip_chord = 2.6
wing_sweep_le = 25.0   # degrees (leading edge)
wing_dihedral = 3.25   # degrees

# common parm names used across many OpenVSP versions; if one fails, safe_set_parm will warn
safe_set_parm(wing, "TotalSpan", "WingGeom", wing_span)
# chord parms are often per XSec; set root/tip using common XSec parms (if present)
safe_set_parm(wing, "Root_Chord", "XSec_1", wing_root_chord)
safe_set_parm(wing, "Tip_Chord", "XSec_2", wing_tip_chord)
safe_set_parm(wing, "Sweep_LE", "WingGeom", wing_sweep_le)
safe_set_parm(wing, "Dihedral", "WingGeom", wing_dihedral)

# place wing in fore/aft on fuselage: set XRelLocation (meters aft of nose or relative system)
# common XForm parms include X_Rel_Location in group "XForm"
safe_set_parm(wing, "X_Rel_Location", "XForm", 12.0)  # attach approx ~12 m aft of nose
safe_set_parm(wing, "Z_Rel_Location", "XForm", 0.0)   # set vertical placement (mid fuselage)

# -- Horizontal tail (modeled as a small wing) --
htail = vsp.AddGeom("WING")
h_span = 13.3
h_root = 3.5
h_tip = 1.8
h_sweep_le = 15.0
safe_set_parm(htail, "TotalSpan", "WingGeom", h_span)
safe_set_parm(htail, "Root_Chord", "XSec_1", h_root)
safe_set_parm(htail, "Tip_Chord", "XSec_2", h_tip)
safe_set_parm(htail, "Sweep_LE", "WingGeom", h_sweep_le)
# move to tail: X near aft fuselage, small z offset
safe_set_parm(htail, "X_Rel_Location", "XForm", 31.8)  # ~ 31.8 m aft of nose
safe_set_parm(htail, "Z_Rel_Location", "XForm", 1.0)   # slightly above fuselage centerline
safe_set_parm(htail, "Dihedral", "WingGeom", 2.0)

# -- Vertical tail (use WING and rotate so it stands up) --
vtail = vsp.AddGeom("WING")
v_height = 7.7
v_root = 4.2
v_tip = 1.6
# We'll create as a wing then rotate 90deg about Y to stand it up
safe_set_parm(vtail, "TotalSpan", "WingGeom", v_height)   # use TotalSpan to control height when rotated
safe_set_parm(vtail, "Root_Chord", "XSec_1", v_root)
safe_set_parm(vtail, "Tip_Chord", "XSec_2", v_tip)
safe_set_parm(vtail, "Sweep_LE", "WingGeom", 10.0)
# position it near tail root
safe_set_parm(vtail, "X_Rel_Location", "XForm", 31.0)
safe_set_parm(vtail, "Z_Rel_Location", "XForm", 0.0)
# Rotate the vertical tail: many API installs let you set X_Rel_Rotation / Y_Rel_Rotation / Z_Rel_Rotation
safe_set_parm(vtail, "Y_Rel_Rotation", "XForm", 90.0)  # rotate so span is in Z instead of Y

# --- Engines as PODs (two) ---
pod_l = vsp.AddGeom("POD")
pod_r = vsp.AddGeom("POD")
# common pod params: Length in "Design" and Radius or Diameter in "Design" or "XSec_1"
safe_set_parm(pod_l, "Length", "Design", 3.4)
safe_set_parm(pod_r, "Length", "Design", 3.4)
safe_set_parm(pod_l, "Radius", "Design", 1.1)
safe_set_parm(pod_r, "Radius", "Design", 1.1)

# Move engine pods under the wing (set X/Y/Z rel location); quarter-chord / nacelle placement
# place around wing X location (12 m from nose), and at about 1/4 span
safe_set_parm(pod_l, "X_Rel_Location", "XForm", 12.0)
safe_set_parm(pod_r, "X_Rel_Location", "XForm", 12.0)
safe_set_parm(pod_l, "Y_Rel_Location", "XForm", wing_span*0.25)
safe_set_parm(pod_r, "Y_Rel_Location", "XForm", -wing_span*0.25)
safe_set_parm(pod_l, "Z_Rel_Location", "XForm", -1.9)
safe_set_parm(pod_r, "Z_Rel_Location", "XForm", -1.9)

# Small rotation for the pods if desired (angle in deg)
safe_set_parm(pod_l, "Z_Rel_Rotation", "XForm", 0.0)
safe_set_parm(pod_r, "Z_Rel_Rotation", "XForm", 0.0)

# Recompute
vsp.Update()

# Save .vsp3 and export STL (or other) - ensure your OpenVSP build supports EXPORT_STL
out_vsp = "737_like.vsp3"
vsp.WriteVSPFile(out_vsp, vsp.SET_ALL)
print("Wrote", out_vsp)

out_stl = "737_like.stl"
try:
    vsp.ExportFile(out_stl, vsp.SET_ALL, vsp.EXPORT_STL)
    print("Exported STL:", out_stl)
except Exception as e:
    print("STL export failed (API may differ). You can open the .vsp3 in OpenVSP GUI and export manually. Error:", e)

print("Done. If any parm warnings occurred above, check the geometry in the GUI and inspect parameter names for differences.")




def plot_stl(filename):
    triangles = []

    with open(filename, "r") as f:
        for line in f:
            if line.strip().startswith("vertex"):
                _, x, y, z = line.split()
                triangles.append([float(x), float(y), float(z)])

    triangles = np.array(triangles).reshape(-1, 3, 3)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    mesh = Poly3DCollection(triangles, alpha=0.5)
    ax.add_collection3d(mesh)

    # Auto-scale axes
    scale = triangles.flatten()
    ax.auto_scale_xyz(scale, scale, scale)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    ax.set_box_aspect([1,1,1])
    ax.view_init(elev=45, azim=45)
    plt.tight_layout()
    plt.show()

plot_stl(out_stl)
