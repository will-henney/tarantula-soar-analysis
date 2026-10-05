"""
2026-10-03 - WJH

- First attempt to run scousepy line component fitting over the SOAR Ha cube
- Based on my experience with the tutorial

"""

from scousepy import scouse
from pathlib import Path
import warnings
import pyspeckit

from 

# Folder that contains this script
script_dir = Path(__file__).resolve().parent
# Assume we are a top-level subfolder of the project
project_dir = script_dir.parent

# Location of data cube: assume that it is in sibling folder of this project
cube_dir = project_dir.parent / "tarantula-soar-cubo"

# This is the continuum-subtracted cube in LSRK velocities
filename = "soar-30dor-ha-cube-lsrk-csub-bin016"
# scousepy wants folders to be strings with trailing slash
datadir = f"{cube_dir}/"
outputdir = f"{script_dir}/"

config_file = scouse.run_setup(filename, datadir, outputdir=outputdir)
s = scouse.stage_1(config=config_file, interactive=True)
s = scouse.stage_2(config=config_file)

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message="Mean of empty slice",
        category=RuntimeWarning,
    )
    warnings.filterwarnings(
        "ignore",
        message="invalid value encountered in divide",
        category=RuntimeWarning,
    )

    s = scouse.stage_3(config=config_file)

# # s = scouse.stage_4(config=config_file, bitesize=True)
# s = scouse.stage_4(config=config_file, bitesize=False)
