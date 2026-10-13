# -*- coding: utf-8 -*-
"""Stage the Shiny calculator: copy the released continuous-curve table into
app/data/. Run from anywhere; paths are relative to the repository root."""
import os
import shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
appdir = os.path.join(ROOT, "app")
datadir = os.path.join(appdir, "data")
os.makedirs(datadir, exist_ok=True)
src = os.path.join(ROOT, "outputs", "M1", "M1c_continuous_curves.csv")
shutil.copy(src, os.path.join(datadir, "continuous_curves.csv"))
print("app dir:", appdir)
print("files:", os.listdir(datadir))
