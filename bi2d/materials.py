from pathlib import Path

from .material import Material, ArrayInterpolate


beam = Material(0, name="beam")
vacuum = Material(0, name="vacuum")
steel_isotropic = Material(0, name="steel_isotropic", sigma=1.4e6)
copper_isotropic = Material(0, name="copper_isotropic", sigma=5.8e7)
titanium_isotropic = Material(0, name="titanium_isotropic", sigma=2.38e6)
_mu8c11_file = Path(__file__).parent / "material_files" / "mu8c11.dat"
_mu8c11_re_interpolate = ArrayInterpolate.from_file(_mu8c11_file, 0, 1)
_mu8c11_im_interpolate = ArrayInterpolate.from_file(_mu8c11_file, 0, 2)
mu8c11 = Material(0, name="mu8c11",
                  mu_r_re=_mu8c11_re_interpolate.interp,
                  mu_im=_mu8c11_im_interpolate.interp)

__all__ = [beam, vacuum, steel_isotropic, copper_isotropic, titanium_isotropic, mu8c11]
