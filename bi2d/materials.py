from pathlib import Path

from .material import Material, ArrayInterpolate


beam = Material(0, name="beam")
vacuum = Material(0, name="vacuum")
steel = Material(0, name="steel", sigma=1.4e6)
copper = Material(0, name="copper", sigma=5.8e7)
titanium = Material(0, name="titanium", sigma=2.38e6)
_mu8c11_file = Path(__file__).parent / "material_files" / "mu8c11.csv"
_mu8c11_re_interpolate = ArrayInterpolate.from_file(_mu8c11_file, 0, 1, delimiter=r",")
_mu8c11_im_interpolate = ArrayInterpolate.from_file(_mu8c11_file, 0, 2, delimiter=r",")
mu8c11 = Material(0, name="mu8c11", eps_r=10,
                  mu_r_re=_mu8c11_re_interpolate.interp,
                  mu_r_im=_mu8c11_im_interpolate.interp)
_mischung43_file = Path(__file__).parent / "material_files" / "mischung43.csv"
_mischung43_re_interpolate = ArrayInterpolate.from_file(_mischung43_file, 0, 1, delimiter=r",")
_mischung43_im_interpolate = ArrayInterpolate.from_file(_mischung43_file, 0, 2, delimiter=r",")
mischung43 = Material(0, name="mischung43", eps_r=10,
                      mu_r_re=_mischung43_re_interpolate.interp,
                      mu_r_im=_mischung43_im_interpolate.interp)

__all__ = [beam, vacuum, steel, copper, titanium, mu8c11, mischung43]
