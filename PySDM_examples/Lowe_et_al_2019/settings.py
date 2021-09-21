import numpy as np
from pystrict import strict
from PySDM.initialisation import spectral_sampling as spec_sampling
from PySDM.physics import si, Formulae, constants as const
from PySDM_examples.Lowe_et_al_2019.aerosol import _Aerosol
from PySDM.dynamics.condensation import default_rtol_x, default_rtol_thd


@strict
class Settings:
    def __init__(self, dt: float, n_sd_per_mode: int,
                 aerosol: _Aerosol,
                 model: str,
                 spectral_sampling: type(spec_sampling.SpectralSampling),
                 w: float = 0.32 * si.m / si.s,
                 rtol_x: float = default_rtol_x,
                 rtol_thd: float = default_rtol_thd
                 ):
        assert model in ('bulk', 'film')
        self.model = model
        self.n_sd_per_mode = n_sd_per_mode
        self.formulae = Formulae(surface_tension='CompressedFilm' if model=='film' else 'Constant')
        self.aerosol = aerosol
        self.spectral_sampling = spectral_sampling
        self.t_max = int(210 / w) * si.m
        self.output_interval = 10 * si.s
        self.dt = dt
        self.rtol_x = rtol_x
        self.rtol_thd = rtol_thd

        self.w = w
        self.g = 9.81 * si.m / si.s**2

        self.p0 = 980 * si.mbar
        self.T0 = 280 * si.K
        pv0 = .99 * self.formulae.saturation_vapour_pressure.pvs_Celsius(self.T0 - const.T0)
        self.q0 = const.eps * pv0 / (self.p0 - pv0)

        self.cloud_radius_range = (
                .5 * si.micrometre,
                np.inf
        )

        self.mass_of_dry_air = 44
        
        self.wet_radius_bins_edges = np.logspace(np.log10(4 * si.um), np.log10(12 * si.um), 128+1, endpoint=True)

    @property
    def rho0(self):
        rhod0 = self.formulae.trivia.p_d(self.p0, self.q0) / self.T0 / const.Rd
        return rhod0 * (1 + self.q0)

    @property
    def nt(self):
        nt = self.t_max / self.dt
        assert nt == int(nt)
        return int(nt)

    @property
    def steps_per_output_interval(self) -> int:
        return int(self.output_interval / self.dt)

    @property
    def output_steps(self) -> np.ndarray:
        return np.arange(0, self.nt + 1, self.steps_per_output_interval)

