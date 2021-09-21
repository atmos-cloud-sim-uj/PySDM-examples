from PySDM.environments import Parcel
from PySDM import Builder
from PySDM.backends import CPU
from PySDM.dynamics import AmbientThermodynamics, Condensation
from PySDM.initialisation.multiplicities import discretise_n
from PySDM.initialisation.r_wet_init import r_wet_init
from PySDM.physics.spectra import Sum
import PySDM.products as PySDM_products
import numpy as np


class Simulation:
    def __init__(self, settings, products=None):
        env = Parcel(dt=settings.dt, mass_of_dry_air=settings.mass_of_dry_air, p0=settings.p0, q0=settings.q0,
                     T0=settings.T0, w=settings.w, g=settings.g)
        n_sd = settings.n_sd_per_mode * len(settings.aerosol.aerosol_modes_per_cc)
        builder = Builder(n_sd=n_sd, backend=CPU, formulae=settings.formulae)
        builder.set_environment(env)

        attributes = {
            'dry volume':np.empty(0),
            'dry volume organic':np.empty(0),
            'kappa times dry volume':np.empty(0),
            'n': np.ndarray(0)
        }
        for mode in settings.aerosol.aerosol_modes_per_cc:
            r_dry, n_in_dv = settings.spectral_sampling(spectrum=mode['spectrum']).sample(settings.n_sd_per_mode)
            n_in_dv /= (settings.rho0 / settings.mass_of_dry_air)
            v_dry = settings.formulae.trivia.volume(radius=r_dry)
            attributes['n'] = np.append(attributes['n'], n_in_dv)
            attributes['dry volume'] = np.append(attributes['dry volume'], v_dry)
            attributes['dry volume organic'] = np.append(attributes['dry volume organic'], mode['f_org'] * v_dry)
            attributes['kappa times dry volume'] = np.append(attributes['kappa times dry volume'], v_dry * mode['kappa'][settings.model])
        for attribute in attributes.values():
            assert attribute.shape[0] == n_sd

        attributes['n'] = discretise_n(attributes['n'])

        dv = settings.mass_of_dry_air / settings.rho0
        np.testing.assert_approx_equal(
            np.sum(attributes['n']) / dv,
            Sum(
                tuple([settings.aerosol.aerosol_modes_per_cc[i]['spectrum'] for i in range(len(settings.aerosol.aerosol_modes_per_cc))])
            ).norm_factor,
            significant=5
        )
        r_wet = r_wet_init(
            r_dry=settings.formulae.trivia.radius(volume=attributes['dry volume']),
            environment=env,
            kappa_times_dry_volume=attributes['kappa times dry volume'],
            f_org=attributes['dry volume organic'] / attributes['dry volume']
        )
        attributes['volume'] = settings.formulae.trivia.volume(radius=r_wet)

        if settings.model == 'bulk':
            del attributes['dry volume organic']

        builder.add_dynamic(AmbientThermodynamics())
        builder.add_dynamic(Condensation(rtol_x=settings.rtol_x, rtol_thd=settings.rtol_thd))

        products = products or (
            PySDM_products.ParcelDisplacement(),
            PySDM_products.Time(),
            PySDM_products.PeakSupersaturation(),
            PySDM_products.CloudDropletConcentration(radius_range=settings.cloud_radius_range),
            PySDM_products.ParticlesWetSizeSpectrum(radius_bins_edges=settings.wet_radius_bins_edges),
        )

        self.particulator = builder.build(attributes=attributes, products=products)
        self.settings = settings

    def _save_scalars(self, output):
        for k, v in self.particulator.products.items():
            if len(v.shape) > 1:
                continue
            value = v.get()
            if isinstance(value, np.ndarray) and value.size == 1:
                value = value[0]
            output[k].append(value)

    def _save_spectrum(self, output):
        value = self.particulator.products['Particles Wet Size Spectrum'].get()
        output['spectrum'] = value

    def run(self):
        output = {k: [] for k in self.particulator.products.keys()}
        for step in self.settings.output_steps:
            self.particulator.run(step - self.particulator.n_steps)
            self._save_scalars(output)
        self._save_spectrum(output)
        return output
