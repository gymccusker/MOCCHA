"""
List of variables from the UKV
==============================

"""

Stash_list = {'m01s16i222': 'air_pressure_at_sea_level',
              'm01s00i024': 'surface_temperature',
              'm01s03i236': 'air_temperature',
              'm01s03i245': 'relative_humidity',
              'm01s03i225': 'eastward_wind',
              'm01s03i226': 'northward_wind',
              'm01s00i025': 'atmosphere_boundary_layer_thickness',
              'm01s03i476': 'atmosphere_boundary_layer_type',
              'm01s03i217': 'surface_upward_sensible_heat_flux',
              'm01s03i234': 'surface_upward_latent_heat_flux',
              'm01s03i202': 'downward_heat_flux_in_soil',
              'm01s09i217': 'cloud_area_fraction',
              'm01s09i203': 'cloud_area_fraction_in_atmosphere_layer',
              'm01s09i204': 'cloud_area_fraction_in_atmosphere_layer',
              'm01s09i205': 'cloud_area_fraction_in_atmosphere_layer',
              'm01s03i026': 'surface_roughness_length',
              'm01s08i223': 'lwe_thickness_of_moisture_content_of_soil_layer',
              'm01s08i225': 'soil_temperature',
              'm01s15i201': 'eastward_wind',
              'm01s15i202': 'northward_wind',
              'm01s03i316': 'surface_temperature_on_tiles',
              'm01s03i290': 'surface_sensible_heat_flux_on_tiles',
              'm01s03i321': 'canopy_water_on_tiles',
              'm01s04i201': 'large_scale_rainfall_amount',
              'm01s04i202': 'large_scale_snowfall_amount',
              'm01s03i241': 'water_evaporation_amount',
              'm01s16i203': 'air_temperature',
              'm01s16i204': 'relative_humidity',
              'm01s00i409': 'surface_air_pressure',
              'm01s16i004': 'air_temperature',      #included
              'm01s00i407': 'air_pressure',
              'm01s00i408': 'air_pressure',
              'm01s00i090': 'total aerosol (for visibility)',
              'm01s00i010': 'specific_humidity',
              'm01s00i002': 'eastward_wind',
              'm01s00i003': 'northward_wind',
              'm01s00i266': 'large_scale_cloud_area_fraction',
              'm01s00i254': 'mass_fraction_of_cloud_liquid_water_in_air',
              'm01s00i272': 'rain',
              'm01s00i012': 'mass_fraction_of_cloud_ice_in_air',
              'm01s09i229': 'relative_humidity',
              'm01s00i150': 'upward_air_velocity',
              'm01s03i216': 'boundary_layer_heat_fluxes',
              'm01s03i222': 'boundary_layer_total_moisture_fluxes',
              'm01s01i235': 'downwelling_SW_radiation',
              'm01s01i201': 'net_SW_radiation',
              'm01s02i207': 'downwelling_LW_radiation',
              'm01s02i201': 'net_LW_radiation',
              'm01s03i248': 'screen_fog_fraction',
              'm01s03i247': 'screen_visibility',
              'm01s00i022': 'canopy_water',
              'm01s03i328': 'screen_temp_on_tiles',
              'm01s00i234': 'roughness_length_on_tiles',
              'm01s03i314': 'surface_net_radiation_on_tiles',
              'm01s03i330': 'surface_latent_heat_on_tiles'
              }


def returnWantedStash():

    return Stash_list.keys()


def findfieldName(stash):

    if stash in Stash_list.keys():
        return Stash_list[stash]
    else:
        return None
