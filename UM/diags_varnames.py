"""
List of variables names to be used
==============================

"""

Var_list = { # paXXX
              'm01s02i207': 'surface_downwelling_LW_radiation',
              'm01s01i235': 'surface_downwelling_SW_radiation',
              'm01s02i201': 'surface_net_LW_radiation',
              'm01s01i201': 'surface_net_SW_radiation',
              'm01s01i207': 'toa_incoming_shortwave_flux',
              'm01s02i205': 'toa_outgoing_longwave_flux',
              'm01s01i208': 'toa_outgoing_shortwave_flux',
               # pbXXX
              'm01s03i241': 'water_evaporation_amount',
              'm01s03i304': 'turbulent_mixing_height_after_BL',
              'm01s03i360': 'h_decoupled_layer_base',
              'm01s03i361': 'h_sc_cloud_base',
              'm01s03i476': 'bl_type',
              'm01s09i216': 'cloud_area_fraction_assuming_random_overlap',
              'm01s09i217': 'cloud_area_fraction_assuming_maximum_random_overlap',
              'm01s09i221': 'wet_bulb_freezing_level_altitude',
              'm01s30i461': 'total_column_q',
              'm01s02i391': 'IWP',
              'm01s02i392': 'LWP',
              'm01s16i222': 'air_pressure_at_sea_level',
              'm01s03i236': 'temp_1.5m',
              'm01s03i025': 'bl_depth',
              'm01s03i250': 'dew_point_temperature_at_1.5m',
              'm01s03i247': 'vis_1.5m',
              'm01s03i248': 'fogfrac_1.5m',
              'm01s09i205': 'high_cloud',
              'm01s09i203': 'low_cloud',
              'm01s09i204': 'medium_cloud',
              'm01s03i245': 'rh_1.5m',
              'm01s03i237': 'specific_humidity_at_1.5m',
              'm01s04i203': 'rainfall_flux',
              'm01s04i204': 'snowfall_flux',
              'm01s00i409': 'sfc_pressure',
              'm01s00i024': 'sfc_temperature',
              'm01s03i234': 'latent_heat_flux',
              'm01s03i217': 'sensible_heat_flux',
              'm01s03i225': 'eastward_wind_at_10m',
              'm01s03i226': 'northward_wind_at_10m',
               # pcXXX  -- CLOUDNET
              'm01s04i118': 'radr_refl',                                        # total_radar_reflectivity
              'm01s00i266': 'cloud_fraction',                                   # large_scale_cloud_area_fraction
              'm01s00i408': 'pressure',                                         # air_pressure
              'm01s16i004': 'temperature',                                      # air_temperature
              'm01s00i012': 'qice',                                             # mass_fraction_of_cloud_ice_in_air
              'm01s00i254': 'qliq',                                             # mass_fraction_of_cloud_liquid_water_in_air
              'm01s00i010': 'q',                                                # specific_humidity
              'm01s00i150': 'wwind',                                            # upward_air_velocity
              'm01s00i002': 'uwind',                                            # eastward_wind
              'm01s00i003': 'vwind',                                            # northward_wind
              'm01s00i267': 'liquid_cloud_fraction',
              'm01s00i268': 'ice_cloud_fraction',
              # pdXXX -- BOUNDARY LAYER
              'm01s03i362': 'entrainment_rate_SML',
              'm01s03i363': 'entrainment_rate_BL',
              'm01s03i464': 'obukhov_length',
              'm01s03i219': 'atmosphere_downward_eastward_stress',              # ro 1-70
              'm01s03i220': 'atmosphere_downward_northward_stress',             # ro 1-70
              'm01s03i473': 'tke',                         # ro 1-70
              # 'm01s00i407': 'air_pressure',                                     # ro 1-70
              'm01s03i460': 'surface_downward_eastward_stress',
              'm01s03i461': 'surface_downward_northward_stress',
              'm01s03i223': 'surface_upward_water_flux',
              'm01s00i026': 'surface_roughness_length',
              'm01s03i208': 'bulk_richardson_number',
              'm01s03i501': 'mixing_length_for_momentum',
              'm01s03i465': 'explicit_friction_velocity',
              'm01s03i471': 'BL_momentum_diffusion',
              'm01s03i469': 'vertical_buoyancy_gradient',
              'm01s00i031': 'sea_ice_fraction',
              # peXXX -- CASIM MICROPHYSICS
              'm01s00i075': 'qnliq',
              'm01s00i078': 'qnice',
              'm01s00i083': 'qnsolaeroliq',
              # 'm01s00i086': 'activated_soluble_aerosol_in_ice',
              'm01s00i084': 'qnsolaerorain',
              'm01s00i088': 'qnsolaero',
              'm01s00i271': 'qicecrystals'
              }


def returnWantedVarnames():

    return Var_list.keys()


def findfieldName(stash):

    if stash in Var_list.keys():
        return Var_list[stash]
    else:
        return None
