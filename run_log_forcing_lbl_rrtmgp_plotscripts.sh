#!/bin/bash
pythondir='/Users/jdelamere/Projects-AER/p1770/python-rrtmgp'
ptemplate='forcing_rrtmgp_lblrtm_config.ini.template'
pdir='/Users/jdelamere/Projects-AER/p1770/plot_share/validation_plots/'
lbldir='/Users/jdelamere/Projects-AER/p1770/runs_lblrtm_validations'
atmvernum='r1176'
rrtmgpvernum='r1403'
lblvernum='v12_6'
nangs='3'
pl_bands='--band 4 5'
y_axis='--log_y'
y_descrip='_log'

run_base=( \
preind-key \
preind-key )

run_forcing=( \
preind-key-minor-dbl_co2 \
preind-key-minor-quad_co2 )

test_forcing_model="RRTMGP_${rrtmgpvernum}"
reference_forcing_model="LBLRTM_${lblvernum}"

for ((i=0;i<${#run_base[@]};++i)); do
  echo "RUNNING...${run_base[$i]} and ${run_forcing[$i]} "
  inifile="forcing_rrtmgp_lblrtm_config_${run_forcing[$i]}-${run_base[$i]}.ini"
  if [ -f "${inifile}" ]; then
      rm "${inifile}"
  fi
  cp forcing_rrtmgp_lblrtm_config.ini.template "${inifile}"
  sed -i.bck "s#test_forcing_model = XYZ#test_forcing_model = ${test_model}#" ${inifile}
  sed -i.bck "s#reference_forcing_model = XYZ#reference_forcing_model = ${reference_model}#" ${inifile}
  sed -i.bck "s#atmosphere = XYZ#atmosphere = ${run_forcing[$i]}-${run_base[$i]}#" ${inifile}
  sed -i.bck "s#test_path = XYZ#test_path = ${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${run_base[$i]}/numangs_${nangs}/concatenate_fluxes/rrtmgp-lw-flux-inputs-outputs-${run_base[$i]}-all.nc#" ${inifile}
  sed -i.bck "s#reference_path = XYZ#reference_path = ${lbldir}/output_lbl_${lblvernum}_atm_${atmvernum}/${run_base[$i]}/numangs_${nangs}/concatenate_fluxes/lblrtm-lw-flux-inputs-outputs-${run_base[$i]}-all.nc#" ${inifile}
  sed -i.bck "s#test_force_path = XYZ#test_force_path = ${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${run_forcing[$i]}/numangs_${nangs}/concatenate_fluxes/rrtmgp-lw-flux-inputs-outputs-${run_forcing[$i]}-all.nc#" ${inifile}
  sed -i.bck "s#reference_force_path = XYZ#reference_force_path = ${lbldir}/output_lbl_${lblvernum}_atm_${atmvernum}/${run_forcing[$i]}/numangs_${nangs}/concatenate_fluxes/lblrtm-lw-flux-inputs-outputs-${run_forcing[$i]}-all.nc#" ${inifile}
  sed -i.bck "s#profiles_prefix = XYZ#profiles_prefix = forcing_profs${y_descrip}_lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${run_forcing[$i]}-${run_base[$i]}#" ${inifile}
  sed -i.bck "s#stats_prefix = XYZ#stats_prefix = forcing_stats_lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${run_forcing[$i]}-${run_base[$i]}#" ${inifile}
  sed -i.bck "s#stats_csv = XYZ#stats_csv = forcing_lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${run_forcing[$i]}-${run_base[$i]}_diff.csv#" ${inifile}
  rm "${inifile}.bck"
  mv "${inifile}" "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${run_forcing[$i]}/numangs_${nangs}/."
  cd "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${run_forcing[$i]}/numangs_${nangs}"
  "${pythondir}/LBLRTM_RRTMGP_compare.py" --config_file "${inifile}" --single_stat --plot_stats --plot_profiles --ref_forcing ${pl_bands} ${y_axis}
  cd "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}"
done
