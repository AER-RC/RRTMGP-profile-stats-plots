#!/bin/bash
pythondir='/Users/jdelamere/Projects-AER/p1770/python-rrtmgp'
ptemplate='rrtmgp_lblrtm_config.ini.template'
pdir='/Users/jdelamere/Projects-AER/p1770/plot_share/validation_plots/'
lbldir='/Users/jdelamere/Projects-AER/p1770/runs_lblrtm_validations'
atmvernum='r1176'
rrtmgpvernum='r1392'
lblvernum='v12_6'
nangs='3'


runs=( \
garand-key \
preind-key \
preind-key-minor-dbl_co2 \
preind-key-minor-quad_co2 )

test_model="RRTMGP_${rrtmgpvernum}"
reference_model="LBLRTM_${lblvernum}"

for srun in "${runs[@]}"; do
  echo "RUNNING...${srun}"
  inifile="rrtmgp_lblrtm_config_${srun}.ini"
  if [ -f "${inifile}" ]; then
      rm "${inifile}"
  fi
  cp rrtmgp_lblrtm_config.ini.template "${inifile}"
  sed -i.bck "s#test_model = XYZ#test_model = ${test_model}#" ${inifile}
  sed -i.bck "s#reference_model = XYZ#reference_model = ${reference_model}#" ${inifile}
  sed -i.bck "s#atmosphere = XYZ#atmosphere = ${srun}_${atmvernum}#" ${inifile}
  sed -i.bck "s#test_path = XYZ#test_path = ${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${srun}/numangs_${nangs}/concatenate_fluxes/rrtmgp-lw-flux-inputs-outputs-${srun}-all.nc#" ${inifile}
  sed -i.bck "s#reference_path = XYZ#reference_path = ${lbldir}/output_lbl_${lblvernum}_atm_${atmvernum}/${srun}/numangs_${nangs}/concatenate_fluxes/lblrtm-lw-flux-inputs-outputs-${srun}-all.nc#" ${inifile}
  sed -i.bck "s#profiles_prefix = XYZ#profiles_prefix = profs_lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${srun}#" ${inifile}
  sed -i.bck "s#stats_prefix = XYZ#stats_prefix = stats_lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${srun}#" ${inifile}
  sed -i.bck "s#stats_csv = XYZ#stats_csv = lbl_${lblvernum}_rrtmgp_${rrtmgpvernum}_${srun}_diff.csv#" ${inifile}
  rm "${inifile}.bck"
  mv "${inifile}" "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${srun}/numangs_${nangs}/."
  cd "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}/${srun}/numangs_${nangs}"
  "${pythondir}/LBLRTM_RRTMGP_compare.py" --config_file "${inifile}" --single_stat --plot_stats --plot_profiles --band 3
  cd "${pdir}/output_rrtmgp_${rrtmgpvernum}_atm_${atmvernum}"
done
