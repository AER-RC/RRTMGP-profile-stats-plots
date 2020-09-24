#!/bin/bash

scriptdir='/Users/jdelamere/Projects-AER/p1770/rebuild_rrtm_nasa_chk_rayl/scripts'
srcdir='/Users/jdelamere/Projects-AER/p1770/rebuild_rrtm_nasa_chk_rayl/src'

cp TAPE5.self.mappee TAPE5
./lblrtm_v12.8_OS_X_intel_dbl
${scriptdir}/mv_files.sh mv ODint NCint

cp TAPE5.all.mapper TAPE5
./lblrtm_v12.8_OS_X_intel_dbl

${srcdir}/kdis_sort_sw_garand sw 6 M 1 10
