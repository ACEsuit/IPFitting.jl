

run(`python ./py_phonons/phonon.py bcc tik_reg.json`)
run(`python ./py_phonons/phonon.py hcp tik_reg.json`)

run(`phonopy-bandplot --legend tik_reg_bcc.yaml CASTEP_bcc.yaml`)
run(`phonopy-bandplot --legend tik_reg_hcp.yaml CASTEP_hcp.yaml`)
