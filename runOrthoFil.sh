#!bin/bash

bin_hmmer=""		# path to hmmer binaries directory
bin_augustus=""		# path to augustus binaries directory
cfg_augustus=""		# path to augustus config directory
scr_augustus=""		# path to augustus scripts directory
loc_orthofinder="" 	# path to directory containing orthofinder.py

export ORTHOFINDER_DIR=$loc_orthofinder
export AUGUSTUS_CONFIG_PATH=$cfg_augustus
PATH="$bin_hmmer:$bin_augustus:$scr_augustus:$PATH"

unset R_HOME


python OrthoFiller.py \
	-i    \ #[sequence locations file]
	-o    \ #[output directory]
	-c 16


python OrthoFiller.py --prep \
	-g    \ #[orthogroups file]
	-s    \ #[singletons file]
	-i    \ #[sequence locations file]
	-o    \ #[output directory]
	-c 16

