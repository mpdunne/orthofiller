#!bin/bash

bin_hmmer=""          #path to hmmer binaries directory
bin_augustus=""	      #path to augustus binaries directory
config_augustus=""    #path to augustus config directory

export AUGUSTUS_CONFIG_PATH=$config_augustus;
PATH="$bin_hmmer:$bin_augustus:$PATH"

python OrthoFiller.py --prep \
	-g    \ #[orthogroups file]
	-s    \ #[singletons file]
	-i    \ #[sequence locations file]
	-o    \ #[output directory]
	-c 16

#python OrthoFiller.py \
#	-i    \ #[sequence locations file]
#	-o    \ #[output directory]
#	-c 16

