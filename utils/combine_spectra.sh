#!/bin/bash

elem=O
for spec in spec0_i.dat spec0_f.dat spec_xas.dat
do
    python ~/mbxaspy/utils/combine_xatm_spec_shirley.py O33/$spec 1.0 O65/$spec 1.0
    mv spec_all.dat $spec
done

