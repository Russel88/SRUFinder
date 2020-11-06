$PYTHON setup.py install

cat >${RECIPE_DIR}/activate.sh <<EOF
#!/bin/sh
export SRUFINDER_DB="${PREFIX}/db"
export TSS_DATA="${PREFIX}/db"
EOF

cat >${RECIPE_DIR}/deactivate.sh <<EOF
#!/bin/sh
unset SRUFINDER_DB
unset TSS_DATA
EOF

for CHANGE in "activate" "deactivate"
do
    mkdir -p "${PREFIX}/etc/conda/${CHANGE}.d"
    cp "${RECIPE_DIR}/${CHANGE}.sh" "${PREFIX}/etc/conda/${CHANGE}.d/${PKG_NAME}_${CHANGE}.sh"
done

rm ${RECIPE_DIR}/activate.sh
rm ${RECIPE_DIR}/deactivate.sh

mkdir -p ${PREFIX}/db
cp ${RECIPE_DIR}/data/* ${PREFIX}/db/

cp ${RECIPE_DIR}/bin/bprom ${PREFIX}/bin/bprom

