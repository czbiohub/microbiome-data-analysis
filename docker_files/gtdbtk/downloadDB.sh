#!/usr/bin/env bash

echo "Downloading the GTDB-Tk database to ${GTDBTK_DATA_PATH}..."

DB_REMOTE_PATH="s3://czbiohub-microbiome/ReferenceDBs/GTDB/20191008_R04-RS89/gtdbtk_r89_data.tar.gz"

DBNAME=$(basename ${DB_REMOTE_PATH})
# GTDBTK_DB_PATH is defined in build.sh, store the db there
aws s3 cp "${DB_REMOTE_PATH}" "${GTDBTK_DATA_PATH}"
tar xvzf "${GTDBTK_DATA_PATH}/${DBNAME}" -C "${GTDBTK_DATA_PATH}" --strip 1
rm "${GTDBTK_DATA_PATH}/${DBNAME}"

echo "GTDB-Tk database has been successfully downloaded."

exit 0